from structlog import get_logger

from ref_builder.ncbi.models import NCBIGenbank
from ref_builder.otu.utils import (
    RefSeqConflictError,
    check_isolate_size,
    check_sequence_length,
    fetch_records_from_accessions,
    group_genbank_records_by_isolate,
    parse_refseq_comment,
)
from ref_builder.plan import MonopartitePlan
from ref_builder.repo import Repo
from ref_builder.resources import RepoOTU, RepoIsolate, RepoSequence
from ref_builder.utils import IsolateName


logger = get_logger("otu.isolate")


def add_genbank_isolate(
    repo: Repo,
    otu: RepoOTU,
    accessions: list[str],
    ignore_cache: bool = False,
) -> RepoIsolate | None:
    """Take a list of accessions that make up a new isolate and check that they make up
     a valid isolate before adding new isolate to the OTU.

    Download the GenBank records, categorize into an isolate bin and pass the isolate name
    and records to the add method.
    """
    otu_logger = logger.bind(taxid=otu.taxid, otu_id=str(otu.id), name=otu.name)

    records = fetch_records_from_accessions(accessions, otu.blocked_accessions, ignore_cache)

    if not records:
        return None

    record_bins = group_genbank_records_by_isolate(records)
    if len(record_bins) != 1:
        otu_logger.error("More than one isolate name found in requested accession.")
        return None

    isolate_name, isolate_records = list(record_bins.items())[0]

    try:
        check_isolate_size(plan=otu.plan, isolate_n=len(isolate_records))
    except ValueError as e:
        otu_logger.warning(
            e,
            isolate_name=str(isolate_name),
            isolate_records=sorted(isolate_records.keys()),
        )
        return None

    if type(otu.plan) is MonopartitePlan:
        if not check_sequence_length(
            records[0].sequence,
            segment_length=otu.plan.length,
            tolerance=repo.settings.default_segment_length_tolerance,
        ):
            otu_logger.error(
                "Sequence does not conform to plan length.", accession=accessions
            )
            return None

    if (isolate_id := otu.get_isolate_id_by_name(isolate_name)) is not None:
        otu_logger.debug(f"{isolate_name} already exists in this OTU")

        for record in isolate_records.values():
            if record.accession.startswith("NC_"):
                raise RefSeqConflictError(
                    f"Potential RefSeq replacement for contents of {isolate_name}",
                    isolate_id=isolate_id,
                    isolate_name=isolate_name,
                    accessions=accessions,
                )
            return None

    return create_isolate_from_records(
        repo,
        otu,
        isolate_name,
        list(isolate_records.values()),
    )


def add_unnamed_isolate(
    repo: Repo,
    otu: RepoOTU,
    accessions: list[str],
    ignore_cache: bool = False,
) -> RepoIsolate | None:
    """Take a list of accessions that make up a new isolate and add a new isolate
    with an empty name field to the OTU.

    Download the GenBank records and pass the isolate name and records to the add method.
    """
    records = fetch_records_from_accessions(accessions, otu.blocked_accessions, ignore_cache)

    if records:
        return create_isolate_from_records(
            repo,
            otu,
            isolate_name=None,
            records=records,
        )

    return None


def add_and_name_isolate(
    repo: Repo,
    otu: RepoOTU,
    accessions: list[str],
    isolate_name: IsolateName,
    ignore_cache: bool = False,
) -> RepoIsolate | None:
    """Take a list of accessions that make up a new isolate and a preferred isolate name,
     then add a new isolate to the OTU.

    Download the GenBank records and pass the isolate name and records to the add method.
    """
    otu_logger = logger.bind(taxid=otu.taxid, otu_id=str(otu.id), name=otu.name)

    records = fetch_records_from_accessions(accessions, otu.blocked_accessions, ignore_cache)

    if records:
        return create_isolate_from_records(
            repo,
            otu,
            isolate_name=isolate_name,
            records=records,
        )

    return None


def create_isolate_from_records(
    repo: Repo,
    otu: RepoOTU,
    isolate_name: IsolateName | None,
    records: list[NCBIGenbank],
) -> RepoIsolate | None:
    """Take a list of GenBank records that make up a new isolate
    and add them to the OTU.

    If a monopartite sequence is outside of recommended length bounds,
    automatically reject the isolate and return None.
    """
    isolate_logger = get_logger("otu.isolate").bind(
        isolate_name=str(isolate_name) if IsolateName is not None else None,
        otu_name=otu.name,
        taxid=otu.taxid,
    )

    if type(otu.plan) is MonopartitePlan:
        if not check_sequence_length(
            records[0].sequence,
            segment_length=otu.plan.length,
            tolerance=otu.plan.length_tolerance,
        ):
            isolate_logger.error(
                "Sequence does not conform to plan length.",
                accession=records[0].accession,
                sequence_length=len(records[0].sequence),
                plan_length=otu.plan.length,
            )

            return None

    try:
        isolate = repo.create_isolate(
            otu.id,
            None,
            isolate_name,
        )
    except ValueError as e:
        if "Isolate name already exists" in str(e):
            isolate_logger.error("OTU already contains isolate with name.")
            return None

        raise

    for record in records:
        sequence = repo.create_sequence(
            otu.id,
            accession=record.accession_version,
            definition=record.definition,
            legacy_id=None,
            segment=record.source.segment,
            sequence=record.sequence,
        )

        repo.link_sequence(otu.id, isolate.id, sequence.id)

        if record.refseq:
            refseq_status, old_accession = parse_refseq_comment(record.comment)
            repo.exclude_accession(
                otu.id,
                old_accession,
            )

    isolate_logger.info(
        f"{isolate.name} created",
        isolate_id=str(isolate.id),
        sequences=sorted([str(record.accession) for record in records]),
    )

    return isolate


def create_monopartite_isolate(
    repo: Repo,
    otu: RepoOTU,
    isolate_name: IsolateName | None,
    record: NCBIGenbank,
) -> RepoIsolate | None:
    """Take a GenBank record that makes up a new isolate
    and add it to the OTU.
    """
    isolate_logger = get_logger("otu.isolate").bind(
        isolate_name=str(isolate_name) if IsolateName is not None else None,
        otu_name=otu.name,
        taxid=otu.taxid,
    )

    if not check_sequence_length(
        record.sequence,
        segment_length=otu.plan.length,
        tolerance=repo.settings.default_segment_length_tolerance,
    ):
        isolate_logger.error("Sequence does not conform to plan length.")
        return None

    try:
        isolate = repo.create_isolate(
            otu.id,
            None,
            isolate_name,
        )
    except ValueError as e:
        if "Isolate name already exists" in str(e):
            logger.error(
                "OTU already contains isolate with name.",
                isolate_name=str(isolate_name) if IsolateName is not None else None,
                otu_name=otu.name,
                taxid=otu.taxid,
            )
            return None

        raise

    sequence = create_sequence_from_record(repo, otu, record)

    if sequence is None:
        raise ValueError("Sequence could not be created")

    repo.link_sequence(otu.id, isolate.id, sequence.id)

    if record.refseq:
        refseq_status, old_accession = parse_refseq_comment(record.comment)
        repo.exclude_accession(
            otu.id,
            old_accession,
        )

    isolate_logger.info(
        f"{isolate.name} created",
        isolate_id=str(isolate.id),
        sequences=str(sequence.accession),
    )

    return isolate


def create_sequence_from_record(
    repo: Repo,
    otu: RepoOTU,
    record: NCBIGenbank,
    segment_name: str | None = None,
) -> RepoSequence | None:
    """Take a NCBI Nucleotide record and create a new sequence."""
    sequence = repo.create_sequence(
        otu.id,
        accession=record.accession_version,
        definition=record.definition,
        legacy_id=None,
        segment=segment_name if segment_name is not None else record.source.segment,
        sequence=record.sequence,
    )

    return sequence

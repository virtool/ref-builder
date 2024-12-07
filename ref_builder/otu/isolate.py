from uuid import UUID

from structlog import get_logger

from ref_builder.ncbi.models import NCBIGenbank
from ref_builder.otu.utils import (
    RefSeqConflictError,
    assign_records_to_segments,
    check_isolate_size,
    check_sequence_length,
    fetch_records_from_accessions,
    group_genbank_records_by_isolate,
    parse_refseq_comment,
)
from ref_builder.plan import extract_segment_name_from_record
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

    records = fetch_records_from_accessions(
        accessions, otu.blocked_accessions, ignore_cache
    )

    if not records:
        return None

    record_bins = group_genbank_records_by_isolate(records)
    if len(record_bins) != 1:
        otu_logger.error("More than one isolate name found in requested accession.")
        return None

    isolate_name, isolate_records = list(record_bins.items())[0]

    try:
        check_isolate_size(plan=otu.plan, isolate_count=len(isolate_records))
    except ValueError as e:
        otu_logger.warning(
            e,
            isolate_name=str(isolate_name),
            isolate_records=sorted(isolate_records.keys()),
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

    if otu.plan.monopartite:
        return create_monopartite_isolate(repo, otu, isolate_name, records[0])

    try:
        return create_multipartite_isolate(
            repo,
            otu,
            isolate_name,
            assigned_records=assign_records_to_segments(records, otu.plan),
        )
    except ValueError:
        otu_logger.exception()

    return None


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
    otu_logger = logger.bind(taxid=otu.taxid, otu_id=str(otu.id), name=otu.name)

    records = fetch_records_from_accessions(
        accessions, otu.blocked_accessions, ignore_cache
    )

    if not records:
        return None

    try:
        check_isolate_size(plan=otu.plan, isolate_count=len(records))
    except ValueError as e:
        otu_logger.warning(
            e,
            accessions=sorted(accessions),
            isolate_name=None,
        )
        return None

    if otu.plan.monopartite:
        return create_monopartite_isolate(
            repo,
            otu,
            isolate_name=None,
            record=records[0],
        )

    try:
        return create_multipartite_isolate(
            repo,
            otu,
            isolate_name=None,
            assigned_records=assign_records_to_segments(records, otu.plan),
        )
    except ValueError as e:
        otu_logger.error(e)

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

    records = fetch_records_from_accessions(
        accessions, otu.blocked_accessions, ignore_cache
    )

    if not records:
        return None

    if otu.plan.monopartite:
        return create_monopartite_isolate(
            repo,
            otu,
            isolate_name=isolate_name,
            record=records[0],
        )

    try:
        return create_multipartite_isolate(
            repo,
            otu,
            isolate_name=isolate_name,
            assigned_records=assign_records_to_segments(records, otu.plan),
        )
    except ValueError as e:
        otu_logger.error(e)

    return None


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
        segment_length=otu.plan.segments[0].length,
        tolerance=otu.plan.segments[0].length_tolerance,
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


def create_multipartite_isolate(
    repo: Repo,
    otu: RepoOTU,
    isolate_name: IsolateName | None,
    assigned_records: dict[UUID, NCBIGenbank],
):
    """Take a dictionary of records keyed by segment name and add them to the OTU."""
    isolate_logger = get_logger("otu.isolate").bind(
        isolate_name=str(isolate_name) if IsolateName is not None else None,
        otu_name=otu.name,
        taxid=otu.taxid,
    )

    for segment_id in assigned_records:
        matching_segment = otu.plan.get_segment_by_id(segment_id)

        record = assigned_records[segment_id]
        if not check_sequence_length(
            record.sequence,
            segment_length=matching_segment.length,
            tolerance=matching_segment.length_tolerance,
        ):
            isolate_logger.warning(
                "Sequence does not conform to recommended segment length.",
                record_accession=record.accession_version,
                record_data_length=len(record.sequence),
                segment_id=str(matching_segment.id),
                segment_name=str(matching_segment.name),
                segment_length=matching_segment.length,
                segment_tolerance=matching_segment.length_tolerance,
            )
            return None

    try:
        isolate = repo.create_isolate(
            otu_id=otu.id,
            legacy_id=None,
            name=isolate_name,
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

        raise e

    for segment_id in assigned_records:
        record = assigned_records[segment_id]

        normalized_segment_name = extract_segment_name_from_record(record)

        sequence = create_sequence_from_record(
            repo,
            otu,
            record,
            str(normalized_segment_name),
        )

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
        f"Isolate created",
        id=str(isolate.id),
        accessions=list(isolate.accessions),
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

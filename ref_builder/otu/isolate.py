from uuid import UUID

from structlog import get_logger

from ref_builder.ncbi.models import NCBIGenbank
from ref_builder.otu.utils import (
    RefSeqConflictError,
    assign_records_to_segments,
    check_sequence_length,
    fetch_records_from_accessions,
    group_genbank_records_by_isolate,
    parse_refseq_comment,
)
from ref_builder.plan import PlanConformationError
from ref_builder.repo import Repo
from ref_builder.resources import RepoIsolate, RepoOTU, RepoSequence
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

    if any(record.refseq for record in records) and not all(
        record.refseq for record in records
    ):
        raise ValueError(
            "Cannot mix RefSeq and non-RefSeq sequences in multipartite isolate."
        )

    binned_records = group_genbank_records_by_isolate(records)

    if len(binned_records) != 1:
        otu_logger.error("More than one isolate name found in requested accession.")
        return None

    isolate_name, isolate_records = next(iter(binned_records.items()))

    if (isolate_id := otu.get_isolate_id_by_name(isolate_name)) is not None:
        otu_logger.warning(
            "Isolate name already exists in this OTU.", name=isolate_name
        )

        for record in isolate_records.values():
            if record.accession.startswith("NC_"):
                raise RefSeqConflictError(
                    f"Potential RefSeq replacement for contents of {isolate_name}",
                    isolate_id=isolate_id,
                    isolate_name=isolate_name,
                    accessions=accessions,
                )
            return None

    with repo.use_transaction() as transaction:
        try:
            isolate = create_isolate(
                repo,
                otu,
                isolate_name,
                records,
            )
        except ValueError:
            otu_logger.exception()

        if isolate:
            return isolate

        transaction.abort()

    return None


def add_unnamed_isolate(
    repo: Repo,
    otu: RepoOTU,
    accessions: list[str],
    ignore_cache: bool = False,
) -> RepoIsolate | None:
    """Create an unnamed isolate from a list of accessions.

    Download the GenBank records and pass the isolate name and records to the add
    method.
    """
    records = fetch_records_from_accessions(
        accessions, otu.blocked_accessions, ignore_cache
    )

    if not records:
        return None

    with repo.use_transaction() as transaction:
        isolate = create_isolate(repo, otu, None, records)

        if isolate:
            return isolate

        transaction.abort()

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
    records = fetch_records_from_accessions(
        accessions, otu.blocked_accessions, ignore_cache
    )

    if not records:
        return None

    with repo.use_transaction() as transaction:
        isolate = create_isolate(repo, otu, isolate_name, records)

        if isolate:
            return isolate

        transaction.abort()

    return None


def create_isolate(
    repo: Repo,
    otu: RepoOTU,
    isolate_name: IsolateName | None,
    records: list[NCBIGenbank],
):
    """Take a dictionary of records keyed by segment name and add them to the OTU."""
    log = get_logger(
        name=str("Unnamed" if isolate_name is None else str(isolate_name)),
        otu_name=otu.name,
        otu_id=str(otu.id),
        taxid=otu.taxid,
    )

    try:
        assigned = assign_records_to_segments(records, otu.plan)

    except PlanConformationError as e:
        log.warning(
            str(e),
            segment_names=[str(segment.name) for segment in otu.plan.segments],
        )

        return None

    for segment_id, record in assigned.items():
        matching_segment = otu.plan.get_segment_by_id(segment_id)

        if not check_sequence_length(
            record.sequence,
            matching_segment.length,
            matching_segment.length_tolerance,
        ):
            log.debug(
                "Sequence does not conform to plan length.",
                accession=record.accession_version,
                length=len(record.sequence),
                segment=(matching_segment.name, matching_segment.id),
                required_length=matching_segment.length,
                tolerance=matching_segment.length_tolerance,
            )

            return None

    try:
        isolate = repo.create_isolate(
            otu.id,
            legacy_id=None,
            name=isolate_name,
        )
    except ValueError as e:
        if "Isolate name already exists" in str(e):
            log.debug("OTU already contains isolate with name.")
            return None

        raise

    for segment_id, record in assigned.items():
        sequence = create_sequence_from_record(repo, otu, record, segment_id)

        if sequence is None:
            raise ValueError("Sequence could not be created.")

        repo.link_sequence(otu.id, isolate.id, sequence.id)

        if record.refseq:
            _, old_accession = parse_refseq_comment(record.comment)

            repo.exclude_accession(
                otu.id,
                old_accession,
            )

    log.info("Isolate created", id=str(isolate.id))

    return isolate


def rename_isolate(
    repo: Repo,
    otu: RepoOTU,
    isolate_id: UUID,
    isolate_name: IsolateName,
) -> RepoIsolate | None:
    isolate_ = otu.get_isolate(isolate_id)

    if isolate_.name == isolate_name:
        logger.warning(
            "Isolate name unchanged.", id=str(isolate_.id), name=str(isolate_.name)
        )

    else:
        old_name = isolate_.name

        with repo.use_transaction():
            new_name = repo.rename_isolate(
                otu.id,
                isolate_id=isolate_id,
                isolate_name=isolate_name,
            )

        logger.info("Isolate renamed", old_name=str(old_name), name=str(new_name))

    return repo.get_otu(otu.id).get_isolate(isolate_.id)


def create_sequence_from_record(
    repo: Repo,
    otu: RepoOTU,
    record: NCBIGenbank,
    segment_id: UUID,
) -> RepoSequence | None:
    """Take a NCBI Nucleotide record and create a new sequence."""
    return repo.create_sequence(
        otu.id,
        accession=record.accession_version,
        definition=record.definition,
        legacy_id=None,
        segment=segment_id,
        sequence=record.sequence,
    )

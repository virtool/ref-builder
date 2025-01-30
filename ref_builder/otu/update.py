from collections import defaultdict
from collections.abc import Collection, Iterable
from uuid import UUID

from structlog import get_logger

from ref_builder.ncbi.client import NCBIClient
from ref_builder.ncbi.models import NCBIGenbank
from ref_builder.otu.isolate import (
    assign_records_to_segments,
    create_isolate,
)
from ref_builder.otu.utils import (
    DeleteRationale,
    get_segments_min_length,
    get_segments_max_length,
    group_genbank_records_by_isolate,
    parse_refseq_comment,
)
from ref_builder.repo import Repo
from ref_builder.resources import RepoIsolate, RepoOTU
from ref_builder.utils import Accession, IsolateName

logger = get_logger("otu.update")

OTU_FEEDBACK_INTERVAL = 100
RECORD_FETCH_LIMIT = 25


def auto_update_otu(
    repo: Repo,
    otu: RepoOTU,
    ignore_cache: bool = False,
) -> None:
    """Fetch new accessions for the OTU and create isolates as possible."""
    ncbi = NCBIClient(False)

    log = logger.bind(taxid=otu.taxid, otu_id=str(otu.id), name=otu.name)

    accessions = ncbi.filter_accessions(ncbi.fetch_accessions_by_taxid(otu.taxid))

    fetch_list = sorted({accession.key for accession in accessions} - otu.blocked_accessions)

    if fetch_list:
        log.info("Syncing OTU with Genbank.")
        update_otu_with_accessions(repo, otu, fetch_list, ignore_cache)
    else:
        log.info("OTU is up to date.")


def batch_update_repo(
    repo: Repo,
    ignore_cache: bool = False,
):
    """Update all OTUS.

    TO BE IMPLEMENTED.
    """
    repo_logger = logger.bind(path=str(repo.path))

    repo_logger.info("Starting batch update...")

    taxid_new_accession_index = batch_fetch_new_accessions(
        repo.iter_otus(), ignore_cache
    )

    if not taxid_new_accession_index:
        logger.info("OTUs are up to date.")
        return None

    fetch_set = set()
    for taxid in taxid_new_accession_index:
        fetch_set.update(taxid_new_accession_index[taxid])

    logger.info("New accessions found.", accession_count=len(taxid_new_accession_index))

    indexed_records = batch_fetch_new_records(fetch_set, ignore_cache)


def batch_fetch_new_accessions(
    otu_iterable: Iterable[RepoOTU],
    ignore_cache: bool = False,
) -> dict[int, set[str]]:
    ncbi = NCBIClient(ignore_cache)

    otu_counter = 0

    taxid_accession_index = {}

    for otu in otu_iterable:
        otu_logger = logger.bind(taxid=otu.taxid, name=otu.name)

        otu_counter += 1

        if otu_counter % OTU_FEEDBACK_INTERVAL == 0:
            otu_logger.info(
                "Fetching accession updates...", otu_counter=otu_counter,
            )

        raw_accessions = ncbi.fetch_accessions_by_taxid(
            otu.taxid,
            sequence_min_length=get_segments_min_length(otu.plan.segments),
            sequence_max_length=get_segments_max_length(otu.plan.segments),
        )

        otu_accessions = ncbi.filter_accessions(raw_accessions)

        otu_fetch_set = {
                            accession.key for accession in otu_accessions
                        } - otu.blocked_accessions

        if otu_fetch_set:
            otu_logger.debug(
                "Potential accessions found.",
                accession_count=len(otu_fetch_set),
                otu_counter=otu_counter,
            )

            taxid_accession_index[otu.taxid] = otu_fetch_set

    return taxid_accession_index


def batch_fetch_new_records(
    accessions: Collection[str],
    ignore_cache: bool = False,
) -> dict[str, NCBIGenbank]:
    """Download a batch of records and return in a dictionary indexed by accession."""
    if not accessions:
        return {}

    ncbi = NCBIClient(ignore_cache)

    fetch_list = list(accessions)

    indexed_records = {}
    for iterator in range(len(fetch_list)):
        logger.info("Fetching records...", page=iterator, page_size=RECORD_FETCH_LIMIT)
        chunked_records = ncbi.fetch_genbank_records(
            fetch_list[
                ((iterator - 1) * RECORD_FETCH_LIMIT):(iterator * RECORD_FETCH_LIMIT)
            ],
        )

        indexed_records.update(
            {record.accession: record for record in chunked_records}
        )

    if indexed_records:
        return indexed_records

    logger.info("No valid records found.")

    return {}


def update_isolate_from_accessions(
    repo: Repo,
    otu: RepoOTU,
    isolate_name: IsolateName,
    accessions: list[str],
    ignore_cache: bool = False,
) -> RepoIsolate | None:
    """Fetch the records attached to a given list of accessions and rebuild the isolate with it."""
    if (isolate_id := otu.get_isolate_id_by_name(isolate_name)) is None:
        logger.error("OTU does not include isolate.", name=isolate_name)
        return None

    ncbi = NCBIClient(ignore_cache)

    return update_isolate_from_records(
        repo, otu, isolate_id, ncbi.fetch_genbank_records(accessions)
    )


def update_isolate_from_records(
    repo: Repo,
    otu: RepoOTU,
    isolate_id: UUID,
    records: list[NCBIGenbank],
) -> RepoIsolate | None:
    """Take a list of GenBank records and replace the existing sequences,
    adding the previous sequence accessions to the excluded accessions list.
    """
    isolate = otu.get_isolate(isolate_id)

    if len(records) == 1:
        assigned = {otu.plan.segments[0].id: records[0]}
    else:
        try:
            assigned = assign_records_to_segments(records, otu.plan)
        except ValueError as exc:
            logger.error(exc)
            return None

    for segment_id, record in assigned.items():
        _, accession = parse_refseq_comment(record.comment)

        if accession in isolate.accessions:
            old_sequence = isolate.get_sequence_by_accession(accession)

            new_sequence = repo.create_sequence(
                otu.id,
                accession=record.accession_version,
                definition=record.definition,
                legacy_id=None,
                segment=segment_id,
                sequence=record.sequence,
            )

            if new_sequence is None:
                logger.error("Isolate update failed when creating new sequence.")
                return None

            repo.replace_sequence(
                otu.id,
                isolate.id,
                new_sequence.id,
                replaced_sequence_id=old_sequence.id,
                rationale=DeleteRationale.REFSEQ,
            )

            repo.exclude_accession(otu.id, old_sequence.accession.key)

    otu = repo.get_otu(otu.id)
    isolate = otu.get_isolate(isolate_id)

    logger.info(
        "Isolate updated",
        name=str(isolate.name),
        id=str(isolate.id),
        accessions=sorted([str(accession) for accession in isolate.accessions]),
    )

    return isolate


def update_otu_with_accessions(
    repo: Repo,
    otu: RepoOTU,
    accessions: list,
    ignore_cache: bool = False,
) -> list[IsolateName]:
    """Take a list of accessions, filter for eligible accessions and
    add new sequences to the OTU.
    """
    log = logger.bind(taxid=otu.taxid, otu_id=str(otu.id), name=otu.name)

    ncbi = NCBIClient(ignore_cache)

    log.info(
        "Fetching records.",
        count=len(accessions),
        fetch_list=sorted(accessions),
    )

    records = ncbi.fetch_genbank_records(accessions)

    if promote_otu_accessions(repo, otu, ignore_cache):
        otu = repo.get_otu(otu.id)

    if records:
        return update_otu_with_records(repo, otu, records)


def update_otu_with_records(
    repo: Repo,
    otu: RepoOTU,
    records: list,
):
    new_isolate_names = []

    for divided_records in (
            [record for record in records if record.refseq],
            [record for record in records if not record.refseq],
    ):
        otu = repo.get_otu(otu.id)

        for isolate_name, isolate_records in group_genbank_records_by_isolate(
                divided_records
        ).items():
            isolate = create_isolate(
                repo, otu, isolate_name, list(isolate_records.values())
            )

            if isolate:
                new_isolate_names.append(isolate.name)

    return new_isolate_names


def promote_otu_accessions(
    repo: Repo, otu: RepoOTU, ignore_cache: bool = False
) -> set | None:
    """Fetch new accessions from NCBI Nucleotide and promote accessions
    with newly added RefSeq equivalents.
    """
    ncbi = NCBIClient(ignore_cache)

    log = logger.bind(otu_id=otu.id, taxid=otu.taxid)

    log.debug("Checking for promotable records.", accessions=sorted(otu.accessions))

    accessions = ncbi.filter_accessions(ncbi.fetch_accessions_by_taxid(otu.taxid))
    fetch_list = sorted({accession.key for accession in accessions} - otu.blocked_accessions)

    if fetch_list:
        records = ncbi.fetch_genbank_records(fetch_list)

        log.debug(
            "New accessions found. Checking for promotable records.",
            fetch_list=fetch_list,
        )

        return promote_otu_accessions_from_records(repo, otu, records)

    log.info("Records are already up to date.")


def promote_otu_accessions_from_records(
    repo: Repo, otu: RepoOTU, records: list[NCBIGenbank]
) -> set:
    """Take a list of records and check them against the contents of an OTU
    for promotable RefSeq sequences. Return a list of promoted accessions.
    """
    otu_logger = logger.bind(otu_id=str(otu.id), taxid=otu.taxid)

    refseq_records = [record for record in records if record.refseq]

    promoted_accessions = set()

    promoted_isolates = defaultdict(list)

    for record in refseq_records:
        _, predecessor_accession = parse_refseq_comment(record.comment)

        if predecessor_accession in otu.accessions:
            otu_logger.debug(
                "Replaceable accession found",
                predecessor_accession=predecessor_accession,
                promoted_accession=record.accession,
            )

            (
                isolate_id,
                sequence_id,
            ) = otu.get_sequence_id_hierarchy_from_accession(predecessor_accession)

            promoted_isolates[isolate_id].append((predecessor_accession, record))

    for isolate_id in promoted_isolates:
        records = [
            record for predecession_accession, record in promoted_isolates[isolate_id]
        ]
        isolate = otu.get_isolate(isolate_id)

        otu_logger.debug(
            f"Promoting {isolate.name}", predecessor_accessions=isolate.accessions
        )

        try:
            promoted_isolate = update_isolate_from_records(repo, otu, isolate_id, records)
        except ValueError as e:
            logger.exception(e)
            continue

        otu_logger.debug(
            f"{isolate.name} sequences promoted using RefSeq data",
            accessions=promoted_isolate.accessions,
        )

        promoted_accessions.update(promoted_isolate.accessions)

    if promoted_accessions:
        otu_logger.debug(
            "Promoted records",
            count=len(promoted_accessions),
            promoted_accessions=sorted(promoted_accessions),
        )
    else:
        otu_logger.debug("All accessions are up to date")

    return promoted_accessions


def _bin_refseq_records(
    records: list[NCBIGenbank],
) -> tuple[list[NCBIGenbank], list[NCBIGenbank]]:
    """Return a list of GenBank records as two lists, RefSeq and non-RefSeq."""
    refseq_records = []
    non_refseq_records = []

    for record in records:
        if record.refseq:
            refseq_records.append(record)
        else:
            non_refseq_records.append(record)

    if len(refseq_records) + len(non_refseq_records) != len(records):
        raise ValueError("Invalid total number of records")

    return refseq_records, non_refseq_records

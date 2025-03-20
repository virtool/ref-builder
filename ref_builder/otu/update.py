import datetime
from collections import defaultdict
from collections.abc import Collection, Iterable, Iterator
from pathlib import Path
from uuid import UUID

import arrow
from pydantic import RootModel
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
from ref_builder.utils import IsolateName
from ref_builder.transaction import AbortTransactionError

logger = get_logger("otu.update")

OTU_FEEDBACK_INTERVAL = 100
"""A default interval for batch OTU feedback."""

RECORD_FETCH_CHUNK_SIZE = 500
"""A default chunk size for NCBI EFetch calls."""

UPDATE_COOLDOWN_INTERVAL_IN_DAYS = 14


BatchFetchIndex = RootModel[dict[int, set[str]]]
"""Assists in reading and writing the fetched accessions by taxid index from file."""


def auto_update_otu(
    repo: Repo,
    otu: RepoOTU,
    fetch_index_path: Path | None = None,
    start_date: datetime.date | None = None,
    ignore_cache: bool = False,
) -> RepoOTU:
    """Fetch new accessions for the OTU and create isolates as possible."""
    ncbi = NCBIClient(False)

    log = logger.bind(taxid=otu.taxid, otu_id=str(otu.id), name=otu.name)

    fetch_set = set()

    if isinstance(fetch_index_path, Path):
        log.info("Loading fetch index...", fetch_index_path=str(fetch_index_path))

        fetch_index = _load_fetch_index(fetch_index_path)

        fetch_set = fetch_index.get(otu.taxid, fetch_set)

    if not fetch_set:
        accessions = ncbi.filter_accessions(
            ncbi.fetch_accessions_by_taxid(
                otu.taxid,
                sequence_min_length=get_segments_min_length(otu.plan.segments),
                sequence_max_length=get_segments_max_length(otu.plan.segments),
                modification_date_start=start_date,
            ),
        )

        fetch_set = {accession.key for accession in accessions} - otu.blocked_accessions

    if fetch_set:
        log.info("Syncing OTU with Genbank.")
        new_isolate_ids = update_otu_with_accessions(repo, otu, fetch_set, ignore_cache)

        if new_isolate_ids:
            log.info("Added new isolates", isolate_ids=new_isolate_ids)

    else:
        log.info("OTU is up to date.")

    repo.write_otu_update_history_entry(otu.id)

    return repo.get_otu(otu.id)


def batch_update_repo(
    repo: Repo,
    start_date: datetime.date | None = None,
    chunk_size: int = RECORD_FETCH_CHUNK_SIZE,
    fetch_index_path: Path | None = None,
    precache_records: bool = False,
    skip_recently_updated: bool = True,
    ignore_cache: bool = False,
) -> set[UUID]:
    """Fetch new accessions for all OTUs in the repo and create isolates as possible."""
    operation_run_timestamp = arrow.utcnow().naive

    updated_otu_ids = set()

    repo_logger = logger.bind(
        path=str(repo.path),
        precache_records=precache_records,
    )
    if start_date is not None:
        repo_logger = repo_logger.bind(start_date.isoformat())

    repo_logger.info("Starting batch update...")

    if fetch_index_path is None:
        if skip_recently_updated:
            otu_iterator = (
                otu
                for otu in repo.iter_otus()
                if _otu_is_cooled(
                    repo,
                    otu.id,
                    timestamp_current=operation_run_timestamp,
                )
            )
        else:
            otu_iterator = repo.iter_otus()

        taxid_new_accession_index = batch_fetch_new_accessions(
            otu_iterator,
            modification_date_start=start_date,
            ignore_cache=ignore_cache,
        )

        fetch_index_cache_path = _cache_fetch_index(
            taxid_new_accession_index, repo.path / ".cache"
        )

        repo_logger.info("Fetch index cached", fetch_index_path=fetch_index_cache_path)

    else:
        repo_logger.info(
            "Loading fetch index...", fetch_index_path=str(fetch_index_path)
        )

        taxid_new_accession_index = _load_fetch_index(fetch_index_path)

    if not taxid_new_accession_index:
        logger.info("OTUs are up to date.")

        return updated_otu_ids

    fetch_set = {
        accession
        for otu_accessions in taxid_new_accession_index.values()
        for accession in otu_accessions
    }

    logger.info(
        "Batch fetch index contains potential new accessions.",
        otu_count=len(taxid_new_accession_index),
    )

    if precache_records:
        logger.info("Precaching records...", accession_count=len(fetch_set))

        indexed_records = batch_fetch_new_records(
            fetch_set,
            chunk_size=chunk_size,
            ignore_cache=ignore_cache,
        )

        if not indexed_records:
            logger.info("No valid accessions found.")
            return updated_otu_ids

        logger.info(
            "Checking new records against OTUs.",
            otu_count=len(taxid_new_accession_index),
        )

        for taxid, accessions in taxid_new_accession_index.items():
            if (otu_id := repo.get_otu_id_by_taxid(taxid)) is None:
                logger.debug("No corresponding OTU found in this repo", taxid=taxid)
                continue

            if skip_recently_updated and not _otu_is_cooled(
                repo,
                otu_id,
                timestamp_current=operation_run_timestamp,
            ):
                logger.info(
                    "This OTU was updated recently. Skipping...",
                    cooldown=UPDATE_COOLDOWN_INTERVAL_IN_DAYS,
                    otu_id=str(otu_id),
                    taxid=str(taxid),
                )
                continue

            otu_records = [
                record
                for accession in accessions
                if (record := indexed_records.get(accession)) is not None
            ]
            if otu_records:
                isolate_ids = _process_records_into_otu(
                    repo, repo.get_otu(otu_id), otu_records
                )

                if isolate_ids:
                    updated_otu_ids.add(otu_id)

            repo.write_otu_update_history_entry(otu_id)

    else:
        ncbi = NCBIClient(ignore_cache)

        for taxid, accessions in taxid_new_accession_index.items():
            if (otu_id := repo.get_otu_id_by_taxid(taxid)) is None:
                logger.debug("No corresponding OTU found in this repo", taxid=taxid)
                continue

            if skip_recently_updated and not _otu_is_cooled(
                repo,
                otu_id,
                timestamp_current=operation_run_timestamp,
            ):
                logger.info(
                    "This OTU was updated recently. Skipping...",
                    cooldown=UPDATE_COOLDOWN_INTERVAL_IN_DAYS,
                    otu_id=str(otu_id),
                    taxid=str(taxid),
                )
                continue

            otu_records = ncbi.fetch_genbank_records(accessions)

            if otu_records:
                isolate_ids = _process_records_into_otu(
                    repo, repo.get_otu(otu_id), otu_records
                )

                if isolate_ids:
                    updated_otu_ids.add(otu_id)

            repo.write_otu_update_history_entry(otu_id)

    repo_logger.info("Batch update complete.")

    return updated_otu_ids


def batch_fetch_new_accessions(
    otus: Iterable[RepoOTU],
    modification_date_start: datetime.date | None = None,
    modification_date_end: datetime.date | None = None,
    ignore_cache: bool = False,
) -> dict[int, set[str]]:
    """Check OTU iterator for new accessions and return results indexed by taxid."""
    ncbi = NCBIClient(ignore_cache)

    otu_counter = 0

    taxid_accession_index = {}

    for otu in otus:
        otu_logger = logger.bind(taxid=otu.taxid, name=otu.name)

        otu_counter += 1

        if otu_counter % OTU_FEEDBACK_INTERVAL == 0:
            otu_logger.info(
                "Fetching accession updates...",
                otu_counter=otu_counter,
            )

        raw_accessions = ncbi.fetch_accessions_by_taxid(
            otu.taxid,
            sequence_min_length=get_segments_min_length(otu.plan.segments),
            sequence_max_length=get_segments_max_length(otu.plan.segments),
            modification_date_start=modification_date_start,
            modification_date_end=modification_date_end,
        )

        accessions_to_fetch = {
            accession.key for accession in ncbi.filter_accessions(raw_accessions)
        } - otu.blocked_accessions

        if accessions_to_fetch:
            otu_logger.debug(
                "Potential accessions found.",
                accession_count=len(accessions_to_fetch),
                otu_counter=otu_counter,
            )

            taxid_accession_index[otu.taxid] = accessions_to_fetch

    return taxid_accession_index


def batch_fetch_new_records(
    accessions: Collection[str],
    chunk_size: int = RECORD_FETCH_CHUNK_SIZE,
    ignore_cache: bool = False,
) -> dict[str, NCBIGenbank]:
    """Download a batch of records and return in a dictionary indexed by accession."""
    fetch_logger = logger.bind(
        accession_count=len(accessions),
        chunk_size=chunk_size,
        ignore_cache=ignore_cache,
    )

    if not accessions:
        return {}

    ncbi = NCBIClient(ignore_cache)

    fetch_list = list(accessions)

    page_counter = 0

    indexed_records = {}
    for fetch_list_chunk in iter_fetch_list(fetch_list, chunk_size):
        fetch_logger.info("Fetching records...", page_counter=page_counter)

        chunked_records = ncbi.fetch_genbank_records(fetch_list_chunk)

        indexed_records.update({record.accession: record for record in chunked_records})

        page_counter += 1

    if indexed_records:
        return indexed_records

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
        except ValueError as e:
            logger.debug(e)
            return None

    for segment_id, record in assigned.items():
        _, accession = parse_refseq_comment(record.comment)

        if accession in isolate.accessions:
            old_sequence = isolate.get_sequence_by_accession(accession)

            # Use one transaction per sequence
            with repo.use_transaction():
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


def _process_records_into_otu(
    repo: Repo,
    otu: RepoOTU,
    records: list[NCBIGenbank],
):
    genbank_records, refseq_records = [], []

    for record in records:
        if record.refseq:
            refseq_records.append(record)

        else:
            genbank_records.append(record)

    if promote_otu_accessions_from_records(
        repo,
        otu=repo.get_otu(otu.id),
        records=refseq_records,
    ):
        otu = repo.get_otu(otu.id)

    new_isolate_ids = update_otu_with_records(
        repo,
        otu=otu,
        records=genbank_records,
    )

    repo.get_otu(otu.id)

    return new_isolate_ids


def update_otu_with_accessions(
    repo: Repo,
    otu: RepoOTU,
    accessions: Collection[str],
    ignore_cache: bool = False,
) -> list[UUID]:
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

    if records:
        return _process_records_into_otu(repo, otu, records)


def update_otu_with_records(
    repo: Repo,
    otu: RepoOTU,
    records: list[NCBIGenbank],
) -> list[UUID]:
    """Take a list of downloaded NCBI Genbank records, filter for eligible records
    and add new sequences to the OTU.
    """
    new_isolate_ids = []

    for divided_records in (
        [record for record in records if record.refseq],
        [record for record in records if not record.refseq],
    ):
        otu = repo.get_otu(otu.id)

        for isolate_name, isolate_records in group_genbank_records_by_isolate(
            divided_records
        ).items():
            with repo.use_transaction():
                isolate = create_isolate(
                    repo, otu, isolate_name, list(isolate_records.values())
                )

                if isolate is None:
                    raise AbortTransactionError()

            if isolate:
                new_isolate_ids.append(isolate.id)

    return new_isolate_ids


def promote_otu_accessions(
    repo: Repo, otu: RepoOTU, ignore_cache: bool = False
) -> set | None:
    """Fetch new accessions from NCBI Nucleotide and promote accessions
    with newly added RefSeq equivalents.
    """
    ncbi = NCBIClient(ignore_cache)

    log = logger.bind(otu_id=otu.id, taxid=otu.taxid)

    log.info("Checking for promotable sequences.")

    accessions = ncbi.filter_accessions(
        ncbi.fetch_accessions_by_taxid(
            otu.taxid,
            sequence_min_length=get_segments_min_length(otu.plan.segments),
            sequence_max_length=get_segments_max_length(otu.plan.segments),
            refseq_only=True,
        ),
    )
    fetch_set = {accession.key for accession in accessions} - otu.blocked_accessions

    if fetch_set:
        records = ncbi.fetch_genbank_records(fetch_set)

        log.debug(
            "New accessions found. Checking for promotable records.",
            fetch_list=sorted(fetch_set),
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

    initial_exceptions = otu.excluded_accessions.copy()

    refseq_records = [record for record in records if record.refseq]

    promoted_accessions = set()

    promoted_isolates = defaultdict(list)

    for record in refseq_records:
        try:
            _, predecessor_accession = parse_refseq_comment(record.comment)

        except ValueError as e:
            logger.debug(e, accession=record.accession_version, comment=record.comment)
            continue

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
            promoted_isolate = update_isolate_from_records(
                repo, otu, isolate_id, records
            )
        except ValueError as e:
            logger.debug(e)
            continue

        if promoted_isolate is not None:
            otu_logger.info(
                f"Isolate promoted using RefSeq data",
                isolate_id=str(isolate.id),
                isolate_name=str(isolate.name) if isolate.name is not None else None,
                accessions=promoted_isolate.accessions,
            )

            promoted_accessions.update(promoted_isolate.accessions)

    if promoted_accessions:
        otu = repo.get_otu(otu.id)

        otu_logger.info(
            "Promoted records",
            count=len(promoted_accessions),
            promoted_accessions=sorted(promoted_accessions),
            new_excluded_accessions=sorted(
                otu.excluded_accessions - initial_exceptions
            ),
        )
    else:
        otu_logger.debug("All accessions are up to date")

    return promoted_accessions


def iter_fetch_list(
    fetch_list: list[str], page_size=RECORD_FETCH_CHUNK_SIZE
) -> Iterator[list[str]]:
    """Divide a list of accessions and yield in pages."""
    page_size = max(1, page_size)
    page_count = len(fetch_list) // page_size + int(len(fetch_list) % page_size != 0)

    for iterator in range(page_count):
        yield fetch_list[iterator * page_size : (iterator + 1) * page_size]


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


def _generate_datestamp_filename():
    """Get the current UTC date and return as a a filename_safe string."""
    timestamp = arrow.utcnow().naive

    return f"{timestamp:%Y}_{timestamp:%m}_{timestamp:%d}"


def _cache_fetch_index(
    fetch_index: dict[int, set[str]],
    cache_path: Path,
) -> Path | None:
    validated_fetch_index = BatchFetchIndex.model_validate(fetch_index)

    fetch_index_path = (
        cache_path / f"fetch_index__{_generate_datestamp_filename()}.json"
    )

    with open(fetch_index_path, "w") as f:
        f.write(validated_fetch_index.model_dump_json())

    if fetch_index_path.exists():
        return fetch_index_path


def _load_fetch_index(path: Path) -> dict[int, set[str]] | None:
    if not path.exists():
        return None

    if path.suffix != ".json":
        return None

    with open(path, "rb") as f:
        fetch_index = BatchFetchIndex.model_validate_json(f.read())

    if fetch_index:
        return fetch_index.model_dump()


def _is_past_cooldown(
    otu_timestamp: datetime.datetime | None,
    current_timestamp: datetime.datetime | None,
    cooldown: int = UPDATE_COOLDOWN_INTERVAL_IN_DAYS,
) -> bool:
    """Return True if ``otu_timestamp`` is None or more than ``cooldown`` days past
    ``current_timestamp``, else False."""
    if otu_timestamp is None:
        return True

    if current_timestamp is None:
        current_timestamp = arrow.utcnow().naive

    return current_timestamp - otu_timestamp > datetime.timedelta(days=cooldown)


def _otu_is_cooled(
    repo: Repo,
    otu_id: UUID,
    timestamp_current: datetime.datetime | None,
    cooldown: int = UPDATE_COOLDOWN_INTERVAL_IN_DAYS,
) -> bool:
    """Return the update cooldown status of an OTU.

    Return True if:
        A) The OTU has been recorded as updated within ``cooldown`` days,
        B) The OTU was created within ``cooldown`` days.
        C) The OTU was created outside and last updated outside ``cooldown`` days.
    """
    cooldown_delta = datetime.timedelta(days=cooldown)

    if timestamp_current is None:
        timestamp_current = arrow.utcnow().naive

    if (timestamp_last_updated := repo.get_otu_last_updated(otu_id)) is not None:
        return timestamp_current - timestamp_last_updated > cooldown_delta

    timestamp_created = repo.get_otu_first_created(otu_id)

    if (timestamp_current - timestamp_created) <= cooldown_delta:
        return True

    timestamp_latest = repo.get_otu_last_modified(otu_id)

    if (timestamp_current - timestamp_latest) > cooldown_delta:
        return True

    logger.debug(
        "Last modified timestamp occured within cooldown",
        timestamp_current=timestamp_current.isoformat(),
        timestamp_otu_created=timestamp_created.isoformat(),
        timestamp_otu_latest=timestamp_latest.isoformat(),
        otu_last_modified_delta=timestamp_current - timestamp_latest,
    )

    return False

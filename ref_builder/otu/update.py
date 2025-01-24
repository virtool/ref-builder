from collections import defaultdict
from uuid import UUID

from structlog import get_logger

from ref_builder.ncbi.client import NCBIClient
from ref_builder.ncbi.models import NCBIGenbank
from ref_builder.otu.isolate import (
    assign_records_to_segments,
    create_monopartite_isolate,
    create_multipartite_isolate,
)
from ref_builder.otu.utils import (
    DeleteRationale,
    check_isolate_size,
    group_genbank_records_by_isolate,
    parse_refseq_comment,
)
from ref_builder.repo import Repo
from ref_builder.resources import RepoIsolate, RepoOTU
from ref_builder.utils import Accession, IsolateName

logger = get_logger("otu.update")


def auto_update_otu(repo: Repo, otu: RepoOTU) -> None:
    """Fetch a full list of Nucleotide accessions associated with the OTU
    and pass the list to the add method.
    """
    ncbi = NCBIClient(False)

    otu_logger = logger.bind(taxid=otu.taxid, otu_id=str(otu.id), name=otu.name)

    accessions = ncbi.fetch_accessions_by_taxid(otu.taxid)

    fetch_list = sorted(
        {Accession.from_string(accession).key for accession in accessions}
        - otu.blocked_accessions
    )

    if fetch_list:
        otu_logger.info("Syncing OTU with Genbank.")
        update_otu_with_accessions(repo, otu, fetch_list)
    else:
        otu_logger.info("OTU is up to date.")


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

    new_records = ncbi.fetch_genbank_records(accessions)

    return update_isolate_from_records(repo, otu, isolate_id, new_records)


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
        assigned = assign_records_to_segments(records, otu.plan)

    for segment_id, record in assigned.items():
        status, accession = parse_refseq_comment(record.comment)

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
) -> None:
    """Take a list of accessions, filter for eligible accessions and
    add new sequences to the OTU.
    """
    otu_logger = logger.bind(taxid=otu.taxid, otu_id=str(otu.id), name=otu.name)

    ncbi = NCBIClient(ignore_cache)

    otu_logger.info(
        "Fetching records",
        count=len(accessions),
        fetch_list=sorted(accessions),
    )

    records = ncbi.fetch_genbank_records(accessions)

    promoted_accessions = promote_otu_accessions(repo, otu, ignore_cache)
    if promoted_accessions:
        otu = repo.get_otu(otu.id)

    refseq_records, non_refseq_records = _bin_refseq_records(records)

    # Add remaining RefSeq records first
    otu = repo.get_otu(otu.id)
    file_records_into_otu(repo, otu, refseq_records)

    # Add non-RefSeq records last.
    otu = repo.get_otu(otu.id)
    file_records_into_otu(repo, otu, non_refseq_records)


def file_records_into_otu(
    repo: Repo,
    otu: RepoOTU,
    records: list[NCBIGenbank],
) -> list[IsolateName]:
    """Create isolates from a list of GenBank records.

    If an isolate candidate does not match the plan, the records are not added.
    """
    otu_logger = logger.bind(taxid=otu.taxid, otu_id=str(otu.id), name=otu.name)

    record_bins = group_genbank_records_by_isolate(records)

    otu_logger.debug(
        "Found potential new isolates",
        count=len(record_bins),
        names=[str(isolate_name) for isolate_name in record_bins],
    )

    new_isolate_names = []

    for isolate_name in record_bins:
        isolate_records = record_bins[isolate_name]

        try:
            check_isolate_size(plan=otu.plan, isolate_count=len(isolate_records))
        except ValueError as e:
            otu_logger.warning(
                e,
                isolate_name=str(isolate_name),
                isolate_records=[str(accession) for accession in isolate_records],
            )
            continue

        if otu.plan.monopartite:
            isolate = create_monopartite_isolate(
                repo,
                otu,
                isolate_name,
                record=list(isolate_records.values())[0],
            )

            if isolate is not None:
                new_isolate_names.append(isolate_name)

        else:
            assigned_records = None
            try:
                assigned_records = assign_records_to_segments(
                    list(isolate_records.values()), otu.plan
                )

            except ValueError as e:
                otu_logger.error(e)

            if assigned_records:
                isolate = create_multipartite_isolate(
                    repo,
                    otu,
                    isolate_name,
                    assigned_records,
                )

                new_isolate_names.append(isolate_name)

    if new_isolate_names:
        otu_logger.info(
            "New isolates added",
            new_isolates=[str(isolate_name) for isolate_name in new_isolate_names],
        )

        return new_isolate_names

    otu_logger.info("No new isolates added.")
    return []


def promote_otu_accessions(
    repo: Repo, otu: RepoOTU, ignore_cache: bool = False
) -> set | None:
    """Fetch new accessions from NCBI Nucleotide and promote accessions
    with newly added RefSeq equivalents.
    """
    ncbi = NCBIClient(ignore_cache)

    otu_logger = logger.bind(otu_id=otu.id, taxid=otu.taxid)

    otu_logger.debug(
        "Checking for promotable records", accessions=sorted(otu.accessions)
    )

    accessions = ncbi.fetch_accessions_by_taxid(otu.taxid)
    fetch_list = sorted(
        set(Accession.from_string(accession).key for accession in accessions)
        - otu.blocked_accessions
    )

    if fetch_list:
        records = ncbi.fetch_genbank_records(fetch_list)

        otu_logger.debug(
            "New accessions found. Checking for promotable records.",
            fetch_list=fetch_list,
        )

        return promote_otu_accessions_from_records(repo, otu, records)

    otu_logger.info("Records are already up to date.")


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
        status, predecessor_accession = parse_refseq_comment(record.comment)

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

        promoted_isolate = update_isolate_from_records(repo, otu, isolate_id, records)

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

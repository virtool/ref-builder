from collections import defaultdict
from uuid import UUID

from structlog import get_logger

from ref_builder.ncbi.client import NCBIClient
from ref_builder.ncbi.models import NCBIGenbank
from ref_builder.otu.utils import (
    DeleteRationale,
    RefSeqConflictError,
    group_genbank_records_by_isolate,
    parse_refseq_comment,
)
from ref_builder.plan import MonopartitePlan, MultipartitePlan
from ref_builder.repo import Repo
from ref_builder.resources import RepoIsolate, RepoOTU, RepoSequence
from ref_builder.utils import IsolateName

logger = get_logger("otu.update")


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
    ncbi = NCBIClient(ignore_cache)

    otu_logger = logger.bind(taxid=otu.taxid, otu_id=str(otu.id), name=otu.name)

    try:
        fetch_list = _create_fetch_list(accessions, otu.blocked_accessions)
    except ValueError:
        otu_logger.error(
            "Could not create a new isolate using the requested accessions.",
            requested_accessions=accessions,
            otu_accessions=otu.accessions,
            excluded_accessions=otu.excluded_accessions,
        )
        return None

    otu_logger.info(
        "Fetching accessions",
        count={len(fetch_list)},
        fetch_list=fetch_list,
    )

    records = ncbi.fetch_genbank_records(fetch_list)

    record_bins = group_genbank_records_by_isolate(records)
    if len(record_bins) != 1:
        otu_logger.error("More than one isolate name found in requested accession.")
        return None

    isolate_name, isolate_records = list(record_bins.items())[0]

    try:
        _check_isolate_size(plan=otu.plan, isolate_n=len(isolate_records))
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
    ncbi = NCBIClient(ignore_cache)

    otu_logger = logger.bind(taxid=otu.taxid, otu_id=str(otu.id), name=otu.name)

    try:
        fetch_list = _create_fetch_list(accessions, otu.blocked_accessions)
    except ValueError:
        otu_logger.error(
            "Could not create a new isolate using the requested accessions.",
            requested_accessions=accessions,
            otu_accessions=otu.accessions,
            excluded_accessions=otu.excluded_accessions,
        )
        return None

    otu_logger.info(
        "Fetching accessions",
        count={len(fetch_list)},
        fetch_list=fetch_list,
    )

    records = ncbi.fetch_genbank_records(fetch_list)

    return create_isolate_from_records(
        repo,
        otu,
        isolate_name=None,
        records=records,
    )


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
    ncbi = NCBIClient(ignore_cache)

    otu_logger = logger.bind(taxid=otu.taxid, otu_id=str(otu.id), name=otu.name)

    try:
        fetch_list = _create_fetch_list(accessions, otu.blocked_accessions)
    except ValueError:
        otu_logger.error(
            "Could not create a new isolate using the requested accessions.",
            requested_accessions=accessions,
            otu_accessions=otu.accessions,
            excluded_accessions=otu.excluded_accessions,
        )
        return None

    otu_logger.info(
        "Fetching accessions",
        count={len(fetch_list)},
        fetch_list=fetch_list,
    )

    records = ncbi.fetch_genbank_records(fetch_list)

    return create_isolate_from_records(
        repo,
        otu,
        isolate_name=isolate_name,
        records=records,
    )


def create_isolate_from_records(
    repo: Repo,
    otu: RepoOTU,
    isolate_name: IsolateName | None,
    records: list[NCBIGenbank],
) -> RepoIsolate | None:
    """Take a list of GenBank records that make up a new isolate
    and add them to the OTU.
    """
    isolate_logger = get_logger("otu.isolate").bind(
        isolate_name=str(isolate_name) if IsolateName is not None else None,
        otu_name=otu.name,
        taxid=otu.taxid,
    )

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
        repo.create_sequence(
            otu.id,
            isolate.id,
            accession=record.accession_version,
            definition=record.definition,
            legacy_id=None,
            segment=record.source.segment,
            sequence=record.sequence,
        )

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


def set_representative_isolate(
    repo: Repo,
    otu: RepoOTU,
    isolate_id: UUID,
) -> UUID | None:
    """Sets an OTU's representative isolate to a given existing isolate ID.

    Returns the isolate ID if successful, else None.
    """
    otu_logger = logger.bind(name=otu.name, taxid=otu.taxid)

    new_representative_isolate = otu.get_isolate(isolate_id)
    if new_representative_isolate is None:
        otu_logger.error("Isolate not found. Please make a new isolate.")
        return None

    if otu.repr_isolate is not None:
        if otu.repr_isolate == new_representative_isolate.id:
            otu_logger.warning("This isolate is already the representative isolate.")
            return otu.repr_isolate

        otu_logger.warning(
            f"Replacing representative isolate {otu.repr_isolate}",
            representative_isolate_id=str(otu.repr_isolate),
        )

    repo.set_repr_isolate(otu.id, new_representative_isolate.id)

    otu_logger.info(
        f"Representative isolate set to {new_representative_isolate.name}.",
        representative_isolate_id=str(new_representative_isolate.id),
    )

    return new_representative_isolate.id


def update_isolate_from_accessions(
    repo: Repo,
    otu: RepoOTU,
    isolate_name: IsolateName,
    accessions: list[str],
    ignore_cache: bool = False,
) -> RepoIsolate | None:
    """Fetch the records attached to a given list of accessions and rebuild the isolate with it."""
    if (isolate_id := otu.get_isolate_id_by_name(isolate_name)) is None:
        logger.error(f"OTU does not include {isolate_name}.")
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

    sequences = isolate.sequences

    isolate_logger = get_logger("otu.isolate").bind(
        isolate_name=str(isolate.name),
        otu_name=otu.name,
        taxid=otu.taxid,
    )

    for record in records:
        status, accession = parse_refseq_comment(record.comment)
        if accession in isolate.accessions:
            old_sequence = isolate.get_sequence_by_accession(accession)

            new_sequence = repo.replace_sequence(
                otu.id,
                isolate.id,
                accession=record.accession_version,
                definition=record.definition,
                legacy_id=None,
                segment=record.source.segment,
                sequence=record.sequence,
                replaced_sequence_id=old_sequence.id,
                rationale=DeleteRationale.REFSEQ,
            )
            if new_sequence is None:
                isolate_logger.error("Isolate could not be refilled.")
                return None

            logger.debug(
                f"{old_sequence.accession} replaced by {new_sequence.accession}"
            )

            repo.exclude_accession(otu.id, old_sequence.accession.key)

            logger.debug(f"{accession} added to exclusion list.")

    otu = repo.get_otu(otu.id)
    isolate = otu.get_isolate(isolate_id)

    isolate_logger.info(
        f"{sorted(isolate.accessions)} added to {isolate.name}, "
        + f"replacing {sorted([str(sequence.accession) for sequence in sequences])}",
        isolate_id=str(isolate.id),
        sequences=sorted([str(record.accession) for record in records]),
    )

    logger.debug(
        "Excluded accessions updated",
        excluded_accessions=sorted(otu.excluded_accessions),
    )

    return isolate


def auto_update_otu(
    repo: Repo,
    otu: RepoOTU,
    ignore_cache: bool = False,
) -> None:
    """Fetch a full list of Nucleotide accessions associated with the OTU
    and pass the list to the add method.
    """
    ncbi = NCBIClient(ignore_cache)

    otu_logger = logger.bind(taxid=otu.taxid, otu_id=str(otu.id), name=otu.name)

    accessions = ncbi.fetch_accessions_by_taxid(otu.taxid)

    if fetch_list := list(set(accessions) - otu.blocked_accessions):
        update_otu_with_accessions(repo, otu, fetch_list)
    else:
        otu_logger.info("OTU is up to date.")


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
        "Fetching accessions",
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

    # Add add non-RefSeq records last
    otu = repo.get_otu(otu.id)
    file_records_into_otu(repo, otu, non_refseq_records)


def file_records_into_otu(
    repo: Repo,
    otu: RepoOTU,
    records: list[NCBIGenbank],
) -> list[IsolateName]:
    """Take a list of GenBank records from NCBI Nucleotide and attempt to create new isolates.
    If an isolate candidate does not match the schema, the constituent records are not added
    to the OTU.
    """
    otu_logger = logger.bind(taxid=otu.taxid, otu_id=str(otu.id), name=otu.name)

    record_bins = group_genbank_records_by_isolate(records)

    otu_logger.debug(
        f"Found {len(record_bins)} potential new isolates",
        potential_isolates=[str(isolate_name) for isolate_name in record_bins.keys()],
    )

    new_isolate_names = []

    for isolate_name in record_bins:
        isolate_records = record_bins[isolate_name]

        try:
            _check_isolate_size(plan=otu.plan, isolate_n=len(isolate_records))
        except ValueError as e:
            otu_logger.warning(
                e,
                isolate_name=str(isolate_name),
                isolate_records=[str(accession) for accession in isolate_records],
            )
            continue

        isolate = create_isolate_from_records(
            repo,
            otu,
            isolate_name,
            list(isolate_records.values()),
        )
        if isolate is not None:
            new_isolate_names.append(isolate_name)

    if new_isolate_names:
        otu_logger.info(
            "New isolates added",
            new_isolates=[str(isolate_name) for isolate_name in new_isolate_names],
        )

        return new_isolate_names

    otu_logger.info("No new isolates added.")
    return []


def exclude_accessions_from_otu(
    repo: Repo,
    otu: RepoOTU,
    accessions: list[str],
) -> None:
    """Take a list of accessions and add them to an OTU's excluded accessions list."""
    excluded_accessions = set()

    # Using a set avoid duplicate accessions.
    for accession in set(accessions):
        excluded_accessions = repo.exclude_accession(otu.id, accession)

    logger.debug(
        f"Accessions currently excluded from fetches: {excluded_accessions}",
        taxid=otu.taxid,
        otu_id=str(otu.id),
        name=otu.name,
    )


def delete_isolate_from_otu(repo: Repo, otu: RepoOTU, isolate_id: UUID) -> None:
    """Remove an isolate from a specified OTU."""
    if (isolate := otu.get_isolate(isolate_id)) is None:
        logger.error("This isolate does not exist in this OTU.")
        return

    repo.delete_isolate(otu.id, isolate.id, rationale=DeleteRationale.USER)

    logger.info(f"{isolate.name} removed.")


def replace_sequence_in_otu(
    repo: Repo,
    otu: RepoOTU,
    new_accession: str,
    replaced_accession: str,
    ignore_cache: bool = False,
) -> RepoSequence | None:
    """Replace a sequence in an OTU."""
    ncbi = NCBIClient(ignore_cache)

    (
        isolate_id,
        sequence_id,
    ) = otu.get_sequence_id_hierarchy_from_accession(replaced_accession)
    if sequence_id is None:
        logger.error("This sequence does not exist in this OTU.")
        return None

    isolate = otu.get_isolate(isolate_id)

    records = ncbi.fetch_genbank_records([new_accession])
    record = records[0]

    new_sequence = repo.replace_sequence(
        otu.id,
        isolate.id,
        accession=record.accession_version,
        definition=record.definition,
        legacy_id=None,
        segment=record.source.segment,
        sequence=record.sequence,
        replaced_sequence_id=sequence_id,
        rationale="Requested by user",
    )

    if new_sequence is not None:
        logger.info(
            f"{replaced_accession} replaced by {new_sequence.accession}.",
            new_sequence_id=new_sequence.id,
        )
        return new_sequence

    logger.error(f"{replaced_accession} could not be replaced.")


def promote_otu_accessions(
    repo: Repo, otu: RepoOTU, ignore_cache: bool = False
) -> set | None:
    """Fetch new accessions from NCBI Nucleotide and promote accessions
    with newly added RefSeq equivalents.
    """
    ncbi = NCBIClient(ignore_cache)

    accessions = ncbi.fetch_accessions_by_taxid(otu.taxid)
    fetch_list = set(accessions) - otu.blocked_accessions

    records = ncbi.fetch_genbank_records(list(fetch_list))

    return promote_otu_accessions_from_records(repo, otu, records)


def promote_otu_accessions_from_records(
    repo: Repo, otu: RepoOTU, records: list[NCBIGenbank]
) -> set:
    """Take a list of records and check them against the contents of an OTU
    for promotable RefSeq sequences. Return a list of promoted accessions.
    """
    refseq_records = [record for record in records if record.refseq]

    promoted_accessions = set()

    promoted_isolates = defaultdict(list)

    for record in refseq_records:
        status, predecessor_accession = parse_refseq_comment(record.comment)

        if predecessor_accession in otu.accessions:
            logger.debug(
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

        logger.debug(
            f"Promoting {isolate.name}", predecessor_accessions=isolate.accessions
        )

        promoted_isolate = update_isolate_from_records(repo, otu, isolate_id, records)

        logger.debug(
            f"{isolate.name} sequences promoted using RefSeq data",
            accessions=promoted_isolate.accessions,
        )

        promoted_accessions.update(promoted_isolate.accessions)

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


def _create_fetch_list(
    requested_accessions: list | set,
    blocked_accessions: set,
) -> list[str]:
    """Return the difference between a set of requested accessions
    and a set of blocked accessions.
    Raise a ValueError if the intersection is complete.
    """
    fetch_list = list(set(requested_accessions).difference(blocked_accessions))
    if fetch_list:
        return fetch_list

    raise ValueError("None of the requested accessions were eligible for inclusion")


def _check_isolate_size(
    plan: MonopartitePlan | MultipartitePlan,
    isolate_n: int,
) -> bool:
    """Return True if the size of the proposed isolate matches the isolate plan."""
    if type(plan) is MonopartitePlan:
        if isolate_n > 1:
            raise ValueError("Too many segments in monopartite isolate.")
        return True

    if isolate_n == len(plan.required_segments):
        return True

    raise ValueError(
        f"The plan requires {len(plan.required_segments)} segments: "
        + f"{[str(segment.name) for segment in plan.segments]}"
    )

from uuid import UUID

from structlog import get_logger

from ref_builder.ncbi.client import NCBIClient
from ref_builder.ncbi.models import NCBIGenbank
from ref_builder.otu.utils import (
    create_schema_from_records,
    group_genbank_records_by_isolate,
    parse_refseq_comment,
)
from ref_builder.repo import Repo
from ref_builder.resources import RepoIsolate, RepoOTU
from ref_builder.utils import IsolateName

logger = get_logger("otu.update")


def add_isolate(
    repo: Repo,
    otu: RepoOTU,
    accessions: list[str],
    ignore_cache: bool = False,
) -> RepoIsolate | None:
    """Take a list of accessions that make up a new isolate and a new isolate to the OTU.

    Download the GenBank records, categorize into an isolate bin and pass the isolate name
    and records to the add method.
    """
    ncbi = NCBIClient(ignore_cache)

    otu_logger = logger.bind(taxid=otu.taxid, otu_id=str(otu.id), name=otu.name)

    fetch_list = list(set(accessions).difference(otu.blocked_accessions))
    if not fetch_list:
        otu_logger.warning(
            "None of the requested accessions were eligible for inclusion.",
            requested_accessions=accessions,
            otu_accessions=otu.accessions,
            excluded_accessions=otu.excluded_accessions,
        )

        otu_logger.error(
            "Could not create a new isolate using the requested accessions.",
            requested_accessions=accessions,
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
        otu_logger.warning("More than one isolate name found in requested accession.")
        # Override later
        return None

    isolate_name, isolate_records = list(record_bins.items())[0]

    if len(isolate_records) != len(otu.schema.segments):
        otu_logger.error(
            f"The schema requires {len(otu.schema.segments)} segments: "
            + f"{[segment.name for segment in otu.schema.segments]}",
        )
        return None

    return create_isolate_from_records(
        repo,
        otu,
        isolate_name,
        list(isolate_records.values()),
    )


def create_isolate_from_records(
    repo: Repo,
    otu: RepoOTU,
    isolate_name: IsolateName,
    records: list[NCBIGenbank],
) -> RepoIsolate | None:
    """Take a list of GenBank records that make up a new isolate
    and add them to the OTU.
    """
    isolate_logger = get_logger("otu.isolate").bind(
        isolate_name=str(isolate_name),
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


def auto_update_otu(
    repo: Repo,
    otu: RepoOTU,
    ignore_cache: bool = False,
):
    """Fetch a full list of Nucleotide accessions associated with the OTU
    and pass the list to the add method.
    """
    ncbi = NCBIClient(ignore_cache)

    otu_logger = logger.bind(taxid=otu.taxid, otu_id=str(otu.id), name=otu.name)

    linked_accessions = ncbi.link_accessions_from_taxid(otu.taxid)

    fetch_list = list(set(linked_accessions).difference(otu.blocked_accessions))
    if not fetch_list:
        otu_logger.info("OTU is up to date.")
        return

    update_otu_with_accessions(repo, otu, fetch_list)


def update_otu_with_accessions(
    repo: Repo,
    otu: RepoOTU,
    accessions: list,
    ignore_cache: bool = False,
):
    """Take a list of accessions, filter for eligible accessions and
    add new sequences to the OTU.
    """
    otu_logger = logger.bind(taxid=otu.taxid, otu_id=str(otu.id), name=otu.name)

    ncbi = NCBIClient(ignore_cache)

    otu_logger.info(
        "Fetching accessions",
        count={len(accessions)},
        fetch_list=accessions,
    )

    records = ncbi.fetch_genbank_records(accessions)

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
    to the OTU."""
    otu_logger = logger.bind(taxid=otu.taxid, otu_id=str(otu.id), name=otu.name)

    record_bins = group_genbank_records_by_isolate(records)

    otu_logger.debug(
        f"Found {len(record_bins)} potential new isolates",
        potential_isolates=[str(isolate_name) for isolate_name in record_bins.keys()],
    )

    new_isolate_names = []

    for isolate_name in record_bins:
        isolate_records = record_bins[isolate_name]
        if len(isolate_records) != len(otu.schema.segments):
            otu_logger.debug(
                f"The schema requires {len(otu.schema.segments)} segments: "
                + f"{[segment.name for segment in otu.schema.segments]}",
                isolate_accessions=[accession.key for accession in isolate_records],
            )
            otu_logger.debug(f"Skipping {isolate_name}")
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
        otu_logger.info("New isolates added", new_isolates=new_isolate_names)

        return new_isolate_names

    otu_logger.info("No new isolates added.")
    return []


def exclude_accessions_from_otu(
    repo: Repo,
    otu: RepoOTU,
    accessions: list[str],
):
    """Take a list of accessions and add them to an OTU's excluded accessions list."""
    excluded_accessions = set()
    for accession in accessions:
        excluded_accessions = repo.exclude_accession(otu.id, accession)

    logger.debug(
        f"Accessions currently excluded from fetches: {excluded_accessions}",
        taxid=otu.taxid,
        otu_id=str(otu.id),
        name=otu.name,
    )


def add_schema_from_accessions(
    repo: Repo,
    taxid: int,
    accessions: list[str],
    ignore_cache: bool = False,
):
    """Take a list of accessions, create an OTU schema based on
    the corresponding Genbank data and add the new schema to the OTU.
    """
    if (otu := repo.get_otu_by_taxid(taxid)) is None:
        logger.fatal(f"OTU not found for {taxid}. Create first.")
        return

    otu_logger = logger.bind(otu_id=otu.id, taxid=taxid)

    if otu.schema is not None:
        logger.warning("OTU already has a schema attached.", schema=otu.schema)

    ncbi = NCBIClient(ignore_cache)

    records = ncbi.fetch_genbank_records(accessions)
    if not records:
        logger.fatal("Records could not be retrieved. Schema cannot be created.")
        return

    schema = create_schema_from_records(records)
    if schema is not None:
        otu_logger.info("Adding schema to OTU", schema=schema)
        repo.create_schema(
            otu_id=otu.id,
            molecule=schema.molecule,
            segments=schema.segments,
        )


def _bin_refseq_records(
    records: list[NCBIGenbank]) -> tuple[list[NCBIGenbank], list[NCBIGenbank]
]:
    """Returns a list of GenBank records as two lists, RefSeq and non-RefSeq."""
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

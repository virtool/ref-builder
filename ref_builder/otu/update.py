from structlog import get_logger

from ref_builder.ncbi.client import NCBIClient
from ref_builder.ncbi.models import NCBIGenbank
from ref_builder.otu.utils import (
    create_schema_from_records,
    group_genbank_records_by_isolate,
)
from ref_builder.repo import Repo
from ref_builder.resources import RepoOTU, RepoIsolate
from ref_builder.utils import IsolateName, Accession


logger = get_logger("otu.update")


def add_isolate(
    repo: Repo,
    otu: RepoOTU,
    accessions: list[str],
    ignore_cache: bool = False
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
            excluded_accessions=otu.excluded_accessions
        )

        otu_logger.error(
            "Could not create a new isolate using the requested accessions.",
            requested_accessions=accessions
        )

        return

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
        return

    isolate_name, isolate_records = list(record_bins.items())[0]

    if len(isolate_records) != len(otu.schema.segments):
        otu_logger.error(
            f"The schema requires {len(otu.schema.segments)} segments: "
            + f"{[segment.name for segment in otu.schema.segments]}"
        )
        return

    return create_isolate_from_records(
        repo, otu, isolate_name, list(isolate_records.values())
    )


def create_isolate_from_records(
    repo: Repo,
    otu: RepoOTU,
    isolate_name: IsolateName,
    records: list[NCBIGenbank],
):
    """Take a list of GenBank records that make up a new isolate
    and add them to the OTU."""
    isolate_logger = get_logger("otu.isolate").bind(
        isolate_name=str(isolate_name),
        otu_name=otu.name,
        taxid=otu.taxid,
    )

    if otu.get_isolate_id_by_name(isolate_name) is not None:
        isolate_logger.error(f"OTU already contains {isolate_name}.")
        return

    isolate = repo.create_isolate(
        otu.id, legacy_id=None, source_name=isolate_name.value, source_type=isolate_name.type
    )

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

    isolate_logger.info(
        f"{isolate.name} created",
        isolate_id=str(isolate.id),
        sequences=sorted([str(record.accession) for record in records]),
    )

    return isolate


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

    record_bins = group_genbank_records_by_isolate(records)

    otu_logger.debug(
        f"Found {len(record_bins)} potential new isolates",
        potential_isolates=[str(isolate_name) for isolate_name in record_bins.keys()],
    )

    new_isolates = []

    for isolate_name in record_bins:
        isolate_records = record_bins[isolate_name]
        if len(isolate_records) != len(otu.schema.segments):
            otu_logger.debug(
                f"The schema requires {len(otu.schema.segments)} segments: "
                + f"{[segment.name for segment in otu.schema.segments]}"
            )
            otu_logger.debug(f"Skipping {isolate_name}")
            continue

        isolate = create_isolate_from_records(
            repo, otu, isolate_name, list(isolate_records.values())
        )
        if isolate is not None:
            new_isolates.append(isolate_name)

    if not new_isolates:
        otu_logger.info("No new isolates added.")


def exclude_accessions_from_otu(
    repo: Repo,
    otu: RepoOTU,
    accessions: list[str],
):
    """Take a list of accessions and add them to an OTU's excluded accessions list."""
    for accession in accessions:
        repo.exclude_accession(otu.id, accession)


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


def _file_and_create_sequences(
    repo: Repo,
    otu: RepoOTU,
    records: list[NCBIGenbank],
) -> list[Accession]:
    otu_logger = logger.bind(taxid=otu.taxid, otu_id=str(otu.id), name=otu.name)

    record_bins = group_genbank_records_by_isolate(records)

    new_accessions = []

    for isolate_key in record_bins:
        record_bin = record_bins[isolate_key]

        isolate_id = otu.get_isolate_id_by_name(isolate_key)
        if isolate_id is None:
            otu_logger.debug(
                f"Creating isolate for {isolate_key.type}, {isolate_key.value}",
            )
            isolate = repo.create_isolate(
                otu_id=otu.id,
                legacy_id=None,
                source_name=isolate_key.value,
                source_type=isolate_key.type,
            )
            isolate_id = isolate.id

        for accession in record_bin:
            record = record_bin[accession]
            if accession.key in otu.accessions:
                extant_sequence = otu.get_sequence_by_accession(accession.key)
                if extant_sequence.accession.version == accession.version:
                    otu_logger.warning(f"{accession} already exists in OTU")
                    continue

                otu_logger.warning(
                    f"New version of {accession.key} found. Replacing {extant_sequence.accession} with {accession}...")
                otu_logger.warning("SEQUENCE REPLACEMENT NOT YET IMPLEMENTED")

            else:
                sequence = repo.create_sequence(
                    otu_id=otu.id,
                    isolate_id=isolate_id,
                    accession=record.accession_version,
                    definition=record.definition,
                    legacy_id=None,
                    segment=record.source.segment,
                    sequence=record.sequence,
                )

                new_accessions.append(sequence.accession)

    return sorted(new_accessions)

import sys
from collections import defaultdict

import structlog

from ref_builder.models import Molecule
from ref_builder.ncbi.client import NCBIClient
from ref_builder.ncbi.models import NCBIGenbank
from ref_builder.repo import Repo
from ref_builder.resources import RepoOTU
from ref_builder.schema import OTUSchema, Segment
from ref_builder.utils import Accession, IsolateName, IsolateNameType

logger = structlog.get_logger("otu")


def create_otu(
    repo: Repo,
    taxid: int,
    accessions: list[str],
    acronym: str,
    ignore_cache: bool = False,
) -> RepoOTU | None:
    """Create a new OTU by taxonomy ID and autogenerate a schema.

    Uses the provided accessions to generate a schema and add a first isolate.

    :param repo: the repository to add the OTU to.
    :param taxid: the taxonomy ID to use.
    :param accessions: a list of accessions to use for the schema.
    :param acronym: an alternative name to use during searches.
    :param ignore_cache: whether to ignore the cache.

    """
    otu_logger = logger.bind(taxid=taxid)

    client = NCBIClient(ignore_cache)

    if repo.get_otu_id_by_taxid(taxid):
        raise ValueError(
            f"Taxonomy ID {taxid} has already been added to this reference.",
        )

    taxonomy = client.fetch_taxonomy_record(taxid)
    if taxonomy is None:
        otu_logger.fatal(f"Could not retrieve {taxid} from NCBI Taxonomy")
        return None

    if not acronym:
        if taxonomy.other_names.acronym:
            acronym = taxonomy.other_names.acronym[0]

    records = client.fetch_genbank_records(accessions)
    if len(records) != len(accessions):
        otu_logger.fatal("Could not retrieve all requested accessions.")
        return None

    binned_records = group_genbank_records_by_isolate(records)

    if len(binned_records) > 1:
        otu_logger.fatal(
            "More than one isolate found. Cannot create schema automatically.",
        )
        return None

    schema = create_schema_from_records(records)

    if schema is None:
        otu_logger.fatal("Could not create schema.")
        return None

    try:
        otu = repo.create_otu(
            acronym=acronym,
            legacy_id=None,
            name=taxonomy.name,
            schema=schema,
            taxid=taxid,
        )

    except ValueError as e:
        otu_logger.fatal(e)
        sys.exit(1)

    isolate_name = list(binned_records.keys())[0]

    isolate = repo.create_isolate(
        otu_id=otu.id,
        legacy_id=None,
        source_name=isolate_name.value,
        source_type=isolate_name.type,
    )

    otu.add_isolate(isolate)

    otu.repr_isolate = repo.set_repr_isolate(otu_id=otu.id, isolate_id=isolate.id)

    for record in records:
        sequence = repo.create_sequence(
            otu_id=otu.id,
            isolate_id=isolate.id,
            accession=record.accession_version,
            definition=record.definition,
            legacy_id=None,
            segment=record.source.segment,
            sequence=record.sequence,
        )
        otu.add_sequence(sequence, isolate_id=isolate.id)

    return otu


def update_otu(
    repo: Repo,
    otu: RepoOTU,
    ignore_cache: bool = False,
):
    """Fetch a full list of Nucleotide accessions associated with the OTU
    and pass the list to the add method.
    """
    ncbi = NCBIClient(ignore_cache)

    linked_accessions = ncbi.link_accessions_from_taxid(otu.taxid)

    add_sequences(repo, otu, linked_accessions)


def add_sequences(
    repo: Repo,
    otu: RepoOTU,
    accessions: list,
    ignore_cache: bool = False,
):
    """Take a list of accessions, filter for eligible accessions and
    add new sequences to the OTU.
    """
    ncbi = NCBIClient(ignore_cache)

    otu_logger = logger.bind(taxid=otu.taxid, otu_id=str(otu.id), name=otu.name)
    fetch_list = list(set(accessions).difference(otu.blocked_accessions))
    if not fetch_list:
        otu_logger.info("OTU is up to date.")
        return

    otu_logger.info(
        "Fetching accessions",
        count={len(fetch_list)},
        fetch_list=fetch_list,
    )

    records = ncbi.fetch_genbank_records(fetch_list)

    new_accessions = _file_and_create_sequences(repo, otu, records)

    if new_accessions:
        otu_logger.info(
            "Added  sequences to OTU",
            count=len(new_accessions),
            new_accessions=new_accessions,
            taxid=otu.taxid,
        )
        return

    otu_logger.info("No new sequences added to OTU")


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


def create_schema_from_records(
    records: list[NCBIGenbank],
    segments: list[Segment] | None = None,
) -> OTUSchema | None:
    molecule = _get_molecule_from_records(records)

    binned_records = group_genbank_records_by_isolate(records)
    if len(binned_records) > 1:
        logger.fatal(
            "More than one isolate found. Cannot create schema automatically.",
            bins=binned_records,
        )
        return None

    if segments is None:
        segments = _get_segments_from_records(records)
        if segments:
            return OTUSchema(molecule=molecule, segments=segments)

    return None


def _get_molecule_from_records(records: list[NCBIGenbank]) -> Molecule:
    """Return relevant molecule metadata from one or more records"""
    if not records:
        raise IndexError("No records given")

    rep_record = None
    for record in records:
        if record.refseq:
            rep_record = record
            break
    if rep_record is None:
        rep_record = records[0]

    return Molecule.model_validate(
        {
            "strandedness": rep_record.strandedness.value,
            "type": rep_record.moltype.value,
            "topology": rep_record.topology.value,
        },
    )


def group_genbank_records_by_isolate(
    records: list[NCBIGenbank],
) -> dict[IsolateName, dict[Accession, NCBIGenbank]]:
    """Indexes Genbank records by isolate name"""
    isolates = defaultdict(dict)

    for record in records:
        if (isolate_name := _get_isolate_name(record)) is not None:
            versioned_accession = Accession.create_from_string(record.accession_version)
            isolates[isolate_name][versioned_accession] = record

    return isolates


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


def _get_isolate_name(record: NCBIGenbank) -> IsolateName | None:
    """Get the isolate name from a Genbank record"""
    record_logger = logger.bind(
        accession=record.accession,
        definition=record.definition,
        source_data=record.source,
    )

    if record.source.model_fields_set.intersection(
        {IsolateNameType.ISOLATE, IsolateNameType.STRAIN, IsolateNameType.CLONE},
    ):
        for source_type in IsolateNameType:
            if source_type in record.source.model_fields_set:
                return IsolateName(
                    type=IsolateNameType(source_type),
                    value=record.source.model_dump()[source_type],
                )

    elif record.refseq:
        record_logger.debug(
            "RefSeq record does not contain sufficient source data. "
            + "Edit before inclusion.",
        )

        return IsolateName(
            type=IsolateNameType(IsolateNameType.REFSEQ),
            value=record.accession,
        )

    record_logger.debug("Record does not contain sufficient source data for inclusion.")

    return None


def _get_segments_from_records(records: list[NCBIGenbank]) -> list[Segment]:
    if len(records) == 1:
        record = records[0]

        if record.source.segment != "":
            segment_name = record.source.segment
        else:
            segment_name = record.source.organism

        return [Segment(name=segment_name, required=True, length=len(record.sequence))]

    segments = []
    for record in sorted(records, key=lambda record: record.accession):
        if record.source.segment:
            segments.append(
                Segment(
                    name=record.source.segment,
                    required=True,
                    length=len(record.sequence),
                ),
            )
        else:
            raise ValueError("No segment name found for multipartite OTU segment.")

    return segments

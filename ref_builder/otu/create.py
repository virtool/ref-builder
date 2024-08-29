import sys

import structlog

from ref_builder.ncbi.client import NCBIClient
from ref_builder.otu.utils import (
    create_schema_from_records,
    group_genbank_records_by_isolate,
    parse_refseq_comment,
)
from ref_builder.repo import Repo
from ref_builder.resources import RepoOTU

logger = structlog.get_logger("otu.create")


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
        name=isolate_name,
    )

    otu.add_isolate(isolate)

    otu.repr_isolate = repo.set_repr_isolate(otu_id=otu.id, isolate_id=isolate.id)

    for record in records:
        repo.create_sequence(
            otu_id=otu.id,
            isolate_id=isolate.id,
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

    return repo.get_otu(otu.id)

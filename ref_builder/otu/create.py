import sys

import structlog

from ref_builder.ncbi.client import NCBIClient
from ref_builder.ncbi.models import NCBIGenbank, NCBITaxonomy
from ref_builder.otu.isolate import create_sequence_from_record
from ref_builder.otu.utils import (
    assign_records_to_segments,
    create_plan_from_records,
    get_molecule_from_records,
    group_genbank_records_by_isolate,
    parse_refseq_comment,
)
from ref_builder.repo import Repo
from ref_builder.resources import RepoOTU
from ref_builder.utils import IsolateName

logger = structlog.get_logger("otu.create")


def create_otu_with_taxid(
    repo: Repo,
    taxid: int,
    accessions: list[str],
    acronym: str,
    ignore_cache: bool = False,
) -> RepoOTU | None:
    """Create a new OTU by taxonomy ID.

    Uses the provided accessions to generate a plan and add a first isolate.

    :param repo: the repository to add the OTU to.
    :param taxid: the taxonomy ID to use.
    :param accessions: accessions to build the new otu from
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

    if not acronym and taxonomy.other_names.acronym:
        acronym = taxonomy.other_names.acronym[0]

    records = client.fetch_genbank_records(accessions)

    if len(records) != len(accessions):
        otu_logger.fatal("Could not retrieve all requested accessions.")
        return None

    binned_records = group_genbank_records_by_isolate(records)

    if len(binned_records) > 1:
        otu_logger.fatal(
            "More than one isolate found. Cannot create plan.",
        )
        return None

    with repo.use_transaction():
        return write_otu(
            repo,
            taxonomy,
            records,
            acronym=acronym,
            isolate_name=next(iter(binned_records.keys())) if binned_records else None,
        )


def create_otu_without_taxid(
    repo: Repo,
    accessions: list[str],
        acronym: str,
        ignore_cache: bool = False,
) -> RepoOTU | None:
    """Create a new OTU from a list of accessions.

    Uses the provided accessions to generate a plan and add a first isolate.

    :param repo: the repository to add the OTU to.
    :param accessions: accessions to build the new otu from
    :param acronym: an alternative name to use during searches.
    :param ignore_cache: whether to ignore the cache.

    """
    ncbi = NCBIClient(ignore_cache)

    records = ncbi.fetch_genbank_records(accessions)

    if len(records) != len(accessions):
        logger.fatal("Could not retrieve all requested accessions.")
        return None

    if len(set(record.source.taxid for record in records)) > 1:
        logger.fatal("Not all records are from the same organism.")

        return None

    taxid = records[0].source.taxid

    binned_records = group_genbank_records_by_isolate(records)

    if len(binned_records) > 1:
        logger.fatal(
            "More than one isolate found. Cannot create plan.",
        )
        return None

    taxonomy = ncbi.fetch_taxonomy_record(taxid)

    if taxonomy is None:
        logger.fatal(f"Could not retrieve {taxid} from NCBI Taxonomy")
        return None

    with repo.use_transaction():
        return write_otu(
            repo,
            taxonomy,
            records,
            acronym,
            isolate_name=next(iter(binned_records.keys())) if binned_records else None,
        )


def write_otu(
    repo: Repo,
    taxonomy: NCBITaxonomy,
    records: list[NCBIGenbank],
    acronym: str,
    isolate_name: IsolateName | None,
) -> RepoOTU | None:
    """Create a new OTU from an NCBI Taxonomy record and a list of Nucleotide records."""
    otu_logger = logger.bind(taxid=taxonomy.id)

    plan = create_plan_from_records(
        records,
        length_tolerance=repo.settings.default_segment_length_tolerance,
    )

    if plan is None:
        logger.fatal("Could not create plan from records.")

    molecule = get_molecule_from_records(records)

    try:
        otu = repo.create_otu(
            acronym=acronym,
            legacy_id=None,
            molecule=molecule,
            name=taxonomy.name,
            plan=plan,
            taxid=taxonomy.id,
        )
    except ValueError as e:
        otu_logger.fatal(e)
        sys.exit(1)

    isolate = repo.create_isolate(
        otu_id=otu.id,
        legacy_id=None,
        name=isolate_name,
    )

    otu.add_isolate(isolate)
    otu.representative_isolate = repo.set_representative_isolate(
        otu_id=otu.id, isolate_id=isolate.id
    )

    if otu.plan.monopartite:
        record = records[0]

        sequence = create_sequence_from_record(repo, otu, record, plan.segments[0].id)

        repo.link_sequence(otu.id, isolate.id, sequence.id)

        if record.refseq:
            _, old_accession = parse_refseq_comment(record.comment)

            repo.exclude_accession(
                otu.id,
                old_accession,
            )

    else:
        for segment_id, record in assign_records_to_segments(records, plan).items():
            sequence = create_sequence_from_record(repo, otu, record, segment_id)

            repo.link_sequence(otu.id, isolate.id, sequence.id)

            if record.refseq:
                _, old_accession = parse_refseq_comment(record.comment)
                repo.exclude_accession(
                    otu.id,
                    old_accession,
                )

    return repo.get_otu(otu.id)

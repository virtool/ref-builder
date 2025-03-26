from uuid import UUID

from structlog import get_logger

from ref_builder.ncbi.client import NCBIClient
from ref_builder.ncbi.models import NCBIGenbank
from ref_builder.repo import Repo
from ref_builder.resources import RepoOTU
from ref_builder.otu.utils import (
    DeleteRationale,
    assign_segment_id_to_record,
    parse_refseq_comment, get_segments_min_length, get_segments_max_length,
)


logger = get_logger("otu.promote")


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

        if promoted_accessions := promote_otu_accessions_from_records(repo, otu, records):
            log.info("Sequences promoted.", new_accessions=sorted(promoted_accessions))

    log.info("Records are already up to date.")


def promote_otu_accessions_from_records(
    repo: Repo, otu: RepoOTU, records: list[NCBIGenbank]
) -> set[str]:
    """Take a list of records and check them against the contents of an OTU
    for promotable RefSeq sequences. Return a list of promoted accessions.
    """
    otu_logger = logger.bind(otu_id=str(otu.id), taxid=otu.taxid)

    refseq_records = [record for record in records if record.refseq]

    records_by_promotable_sequence_id = {}

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

            predecessor_sequence = otu.get_sequence_by_accession(predecessor_accession)

            records_by_promotable_sequence_id[predecessor_sequence.id] = record

    promoted_sequence_ids = replace_accessions_from_records(
        repo, otu, records_by_promotable_sequence_id
    )

    otu = repo.get_otu(otu.id)

    return set(
        otu.get_sequence_by_id(sequence_id).accession.key
        for sequence_id in promoted_sequence_ids
    )


def replace_accessions_from_records(
    repo: Repo,
    otu: RepoOTU,
    record_by_replaceable_sequence_id: dict[UUID, NCBIGenbank],
):
    initial_exceptions = otu.excluded_accessions.copy()

    replacement_sequence_ids = set()

    for sequence_id in record_by_replaceable_sequence_id:
        predecessor_sequence = otu.get_sequence_by_id(sequence_id)
        if predecessor_sequence is None:
            logger.warning("Predecessor sequences not found")

        containing_isolate_ids = otu.get_isolate_ids_containing_sequence_id(
            predecessor_sequence.id
        )
        if not containing_isolate_ids:
            logger.info("Sequence id not found in any isolates.")
            continue

        logger.debug(
            "Isolates containing sequence found",
            replaceable_sequence=str(predecessor_sequence.id),
            isolate_ids=[str(isolate_id) for isolate_id in containing_isolate_ids],
        )

        replacement_record = record_by_replaceable_sequence_id[sequence_id]

        segment_id = assign_segment_id_to_record(replacement_record, otu.plan)
        if segment_id is None:
            logger.error("This segment does not match the plan.")
            continue

        with repo.use_transaction():
            if otu.get_sequence_by_accession(replacement_record.accession) is None:
                replacement_sequence = repo.create_sequence(
                    otu.id,
                    accession=replacement_record.accession_version,
                    definition=replacement_record.definition,
                    legacy_id=None,
                    segment=segment_id,
                    sequence=replacement_record.sequence,
                )

                if replacement_sequence is None:
                    logger.error("Isolate update failed when creating new sequence.")
                    return None

            else:
                replacement_sequence = otu.get_sequence_by_accession(
                    replacement_record.accession
                )

            for isolate_id in containing_isolate_ids:
                repo.replace_sequence(
                    otu.id,
                    isolate_id,
                    replacement_sequence.id,
                    replaced_sequence_id=predecessor_sequence.id,
                    rationale=DeleteRationale.REFSEQ,
                )

            repo.exclude_accession(otu.id, predecessor_sequence.accession.key)

            replacement_sequence_ids.add(replacement_sequence.id)

        otu = repo.get_otu(otu.id)

    if replacement_sequence_ids:
        otu = repo.get_otu(otu.id)

        replaced_sequence_index = (
            {
                str(otu.get_sequence_by_id(sequence_id).accession): str(sequence_id)
                for sequence_id in replacement_sequence_ids
            },
        )

        logger.info(
            "Replaced sequences",
            count=len(replacement_sequence_ids),
            replaced_sequences=replaced_sequence_index,
            new_excluded_accessions=sorted(
                otu.excluded_accessions - initial_exceptions
            ),
        )

    return replacement_sequence_ids

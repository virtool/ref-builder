from uuid import UUID

from structlog import get_logger

from ref_builder.ncbi.models import NCBIGenbank
from ref_builder.repo import Repo
from ref_builder.resources import RepoOTU
from ref_builder.otu.utils import DeleteRationale, assign_segment_id_to_record


logger = get_logger("otu.promote")


def replace_accessions_from_records(
    repo: Repo,
    otu: RepoOTU,
    record_by_replaceable_sequence_id: dict[UUID, NCBIGenbank],
):
    initial_exceptions = otu.excluded_accessions.copy()

    replacement_sequence_ids = set()

    for sequence_id in record_by_replaceable_sequence_id:
        predecessor_sequence = otu.get_sequence_by_id(sequence_id)

        containing_isolate_ids = otu.get_isolate_ids_containing_sequence_id(sequence_id)

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
                replacement_sequence = otu.get_sequence_by_accession(replacement_record.accession)

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

    if replacement_sequence_ids:
        otu = repo.get_otu(otu.id)

        logger.info(
            "Replaced sequences",
            count=len(replacement_sequence_ids),
            replaced_sequence_ids=[str(sequence_id) for sequence_id in replacement_sequence_ids],
            new_excluded_accessions=sorted(otu.excluded_accessions - initial_exceptions),
        )

    return replacement_sequence_ids






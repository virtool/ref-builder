from uuid import UUID

from pydantic import ValidationError
from structlog import get_logger

from ref_builder.ncbi.client import NCBIClient
from ref_builder.otu.utils import DeleteRationale, create_segments_from_records
from ref_builder.plan import (
    Plan,
    Segment,
    SegmentName,
    SegmentRule,
)
from ref_builder.repo import Repo
from ref_builder.resources import RepoOTU, RepoSequence

logger = get_logger("otu.modify")


def exclude_accessions_from_otu(
    repo: Repo,
    otu: RepoOTU,
    accessions: list[str],
) -> None:
    """Take a list of accessions and add them to an OTU's excluded accessions list."""
    accession_set = set(accessions)

    for accession in accession_set:
        repo.exclude_accession(otu.id, accession)


def delete_isolate_from_otu(repo: Repo, otu: RepoOTU, isolate_id: UUID) -> None:
    """Remove an isolate from a specified OTU."""
    otu_logger = logger.bind(otu_id=str(otu.id), taxid=otu.taxid)

    if isolate := otu.get_isolate(isolate_id):
        repo.delete_isolate(otu.id, isolate.id, rationale=DeleteRationale.USER)

        otu_logger.info(
            f"{isolate.name} removed.",
            removed_isolate_id=isolate_id,
            current_isolate_ids=list[otu.isolate_ids],
        )

    logger.error(
        "This isolate does not exist in this OTU.",
        isolate_id=isolate_id,
        current_isolate_ids=list[otu.isolate_ids],
    )


def set_plan(
    repo: Repo,
    otu: RepoOTU,
    plan: Plan,
) -> Plan | None:
    """Sets an OTU's isolate plan to a new plan."""
    otu_logger = logger.bind(name=otu.name, taxid=otu.taxid, plan=plan.model_dump())

    try:
        repo.set_plan(otu.id, plan)
    except ValueError as e:
        otu_logger.error(e)
        return None

    return repo.get_otu(otu.id).plan


def set_plan_length_tolerances(
    repo: Repo,
    otu: RepoOTU,
    tolerance: float,
) -> Plan | None:
    """Sets a plan's length tolerances to a new float value."""
    if otu.plan.monopartite:
        try:
            repo.set_plan(
                otu.id,
                plan=Plan.new(
                    segments=[
                        Segment.new(
                            length=otu.plan.segments[0].length,
                            length_tolerance=tolerance,
                            name=None,
                            required=SegmentRule.REQUIRED,
                        )
                    ]
                ),
            )
        except (ValueError, ValidationError) as e:
            logger.error(
                e,
                name=otu.name,
                taxid=otu.taxid,
                tolerance=otu.plan.segments[0].length_tolerance,
                new_tolerance=tolerance,
            )
            return None

    return repo.get_otu(otu.id).plan


def add_segments_to_plan(
    repo: Repo,
    otu: RepoOTU,
    rule: SegmentRule,
    accessions: list[str],
    ignore_cache: bool = False,
) -> Plan | None:
    """Add new segments to a multipartite plan."""
    expand_logger = logger.bind(
        name=otu.name, taxid=otu.taxid, accessions=accessions, rule=rule
    )

    sequence_length_tolerance = repo.settings.default_segment_length_tolerance

    if otu.plan.monopartite and otu.plan.segments[0].name is None:
        raise ValueError("Current monopartite segment is unnamed.")

    client = NCBIClient(ignore_cache)

    expand_logger.info(
        f"Adding {len(accessions)} sequences to plan as {rule} segments",
        current_plan=otu.plan.model_dump(),
    )

    records = client.fetch_genbank_records(accessions)
    if not records:
        expand_logger.warning("Could not fetch records.")
        return None

    new_segments = create_segments_from_records(
        records,
        rule,
        length_tolerance=sequence_length_tolerance,
    )
    if not new_segments:
        expand_logger.warning("No segments can be added.")
        return None

    new_plan = Plan.new(
        segments=otu.plan.segments + new_segments,
    )

    return set_plan(repo, otu, new_plan)


def rename_plan_segment(
    repo: Repo,
    otu: RepoOTU,
    segment_id: UUID,
    segment_name: SegmentName,
) -> Plan | None:
    """Set a new name for the segment matching the given ID if possible, else return None."""
    rename_logger = logger.bind(
        name=otu.name,
        taxid=otu.taxid,
        segment_id=str(segment_id),
        segment_name=str(segment_name),
    )

    modified_plan = otu.plan.model_copy()

    for segment in modified_plan.segments:
        if segment.id == segment_id:
            segment.name = segment_name

            return set_plan(repo, otu, modified_plan)

    rename_logger.error("Segment with requested ID not found.")

    return None


def resize_monopartite_plan(
    repo: Repo,
    otu: RepoOTU,
    name: SegmentName,
    rule: SegmentRule,
    accessions: list[str],
    ignore_cache: bool = False,
) -> Plan | None:
    """Convert a monopartite plan to a multipartite plan and add segments."""
    expand_logger = logger.bind(
        name=otu.name, taxid=otu.taxid, accessions=accessions, rule=rule
    )

    if not otu.plan.monopartite:
        expand_logger.warning("OTU plan is already a multipartite plan.")
        return None

    sequence_length_tolerance = repo.settings.default_segment_length_tolerance

    client = NCBIClient(ignore_cache)

    records = client.fetch_genbank_records(accessions)
    if not records:
        expand_logger.warning("Could not fetch records.")
        return None

    new_segments = create_segments_from_records(
        records,
        rule,
        length_tolerance=sequence_length_tolerance,
    )
    if not new_segments:
        expand_logger.warning("No segments can be added.")
        return None

    new_plan = Plan.new(
        segments=[
            Segment.new(
                length=otu.plan.segments[0].length,
                length_tolerance=sequence_length_tolerance,
                name=name,
                required=SegmentRule.REQUIRED,
            )
        ]
        + new_segments,
    )

    return set_plan(repo, otu, new_plan)


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

    new_sequence = repo.create_sequence(
        otu.id,
        accession=record.accession_version,
        definition=record.definition,
        legacy_id=None,
        segment=record.source.segment,
        sequence=record.sequence,
    )

    repo.replace_sequence(
        otu.id,
        isolate_id=isolate.id,
        sequence_id=new_sequence.id,
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

    if otu.representative_isolate is not None:
        if otu.representative_isolate == new_representative_isolate.id:
            otu_logger.warning("This isolate is already the representative isolate.")
            return otu.representative_isolate

        otu_logger.warning(
            f"Replacing representative isolate {otu.representative_isolate}",
            representative_isolate_id=str(otu.representative_isolate),
        )

    repo.set_representative_isolate(otu.id, new_representative_isolate.id)

    otu_logger.info(
        f"Representative isolate set to {new_representative_isolate.name}.",
        representative_isolate_id=str(new_representative_isolate.id),
    )

    return new_representative_isolate.id

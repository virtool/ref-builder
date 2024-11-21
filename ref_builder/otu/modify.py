from uuid import UUID

from pydantic import ValidationError
from structlog import get_logger

from ref_builder.ncbi.client import NCBIClient
from ref_builder.otu.utils import DeleteRationale, create_segments_from_records
from ref_builder.plan import (
    MonopartitePlan,
    MultipartitePlan,
    SegmentRule,
    SegmentName,
    Segment,
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


def set_isolate_plan(
    repo: Repo,
    otu: RepoOTU,
    plan: MonopartitePlan | MultipartitePlan,
) -> MonopartitePlan | MultipartitePlan | None:
    """Sets an OTU's isolate plan to a new plan."""
    otu_logger = logger.bind(name=otu.name, taxid=otu.taxid, plan=plan.model_dump())

    try:
        repo.set_isolate_plan(otu.id, plan)
    except ValueError as e:
        otu_logger.error(e)
        return None

    return repo.get_otu(otu.id).plan


def set_plan_length_tolerances(
    repo: Repo,
    otu: RepoOTU,
    tolerance: float,
) -> MonopartitePlan | MultipartitePlan | None:
    """Sets a plan's length tolerances to a new float value."""
    if otu.plan.plan_type == "monopartite":
        try:
            repo.set_isolate_plan(
                otu.id,
                MonopartitePlan.new(
                    length=otu.plan.length,
                    length_tolerance=tolerance,
                ),
            )
        except (ValueError, ValidationError) as e:
            logger.error(
                e,
                name=otu.name,
                taxid=otu.taxid,
                tolerance=otu.plan.length_tolerance,
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
) -> MultipartitePlan | None:
    """Add new segments to a multipartite plan."""
    expand_logger = logger.bind(
        name=otu.name, taxid=otu.taxid, accessions=accessions, rule=rule
    )

    sequence_length_tolerance = repo.settings.default_segment_length_tolerance

    if type(otu.plan) is MonopartitePlan:
        raise ValueError("Cannot add segments to a monopartite plan.")

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

    new_plan = MultipartitePlan.new(
        segments=otu.plan.segments + new_segments,
    )

    return set_isolate_plan(repo, otu, new_plan)


def resize_monopartite_plan(
    repo: Repo,
    otu: RepoOTU,
    name: SegmentName,
    rule: SegmentRule,
    accessions: list[str],
    # sequence_length_tolerance: float | None = None,
    ignore_cache: bool = False,
) -> MultipartitePlan | None:
    """Convert a monopartite plan to a multipartite plan and add segments."""
    expand_logger = logger.bind(
        name=otu.name, taxid=otu.taxid, accessions=accessions, rule=rule
    )

    if type(otu.plan) is MultipartitePlan:
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

    new_plan = MultipartitePlan.new(
        segments=[
            Segment.new(
                length=otu.plan.length,
                length_tolerance=sequence_length_tolerance,
                name=name,
                required=SegmentRule.REQUIRED,
            )
        ]
        + new_segments,
    )

    return set_isolate_plan(repo, otu, new_plan)


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

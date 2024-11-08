from uuid import uuid4

from ref_builder.ncbi.client import NCBIClient
from ref_builder.otu.update import logger
from ref_builder.otu.utils import create_segments_from_records
from ref_builder.plan import (
    MonopartitePlan,
    MultipartitePlan,
    SegmentRule,
    SegmentName,
    Segment,
)
from ref_builder.repo import Repo
from ref_builder.resources import RepoOTU


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
    )
    if not new_segments:
        expand_logger.warning("No segments can be added.")
        return None

    new_plan = MultipartitePlan.new(
        segments=otu.plan.segments + new_segments,
    )

    try:
        repo.set_isolate_plan(otu.id, new_plan)
    except ValueError as e:
        expand_logger.error(e)
        return None

    expand_logger.info(
        f"Segments {[str(segment.name) for segment in new_segments]} added to plan",
        expanded_plan=new_plan.model_dump(),
    )

    return repo.get_otu(otu.id).plan


def resize_monopartite_plan(
    repo: Repo,
    otu: RepoOTU,
    name: SegmentName,
    rule: SegmentRule,
    accessions: list[str],
    sequence_length_tolerance: float | None = None,
    ignore_cache: bool = False,
) -> MultipartitePlan | None:
    """Convert a monopartite plan to a multipartite plan and add segments."""
    expand_logger = logger.bind(
        name=otu.name, taxid=otu.taxid, accessions=accessions, rule=rule
    )

    if type(otu.plan) is MultipartitePlan:
        expand_logger.warning("OTU plan is already a multipartite plan.")
        return None

    if sequence_length_tolerance is None:
        sequence_length_tolerance = repo.settings.default_segment_length_tolerance

    client = NCBIClient(ignore_cache)

    records = client.fetch_genbank_records(accessions)
    if not records:
        expand_logger.warning("Could not fetch records.")
        return None

    new_segments = create_segments_from_records(
        records,
        rule,
    )
    if not new_segments:
        expand_logger.warning("No segments can be added.")
        return None

    new_plan = MultipartitePlan.new(
        segments=[
            Segment(
                id=uuid4(),
                length=otu.plan.length,
                name=name,
                required=SegmentRule.REQUIRED,
            )
        ]
        + new_segments,
    )

    try:
        repo.set_isolate_plan(otu.id, new_plan)
    except ValueError as e:
        expand_logger.error(e)
        return None

    expand_logger.info(
        f"Segments {[str(segment.name) for segment in new_segments]} added to plan",
        expanded_plan=new_plan.model_dump(),
    )

    return repo.get_otu(otu.id).plan

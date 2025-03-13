from collections.abc import Collection
from uuid import UUID

from pydantic import ValidationError
from structlog import get_logger

from ref_builder.ncbi.client import NCBIClient
from ref_builder.otu.utils import (
    DeleteRationale,
    assign_records_to_segments,
    create_segments_from_records,
)
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
    accessions: Collection[str],
) -> None:
    """Exclude accessions from future addition to an OTU."""
    original_excluded_accessions = otu.excluded_accessions.copy()

    with repo.use_transaction():
        excluded_accessions = repo.exclude_accessions(
            otu_id=otu.id, accessions=accessions
        )

    if excluded_accessions != original_excluded_accessions:
        logger.info(
            "Updated excluded accession list.",
            otu_id=str(otu.id),
            excluded_accessions=sorted(excluded_accessions),
        )

    if excluded_accessions == original_excluded_accessions:
        logger.info(
            "Excluded accession list already up to date.",
            excluded_accessions=excluded_accessions,
        )


def allow_accessions_into_otu(
    repo: Repo,
    otu: RepoOTU,
    accessions: Collection[str],
) -> None:
    """Allow accessions for future addition to an OTU.

    This reverses the effect of exclude_accessions_from_otu.
    """
    original_excluded_accessions = otu.excluded_accessions.copy()

    with repo.use_transaction():
        excluded_accessions = repo.allow_accessions(
            otu_id=otu.id, accessions=accessions
        )

    if excluded_accessions != original_excluded_accessions:
        logger.info(
            "Updated excluded accession list.",
            otu_id=str(otu.id),
            excluded_accessions=sorted(excluded_accessions),
        )

    if excluded_accessions == original_excluded_accessions:
        logger.info(
            "Excluded accession list already up to date.",
            excluded_accessions=excluded_accessions,
        )


def delete_isolate_from_otu(repo: Repo, otu: RepoOTU, isolate_id: UUID) -> bool:
    """Remove an isolate from a specified OTU."""
    otu_logger = logger.bind(otu_id=str(otu.id), taxid=otu.taxid)

    if isolate_id == otu.representative_isolate:
        otu_logger.error(
            "The representative isolate cannot be deleted from the OTU.",
            isolate_id=str(isolate_id),
        )
        return False

    isolate = otu.get_isolate(isolate_id)

    if not isolate:
        otu_logger.error("Isolate not found.", isolate_id=str(isolate_id))
        return False

    with repo.use_transaction():
        repo.delete_isolate(otu.id, isolate.id, rationale=DeleteRationale.USER)

    otu_logger.info(
        "Isolate removed.",
        name=isolate.name,
        removed_isolate_id=str(isolate_id),
        current_isolate_ids=[str(isolate_id_) for isolate_id_ in otu.isolate_ids],
    )

    repo.get_otu(otu.id)

    return True


def set_plan(
    repo: Repo,
    otu: RepoOTU,
    plan: Plan,
) -> Plan | None:
    """Set an OTU's plan to the passed ``plan``."""
    log = logger.bind(name=otu.name, taxid=otu.taxid, plan=plan.model_dump())

    try:
        with repo.use_transaction():
            repo.set_plan(otu.id, plan)
    except ValueError:
        log.exception()
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
            with repo.use_transaction():
                repo.set_plan(
                    otu.id,
                    plan=Plan.new(
                        segments=[
                            Segment.new(
                                length=otu.plan.segments[0].length,
                                length_tolerance=tolerance,
                                name=None,
                                rule=SegmentRule.REQUIRED,
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
    log = logger.bind(name=otu.name, taxid=otu.taxid, accessions=accessions, rule=rule)

    sequence_length_tolerance = repo.settings.default_segment_length_tolerance

    if otu.plan.monopartite and otu.plan.segments[0].name is None:
        raise ValueError("Current monopartite segment is unnamed.")

    client = NCBIClient(ignore_cache)

    log.info("Adding sequences to plan as segments", count=len(accessions), rule=rule)

    records = client.fetch_genbank_records(accessions)
    if not records:
        log.warning("Could not fetch records.")
        return None

    new_segments = create_segments_from_records(
        records,
        rule,
        length_tolerance=sequence_length_tolerance,
    )
    if not new_segments:
        log.warning("No segments can be added.")
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
    log = logger.bind(
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

    log.error("Segment with requested ID not found.")

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
    log = logger.bind(name=otu.name, taxid=otu.taxid, accessions=accessions, rule=rule)

    if not otu.plan.monopartite:
        log.warning("OTU plan is already a multipartite plan.")
        return None

    sequence_length_tolerance = repo.settings.default_segment_length_tolerance

    client = NCBIClient(ignore_cache)

    records = client.fetch_genbank_records(accessions)
    if not records:
        log.warning("Could not fetch records.")
        return None

    new_segments = create_segments_from_records(
        records,
        rule,
        length_tolerance=sequence_length_tolerance,
    )
    if not new_segments:
        log.warning("No segments can be added.")
        return None

    new_plan = Plan.new(
        segments=[
            Segment.new(
                length=otu.plan.segments[0].length,
                length_tolerance=sequence_length_tolerance,
                name=name,
                rule=SegmentRule.REQUIRED,
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

    record = ncbi.fetch_genbank_records([new_accession])[0]

    if record.source.segment:
        segment_id, _ = next(iter(assign_records_to_segments([record], otu.plan)))
    else:
        segment_id = otu.plan.segments[0].id

    with repo.use_transaction():
        new_sequence = repo.create_sequence(
            otu.id,
            accession=record.accession_version,
            definition=record.definition,
            legacy_id=None,
            segment=segment_id,
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

    with repo.use_transaction():
        repo.set_representative_isolate(otu.id, new_representative_isolate.id)

    otu_logger.info(
        f"Representative isolate set to {new_representative_isolate.name}.",
        representative_isolate_id=str(new_representative_isolate.id),
    )

    return new_representative_isolate.id

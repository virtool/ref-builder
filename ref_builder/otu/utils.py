import re
from collections import defaultdict
from enum import StrEnum
from uuid import UUID

import structlog

from ref_builder.models import Molecule
from ref_builder.ncbi.client import NCBIClient
from ref_builder.ncbi.models import NCBIGenbank
from ref_builder.plan import (
    Plan,
    Segment,
    SegmentRule,
    extract_segment_name_from_record,
)
from ref_builder.utils import Accession, IsolateName, IsolateNameType

logger = structlog.get_logger("otu.utils")


class DeleteRationale(StrEnum):
    """Default strings delineating reasons for resource deletion."""

    USER = "Requested by user"
    REFSEQ = "Superceded by RefSeq"


class RefSeqConflictError(ValueError):
    """Raised when a potential RefSeq replacement is found."""

    def __init__(
        self,
        message: str,
        isolate_id: UUID,
        isolate_name: IsolateName,
        accessions: list[str],
    ):
        super().__init__(message)

        self.isolate_id = isolate_id

        self.isolate_name = isolate_name

        self.accessions = accessions


def check_isolate_size(
    plan: Plan,
    isolate_count: int,
) -> bool:
    """Return True if the size of the proposed isolate matches the isolate plan."""
    if plan.monopartite:
        if isolate_count > 1:
            raise ValueError("Too many segments in monopartite isolate.")
        return True

    if isolate_count == len(plan.required_segments):
        return True

    raise ValueError(
        f"The plan requires {len(plan.required_segments)} segments: "
        + f"{[str(segment.name) for segment in plan.segments]}"
    )


def check_sequence_length(sequence: str, segment_length: int, tolerance: float) -> bool:
    """Check if the sequence length is within acceptable segment length tolerance."""
    min_length = segment_length * (1.0 - tolerance)
    max_length = segment_length * (1.0 + tolerance)

    return min_length <= len(sequence) <= max_length


def create_segments_from_records(
    records: list[NCBIGenbank], rule: SegmentRule, length_tolerance: float
) -> list[Segment]:
    """Return a list of SegmentPlans."""
    segments = []

    for record in sorted(records, key=lambda record: record.accession):
        if not record.source.segment:
            raise ValueError("No segment name found for multipartite OTU segment.")

        segments.append(
            Segment.from_record(
                record=record,
                length_tolerance=length_tolerance,
                required=rule,
            ),
        )

    return segments


def create_plan_from_records(
    records: list[NCBIGenbank],
    length_tolerance: float,
    segments: list[Segment] | None = None,
) -> Plan | None:
    """Return a plan from a list of records representing an isolate."""
    if len(records) == 1:
        record = records[0]

        segment_name = None
        try:
            segment_name = extract_segment_name_from_record(record)
        except ValueError:
            pass

        return Plan.new(
            segments=[
                Segment.new(
                    length=len(record.sequence),
                    length_tolerance=length_tolerance,
                    name=segment_name,
                    required=SegmentRule.REQUIRED,
                )
            ]
        )

    binned_records = group_genbank_records_by_isolate(records)
    if len(binned_records) > 1:
        logger.fatal(
            "More than one isolate found. Cannot create schema automatically.",
            bins=binned_records,
        )
        return None

    if segments is not None:
        return Plan.new(segments=segments)

    segments = create_segments_from_records(
        records,
        rule=SegmentRule.REQUIRED,
        length_tolerance=length_tolerance,
    )

    if segments:
        return Plan.new(segments=segments)

    return None


def fetch_records_from_accessions(
    requested_accessions: list | set,
    blocked_accessions: set,
    ignore_cache: bool = False,
) -> list[NCBIGenbank]:
    """Take a list of requested accessions, remove blocked accessions
    and fetch records using the remaining accessions.

    Return a list of records.
    """
    fetch_logger = logger.bind(
        requested=sorted(requested_accessions),
        blocked=sorted(blocked_accessions),
    )

    ncbi = NCBIClient(ignore_cache)

    try:
        fetch_set = set(requested_accessions) - blocked_accessions
        if not fetch_set:
            fetch_logger.error(
                "None of the requested accessions were eligible for inclusion"
            )

            return []

    except ValueError:
        fetch_logger.error(
            "Could not create a new isolate using the requested accessions",
        )
        return []

    fetch_logger.info(
        "Fetching accessions",
        count=len(fetch_set),
        fetch_list=sorted(fetch_set),
    )

    records = ncbi.fetch_genbank_records(fetch_set)

    return records


def get_molecule_from_records(records: list[NCBIGenbank]) -> Molecule:
    """Return relevant molecule metadata from one or more records.
    Molecule metadata is retrieved from the first RefSeq record in the list.
    If no RefSeq record is found in the list, molecule metadata is retrieved
    from record[0].
    """
    if not records:
        raise IndexError("No records given")

    # Assign first record as benchmark to start
    representative_record = records[0]

    if not representative_record.refseq:
        for record in records:
            if record.refseq:
                # Replace representative record with first RefSeq record found
                representative_record = record
                break

    return Molecule.model_validate(
        {
            "strandedness": representative_record.strandedness.value,
            "type": representative_record.moltype.value,
            "topology": representative_record.topology.value,
        },
    )


def group_genbank_records_by_isolate(
    records: list[NCBIGenbank],
) -> dict[IsolateName, dict[Accession, NCBIGenbank]]:
    """Indexes Genbank records by isolate name."""
    isolates = defaultdict(dict)

    for record in records:
        try:
            isolate_name = _get_isolate_name(record)

            if isolate_name is None:
                logger.debug(
                    "RefSeq record does not contain sufficient source data "
                    " for automatic inclusion. Add this record manually.",
                    accession=record.accession,
                    definition=record.definition,
                    source_data=record.source,
                )
                continue

            versioned_accession = Accession.from_string(record.accession_version)
            isolates[isolate_name][versioned_accession] = record

        except ValueError:
            continue

    return isolates


def parse_refseq_comment(comment: str) -> tuple[str, str]:
    """Parse a standard RefSeq comment."""
    if not comment:
        raise ValueError("Empty comment")

    refseq_pattern = re.compile(r"^(\w+ REFSEQ): [\w ]+. [\w ]+ ([\w]+).")

    match = refseq_pattern.search(comment)

    if match:
        return match.group(1), match.group(2)

    raise ValueError("Invalid RefSeq comment")


def _get_isolate_name(record: NCBIGenbank) -> IsolateName | None:
    """Get the isolate name from a Genbank record."""
    if record.source.model_fields_set.intersection(
        {IsolateNameType.ISOLATE, IsolateNameType.STRAIN, IsolateNameType.CLONE},
    ):
        for source_type in IsolateNameType:
            if source_type in record.source.model_fields_set:
                return IsolateName(
                    type=IsolateNameType(source_type),
                    value=record.source.model_dump()[source_type],
                )

    if record.refseq:
        return None

    raise ValueError("Record does not contain sufficient source data for inclusion.")


def assign_records_to_segments(
    records: list[NCBIGenbank], plan: Plan
) -> dict[UUID, NCBIGenbank]:
    """Return a dictionary of records keyed by segment UUID, else raise a ValueError."""
    assigned_records = {}

    unassigned_segments = plan.segments.copy()

    for record in records:
        normalized_segment_name = extract_segment_name_from_record(record)

        for segment in unassigned_segments:
            if segment.name == normalized_segment_name:
                assigned_records[segment.id] = record
                unassigned_segments.remove(segment)
                break

    for segment in plan.required_segments:
        if segment.id not in assigned_records:
            raise ValueError("Missing one or more required segments")

    return assigned_records

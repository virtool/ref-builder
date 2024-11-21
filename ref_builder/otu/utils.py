import re
from collections import defaultdict
from enum import StrEnum
from uuid import UUID

import structlog

from ref_builder.models import Molecule
from ref_builder.ncbi.client import NCBIClient
from ref_builder.ncbi.models import NCBIGenbank
from ref_builder.plan import (
    MonopartitePlan,
    MultipartitePlan,
    Segment,
    SegmentRule,
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
    plan: MonopartitePlan | MultipartitePlan,
    isolate_n: int,
) -> bool:
    """Return True if the size of the proposed isolate matches the isolate plan."""
    if type(plan) is MonopartitePlan:
        if isolate_n > 1:
            raise ValueError("Too many segments in monopartite isolate.")
        return True

    if isolate_n == len(plan.required_segments):
        return True

    raise ValueError(
        f"The plan requires {len(plan.required_segments)} segments: "
        + f"{[str(segment.name) for segment in plan.segments]}"
    )


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


def create_isolate_plan_from_records(
    records: list[NCBIGenbank],
    length_tolerance: float,
    segments: list[Segment] | None = None,
) -> MonopartitePlan | MultipartitePlan | None:
    """Return a plan from a list of records representing an isolate."""
    if len(records) == 1:
        return MonopartitePlan.new(
            length=len(records[0].sequence),
            length_tolerance=length_tolerance,
        )

    binned_records = group_genbank_records_by_isolate(records)
    if len(binned_records) > 1:
        logger.fatal(
            "More than one isolate found. Cannot create schema automatically.",
            bins=binned_records,
        )
        return None

    if segments is not None:
        return MultipartitePlan.new(segments=segments)

    segments = create_segments_from_records(
        records,
        rule=SegmentRule.REQUIRED,
        length_tolerance=length_tolerance,
    )

    if segments:
        return MultipartitePlan.new(segments=segments)

    return None


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


def get_molecule_from_records(records: list[NCBIGenbank]) -> Molecule:
    """Return relevant molecule metadata from one or more records."""
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


def check_sequence_length(sequence: str, segment_length: int, tolerance: float) -> bool:
    if len(sequence) < segment_length * (1.0 - tolerance) or len(
        sequence
    ) > segment_length * (1.0 + tolerance):
        return False

    return True


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


def create_sequence_from_record(
    repo: Repo,
    otu: RepoOTU,
    record: NCBIGenbank,
    segment_name: str | None = None,
) -> RepoSequence | None:
    """Take a NCBI Nucleotide record and create a new sequence."""
    sequence = repo.create_sequence(
        otu.id,
        accession=record.accession_version,
        definition=record.definition,
        legacy_id=None,
        segment=segment_name if segment_name is not None else record.source.segment,
        sequence=record.sequence,
    )

    return sequence


def fetch_records_from_accessions(
    requested_accessions: list | set,
    blocked_accessions: set,
    ignore_cache: bool = False,
) -> list[NCBIGenbank]:
    fetch_logger = logger.bind()

    ncbi = NCBIClient(ignore_cache)

    try:
        fetch_list = list(set(requested_accessions).difference(blocked_accessions))
        if not fetch_list:
            raise ValueError(
                "None of the requested accessions were eligible for inclusion"
            )
    except ValueError:
        fetch_logger.error(
            "Could not create a new isolate using the requested accessions.",
            requested=sorted(requested_accessions),
            blocked=sorted(blocked_accessions),
        )
        return []

    fetch_logger.info(
        "Fetching accessions",
        count=len(fetch_list),
        fetch_list=fetch_list,
    )

    records = ncbi.fetch_genbank_records(fetch_list)

    return records

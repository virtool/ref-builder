import re
from contextlib import suppress
from enum import StrEnum
from typing import Optional
from uuid import UUID, uuid4

from pydantic import UUID4, BaseModel, ConfigDict, Field, model_validator
from pydantic.dataclasses import dataclass

from ref_builder.ncbi.models import NCBIGenbank

SIMPLE_NAME_PATTERN = re.compile(r"([A-Za-z0-9])+")
"""Regex pattern for parsing segment name strings with no prefix."""

COMPLEX_NAME_PATTERN = re.compile(r"([A-Za-z]+)[-_ ]+(.*)")
"""Regex pattern for parsing segment name strings consisting of a prefix and a key."""


class PlanWarning(UserWarning):
    pass


class PlanConformationError(ValueError):
    """Raised when potential new sequences do not pass validation against the OTU plan."""


class SegmentRule(StrEnum):
    """Mark the importance of a particular segment."""

    REQUIRED = "required"
    """Segment is always required."""

    RECOMMENDED = "recommended"
    """Risk has to be acknowledged to add violating isolate."""

    OPTIONAL = "optional"
    """Segment is entirely optional."""


@dataclass(frozen=True)
class SegmentName:
    """A normalized segment name. Can be used as a key."""

    prefix: str
    """The prefix of the segment name."""

    key: str
    """The identifying key portion of the segment name."""

    def __str__(self) -> str:
        """Return the segment name as a formatted string."""
        return f"{self.prefix} {self.key}"

    @classmethod
    def from_string(cls: type, string: str) -> Optional["SegmentName"]:
        """Create a SegmentName from a string.

        Return None if the string does not match the expected format.
        """
        segment_name_parse = COMPLEX_NAME_PATTERN.fullmatch(string)

        if segment_name_parse:
            return SegmentName(
                prefix=segment_name_parse.group(1),
                key=segment_name_parse.group(2),
            )

        return None


class Segment(BaseModel):
    """A segment in a multipartite plan."""

    model_config = ConfigDict(use_enum_values=True, validate_assignment=True)

    id: UUID4
    """The unique ID of the segment or monopartite plan."""

    length: int
    """The expected length of the sequence"""

    length_tolerance: float = Field(ge=0.0, le=1.0)
    """The acceptable deviation from the recommended sequence length."""

    name: SegmentName | None
    """The name of the segment"""

    rule: SegmentRule
    """Whether this segment must be present in all additions."""

    @classmethod
    def new(
        cls,
        length: int,
        length_tolerance: float,
        name: SegmentName | None,
        rule: SegmentRule = SegmentRule.REQUIRED,
    ) -> "Segment":
        """Return a new segment."""
        return Segment(
            id=uuid4(),
            length=length,
            length_tolerance=length_tolerance,
            name=name,
            rule=rule,
        )

    @classmethod
    def from_record(
        cls,
        record: NCBIGenbank,
        length_tolerance: float,
        rule: SegmentRule,
    ) -> "Segment":
        """Return a new segment from an NCBI Genbank record."""
        return Segment(
            id=uuid4(),
            length=len(record.sequence),
            length_tolerance=length_tolerance,
            name=extract_segment_name_from_record(record),
            rule=rule,
        )


class Plan(BaseModel):
    """The segments required for an isolate in a given OTU."""

    id: UUID4
    """The unique id number of the multipartite plan"""

    segments: list[Segment]
    """A list of segments that define the plan."""

    @property
    def segment_ids(self) -> set[UUID]:
        """Return all segment IDs as a set."""
        return {segment.id for segment in self.segments}

    @property
    def monopartite(self) -> bool:
        """Whether the plan is monopartite."""
        return len(self.segments) == 1

    @property
    def required_segments(self) -> list[Segment]:
        """Return a list of segments that are required by all additions."""
        return [
            segment for segment in self.segments if segment.rule == SegmentRule.REQUIRED
        ]

    @property
    def required_segment_ids(self) -> set[UUID]:
        """Return a set of segments that are required by all additions."""
        return {
            segment.id
            for segment in self.segments
            if segment.rule == SegmentRule.REQUIRED
        }

    @property
    def not_required_segments(self) -> list[Segment]:
        """Return a list of segments that are not always required for inclusion."""
        return [
            segment for segment in self.segments if segment.rule != SegmentRule.REQUIRED
        ]

    @classmethod
    def new(cls, segments: list[Segment]) -> "Plan":
        """Initialize a new Plan from a list of segments."""
        return Plan(id=uuid4(), segments=segments)

    @model_validator(mode="after")
    def check_naming(self) -> "Plan":
        """Check that all segments have non-None, unique names if the plan is
        multipartite.
        """
        if self.monopartite:
            return self

        names = [segment.name for segment in self.segments]

        if any(name is None for name in names):
            raise ValueError("All segments must have a name in a multipartite plan.")

        if len(names) != len(set(names)):
            raise ValueError("Segment names must be unique within a plan.")

        return self

    def get_segment_by_id(self, segment_id: UUID) -> Segment | None:
        """Get the segment with the given ``segment_id`.

        Return ``None`` if the segment is not found.

        :param segment_id: the ID of the segment to retrieve
        :return: the segment with the given ID, or None
        """
        for segment in self.segments:
            if segment.id == segment_id:
                return segment

        return None

    def get_segment_by_name_key(self, name_key: str) -> Segment | None:
        """Return the segment with the given ``name.key`` if it exists."""
        for segment in self.segments:
            if segment.name is None:
                continue

            if segment.name.key == name_key:
                return segment

        return None


def extract_segment_name_from_record(record: NCBIGenbank) -> SegmentName | None:
    """Extract a segment name from a Genbank record.

    This parsing function should be used when creating new OTUs where no plan is
    available. Parsing segment names for existing OTUs should use plans to better infer
    segment names.

    First, try to parse the segment name from the ``segment`` field in the source table.
    If that fails, use the record moltype as a prefix and the segment as a key.

    If no segment name can be extracted, `None` is returned.

    :param record: A Genbank record.
    :return: A segment name or `None`.
    """
    if (name := SegmentName.from_string(record.source.segment)) is not None:
        return name

    # Handles common cases without delimiters
    for moltype_prefix in ["DNA", "RNA"]:
        if record.source.segment.startswith(moltype_prefix):
            return SegmentName(
                prefix=record.moltype,
                key=record.source.segment[3:].strip(),
            )

    if SIMPLE_NAME_PATTERN.fullmatch(record.source.segment):
        return SegmentName(
            prefix=record.moltype,
            key=record.source.segment,
        )

    return None


def extract_segment_name_from_record_with_plan(
    record: NCBIGenbank, plan: Plan
) -> SegmentName | None:
    """Extract a segment name from a Genbank record using an OTU plan.

    If the segment name can be parsed normally (eg. DNA A), it is returned. If there is
    no delimiter or the prefix is missing, the function tries to match the segment name
    with a segment in the plan.

    If no segment name can be extracted, `None` is returned.

    :param record: A Genbank record.
    :param plan: A plan.
    :return: A segment name or `None`.
    """
    if not record.source.segment:
        return None

    if (segment_name := SegmentName.from_string(record.source.segment)) is not None:
        return segment_name

    if not plan.monopartite:
        try:
            plan_keys_and_prefixes = {
                segment.name.key: segment.name.prefix for segment in plan.segments
            }
        except AttributeError:
            raise ValueError("Multipartite plan contains unnamed segments")

        # Handle no delimiter.
        for prefix in plan_keys_and_prefixes.values():
            if record.source.segment.casefold().startswith(prefix.casefold()):
                return SegmentName(
                    prefix=prefix, key=record.source.segment[len(prefix) :].strip()
                )

        # Handle no prefix.
        with suppress(KeyError):
            return SegmentName(
                prefix=plan_keys_and_prefixes[record.source.segment],
                key=record.source.segment,
            )

    return None

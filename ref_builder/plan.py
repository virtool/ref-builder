import re
from enum import StrEnum
from uuid import UUID, uuid4

from pydantic import UUID4, BaseModel, ConfigDict, Field, model_validator
from pydantic.dataclasses import dataclass

from ref_builder.ncbi.models import NCBIGenbank

SIMPLE_NAME_PATTERN = re.compile(r"([A-Za-z0-9])+")
"""Regex pattern for parsing segment name strings with no prefix."""

COMPLEX_NAME_PATTERN = re.compile(r"([A-Za-z]+)[-_ ]+([A-Za-z0-9]+)")
"""Regex pattern for parsing segment name strings consisting of a prefix and a key."""


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


class Segment(BaseModel):
    """A segment in a multipartite plan."""

    model_config = ConfigDict(use_enum_values=True)

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


def parse_segment_name(string: str) -> SegmentName | None:
    """Parse a SegmentName from a string.

    Returns none if no segment name is found in the string.
    """
    segment_name_parse = COMPLEX_NAME_PATTERN.fullmatch(string)

    if segment_name_parse:
        return SegmentName(
            prefix=segment_name_parse.group(1),
            key=segment_name_parse.group(2),
        )

    return None


def extract_segment_name_from_record(record: NCBIGenbank) -> SegmentName | None:
    """Get a segment name from a Genbank record."""
    if SIMPLE_NAME_PATTERN.fullmatch(record.source.segment):
        return SegmentName(
            prefix=record.moltype,
            key=record.source.segment,
        )

    return parse_segment_name(record.source.segment)

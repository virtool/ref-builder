import re
from enum import StrEnum
from uuid import UUID, uuid4

from pydantic import UUID4, BaseModel, ConfigDict, Field
from pydantic.dataclasses import dataclass

from ref_builder.ncbi.models import NCBIGenbank, NCBISourceMolType

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


class SegmentMetadata(BaseModel):
    """Metadata and expected properties for an included sequence."""

    id: UUID4
    """The unique ID of the segment or monopartite plan."""

    length: int
    """The expected length of the sequence"""


class Segment(SegmentMetadata):
    """A segment in a multipartite plan."""

    model_config = ConfigDict(use_enum_values=True)

    length_tolerance: float = Field(ge=0.0, le=1.0)
    """The acceptable deviation from the recommended sequence length."""

    name: SegmentName | None
    """The name of the segment"""

    required: SegmentRule
    """Whether this segment must be present in all additions."""

    @classmethod
    def new(
        cls,
        length: int,
        length_tolerance: float,
        name: SegmentName | None,
        required: SegmentRule,
    ) -> "Segment":
        """Return a new segment."""
        return Segment(
            id=uuid4(),
            length=length,
            length_tolerance=length_tolerance,
            name=name,
            required=required,
        )

    @classmethod
    def from_record(
        cls,
        record: NCBIGenbank,
        length_tolerance: float,
        required: SegmentRule,
    ) -> "Segment":
        """Return a new segment from an NCBI Genbank record."""
        return Segment(
            id=uuid4(),
            length=len(record.sequence),
            length_tolerance=length_tolerance,
            name=get_multipartite_segment_name(record),
            required=required,
        )


class Plan(BaseModel):
    """The segments required for an isolate in a given OTU."""

    id: UUID4
    """The unique id number of the multipartite plan"""

    segments: list[Segment]
    """A list of segments"""

    @property
    def required_segments(self) -> list[Segment]:
        """Return a list of segments that are required by all additions."""
        return [
            segment
            for segment in self.segments
            if segment.required == SegmentRule.REQUIRED
        ]

    @property
    def not_required_segments(self) -> list[Segment]:
        """Return a list of segments that are not always required for inclusion."""
        return [
            segment
            for segment in self.segments
            if not segment.required == SegmentRule.REQUIRED
        ]

    @classmethod
    def new(cls, segments: list["Segment"]) -> "Plan":
        """Initialize a new Plan from a list of segments."""
        return Plan(id=uuid4(), segments=segments)

    def get_segment_by_id(self, segment_id: UUID) -> Segment | None:
        """Return the segment with the matching segment ID if it exists,
        otherwise return None.
        """
        for segment in self.segments:
            if segment.id == segment_id:
                return segment

        return None

    def get_segment_by_key(self, name_key: str) -> Segment | None:
        """Return the segment with the matching SegmentName.key it exists,
        otherwise return None.
        """
        for segment in self.segments:
            if segment.name.key == name_key:
                return segment

        return None


def determine_segment_prefix(moltype: NCBISourceMolType) -> str:
    """Return an acceptable SegmentName prefix corresponding to
    the given NCBISourceMolType.
    """
    if moltype in (
        NCBISourceMolType.GENOMIC_DNA,
        NCBISourceMolType.OTHER_DNA,
        NCBISourceMolType.UNASSIGNED_DNA,
    ):
        return "DNA"

    prefix = None

    match moltype:
        case NCBISourceMolType.GENOMIC_RNA:
            prefix = "RNA"
        case NCBISourceMolType.MRNA:
            prefix = "mRNA"
        case NCBISourceMolType.TRNA:
            prefix = "tRNA"
        case NCBISourceMolType.TRANSCRIBED_RNA:
            prefix = "RNA"
        case NCBISourceMolType.VIRAL_CRNA:
            prefix = "cRNA"
        case NCBISourceMolType.OTHER_RNA:
            prefix = "RNA"

    if prefix:
        return prefix

    raise ValueError(f"{moltype} may not be a valid NCBISourceMolType.")


def parse_segment_name(raw: str) -> SegmentName:
    """Parse a SegmentName from a raw string."""
    segment_name_parse = COMPLEX_NAME_PATTERN.fullmatch(raw)
    if segment_name_parse:
        return SegmentName(
            prefix=segment_name_parse.group(1),
            key=segment_name_parse.group(2),
        )

    raise ValueError(f"{raw} is not a valid segment name")


def get_multipartite_segment_name(record: NCBIGenbank) -> SegmentName:
    """Get a multipartite segment name from the record."""
    if SIMPLE_NAME_PATTERN.fullmatch(record.source.segment):
        return SegmentName(
            prefix=record.moltype,
            key=record.source.segment,
        )

    return parse_segment_name(record.source.segment)

import re
from enum import StrEnum
from typing import Annotated, Literal, Union
from uuid import uuid4

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


class MonopartitePlan(BaseModel):
    """Expected properties for an acceptable monopartite isolate."""

    plan_type: Literal["monopartite"]
    """The plan type identifier."""

    id: UUID4
    """The unique ID of the monopartite plan."""

    length: int
    """The expected length of the sequence"""

    length_tolerance: float = Field(ge=0.0, le=1.0)
    """The acceptable deviation from the recommended sequence length."""

    name: SegmentName | None = None
    """The name of the monopartite plan"""

    @property
    def segments(self) -> list[Segment]:
        """Return a simulated single-segment list of segments."""
        return [
            Segment(
                id=self.id,
                length=self.length,
                length_tolerance=self.length_tolerance,
                name=self.name,
                required=SegmentRule.REQUIRED,
            )
        ]

    @property
    def required_segments(self) -> list[Segment]:
        """Return a simulated single-segment list of segments."""
        return self.segments

    @classmethod
    def new(
        cls,
        length: int,
        length_tolerance: float,
        name: SegmentName | None = None,
    ) -> "MonopartitePlan":
        """Initialize a MonopartitePlan from a list of segments."""
        return MonopartitePlan(
            id=uuid4(),
            plan_type="monopartite",
            length_tolerance=length_tolerance,
            length=length,
            name=name,
        )


class MultipartitePlan(BaseModel):
    """Expected segments for an acceptable multipartite isolate."""

    plan_type: Literal["multipartite"]
    """The plan type identifier."""

    id: UUID4
    """The unique id number of the multipartite plan"""

    segments: list[Segment]

    @property
    def required_segments(self) -> list[Segment]:
        """Return a list of segments that are required by all additions."""
        return [segment for segment in self.segments if segment.required]

    @classmethod
    def new(cls, segments: list["Segment"]) -> "MultipartitePlan":
        """Initialize a MultipartitePlan from a list of segments."""
        return MultipartitePlan(
            id=uuid4(),
            plan_type="multipartite",
            segments=segments,
        )


Plan = Annotated[
    Union[MultipartitePlan, MonopartitePlan], Field(discriminator="plan_type")
]
"""A discriminated union representing both plan types, monopartite and multipartite."""


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

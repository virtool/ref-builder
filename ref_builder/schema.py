import re
from typing import NamedTuple

from pydantic import BaseModel, UUID4, computed_field

from ref_builder.models import Molecule
from ref_builder.ncbi.models import NCBISourceMolType


class SegmentName(NamedTuple):
    prefix: str
    key: str


class Segment(BaseModel):
    """The metadata of the segment"""
    id: UUID4
    """The unique id number of this segment"""

    name: str
    """The name of the segment"""

    required: bool
    """Whether this segment must be present in all additions."""

    length: int | None = None


class OTUSchema(BaseModel):
    """A schema for the intended data"""

    molecule: Molecule
    """The molecular metadata for this OTU."""

    segments: list[Segment]
    """The segments contained in this OTU."""

    @computed_field
    def multipartite(self) -> bool:
        return len(self.segments) > 1


def set_segment_prefix(moltype: NCBISourceMolType):
    if moltype in (
        NCBISourceMolType.GENOMIC_DNA,
        NCBISourceMolType.OTHER_DNA,
        NCBISourceMolType.UNASSIGNED_DNA,
    ):
        return "DNA"

    match moltype:
        case NCBISourceMolType.GENOMIC_RNA:
            return "RNA"
        case NCBISourceMolType.MRNA:
            return "mRNA"
        case NCBISourceMolType.TRNA:
            return "tRNA"
        case NCBISourceMolType.TRANSCRIBED_RNA:
            return "RNA"
        case NCBISourceMolType.VIRAL_CRNA:
            return "cRNA"
        case NCBISourceMolType.OTHER_RNA:
            return "RNA"


simple_name_pattern = re.compile(r"([A-Za-z0-9])+")

complex_name_pattern = re.compile(r"([A-Za-z]+)[-_ ]+([A-Za-z0-9]+)")


def parse_segment_name(raw: str) -> SegmentName:
    segment_name_parse = complex_name_pattern.fullmatch(raw)
    if segment_name_parse:
        return SegmentName(
            prefix=segment_name_parse.group(1),
            key=segment_name_parse.group(2)
        )

    raise ValueError(f"{raw} is not a valid segment name")

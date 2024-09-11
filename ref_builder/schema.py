from typing import NamedTuple

from pydantic import BaseModel, UUID4, computed_field

from ref_builder.models import Molecule
from ref_builder.ncbi.models import NCBISourceMolType

"""Regex pattern for parsing segment name strings containing both a prefix and a key."""


class SegmentName(NamedTuple):
    prefix: str
    key: str

    def __str__(self) -> str:
        """Return the segment name as a formatted string."""
        return f"{self.prefix} {self.key}"


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


def determine_segment_prefix(moltype: NCBISourceMolType) -> str:
    """Returns an acceptable SegmentName prefix
    corresponding to the given NCBISourceMolType.
    """
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

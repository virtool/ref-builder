from pydantic import BaseModel, computed_field

from ref_builder.models import Molecule


class Segment(BaseModel):
    """The metadata of the segment"""

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

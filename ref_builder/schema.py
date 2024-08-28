from pydantic import BaseModel, computed_field, model_validator

from ref_builder.models import Molecule


class Segment(BaseModel):
    """The metadata of the segment"""

    id: int
    """The numerical identifier of the segment."""

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

    @model_validator(mode="after")
    def check_segment_ids(self):
        last_segment_id = 0
        for segment in self.segments:
            if segment.id != last_segment_id + 1:
                raise ValueError("Segment IDs are not sequential.")
            last_segment_id += 1

        return self



import re

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


simple_name_pattern = re.compile(r"([A-Za-z0-9])+")

complex_name_pattern = re.compile(r"([A-Za-z]+)[-_ ]+([A-Za-z0-9]+)")


def parse_segment_name(raw: str):
    """Takes a raw segment name from a NCBI Genbank source table
    and standardizes the relevant identifier"""
    if simple_name_pattern.fullmatch(raw):
        return raw

    segment_name_parse = complex_name_pattern.fullmatch(raw)
    if segment_name_parse:
        return segment_name_parse.group(2)

    raise ValueError(f"{raw} is not a valid segment name")



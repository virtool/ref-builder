from dataclasses import dataclass
from enum import Enum
from typing import Annotated

from pydantic import (
    BaseModel,
    Field,
    StringConstraints,
    field_validator,
    model_validator,
)


def transform_to_uppercase(string: str) -> str:
    return string.upper()


class LegacySourceType(str, Enum):
    CLONE = "clone"
    CULTURE = "culture"
    ISOLATE = "isolate"
    GENBANK = "genbank"
    GENOTYPE = "genotype"
    SEROTYPE = "serotype"
    STRAIN = "strain"
    UNKNOWN = "unknown"
    VARIANT = "variant"


@dataclass
class LegacyIsolateSource:
    name: str
    type: LegacySourceType


class LegacySequence(BaseModel):
    id: str = Field(alias="_id")
    accession: Annotated[
        str,
        StringConstraints(
            min_length=4,
            pattern=r"^\w+\.\d+$",
            strip_whitespace=True,
            to_upper=True,
        ),
    ]
    definition: Annotated[str, StringConstraints(min_length=10, strip_whitespace=True)]
    host: Annotated[str, StringConstraints(strip_whitespace=True)]
    segment: str
    sequence: Annotated[
        str,
        StringConstraints(
            min_length=10,
            pattern="^[ATCGRYKMSWBDHVNatcgrykmswbdhvn]+$",
            strip_whitespace=True,
            to_upper=True,
        ),
    ]


class LegacyIsolate(BaseModel):
    id: str
    default: bool
    sequences: Annotated[list[LegacySequence], StringConstraints(min_length=1)]
    source_type: LegacySourceType
    source_name: Annotated[str, StringConstraints(strip_whitespace=True)]

    @model_validator(mode="after")
    def check_source_name(self: "LegacyIsolate") -> "LegacyIsolate":
        """Check if the source name is not empty."""
        if not self.source_name and self.source_type is not LegacySourceType.UNKNOWN:
            raise ValueError(
                "Source name cannot be empty unless source type is unknown"
            )

        return self


class LegacySchemaSegment(BaseModel):
    molecule: str
    name: str
    required: bool


class LegacyOTU(BaseModel):
    id: Annotated[str, Field(alias="_id")]
    abbreviation: Annotated[str, StringConstraints(strip_whitespace=True)]
    isolates: Annotated[list[LegacyIsolate], Field(min_length=1)]
    name: Annotated[str, StringConstraints(min_length=5, strip_whitespace=True)]
    otu_schema: Annotated[
        list[LegacySchemaSegment],
        Field(min_length=1, alias="schema"),
    ]
    taxid: int

    @field_validator("isolates")
    @classmethod
    def check_default_isolate(cls, v):
        """Check if there is only one default isolate."""
        default_isolates = [isolate for isolate in v if isolate.default]

        if len(default_isolates) == 0:
            raise ValueError("At least one isolate must be default")

        if len(default_isolates) > 1:
            raise ValueError("Only one isolate can be default")

        return v

    @field_validator("otu_schema")
    @classmethod
    def check_schema_molecule(cls, otu_schema: list[LegacySchemaSegment]):
        """Check if all segments in the schema have the same molecule."""
        molecules = {segment.molecule for segment in otu_schema}

        if len(molecules) > 1:
            raise ValueError("All segments in a schema must have the same molecule")

        return otu_schema

    @field_validator("otu_schema")
    @classmethod
    def check_schema_names(cls, otu_schema: list[LegacySchemaSegment]):
        """Check if there are duplicate schema segment names."""
        names = {segment.name for segment in otu_schema}

        if len(names) != len(otu_schema):
            raise ValueError("All schema segments must have a unique name")

        return otu_schema

    @model_validator(mode="after")
    def check_schema_and_sequences(self) -> "LegacyOTU":
        """Make sure all sequences in the OTU have a corresponding schema segment."""
        schema_segments = {segment.name for segment in self.otu_schema}

        required_schema_segments = {
            segment.name for segment in self.otu_schema if segment.required
        }

        for isolate in self.isolates:
            for sequence in isolate.sequences:
                if sequence.segment not in schema_segments:
                    raise ValueError(
                        "sequence contains invalid segment name",
                    )

            if required_schema_segments - {
                sequence.segment for sequence in isolate.sequences
            }:
                raise ValueError(
                    "isolate does not contain all required schema segments",
                )

        sequence_segments = {
            sequence.segment
            for isolate in self.isolates
            for sequence in isolate.sequences
        }

        if schema_segments != sequence_segments:
            raise ValueError("All schema segments must have a corresponding sequence")

        return self

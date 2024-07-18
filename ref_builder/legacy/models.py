"""Models for the legacy data format used by Virtool and in virtool-cli."""

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


class LegacySourceType(str, Enum):
    """Possible source types for legacy isolates."""

    CLONE = "clone"
    CULTURE = "culture"
    ISOLATE = "isolate"
    GENBANK = "genbank"
    GENOTYPE = "genotype"
    REFSEQ = "refseq"
    SEROTYPE = "serotype"
    STRAIN = "strain"
    UNKNOWN = "unknown"
    VARIANT = "variant"


@dataclass
class LegacyIsolateSource:
    """The source of a legacy isolate."""

    name: str
    type: LegacySourceType


class LegacySequence(BaseModel):
    """A sequence as stored in a legacy isolate."""

    id: str = Field(alias="_id")
    """The unique identifier for the sequence.

    Originated from the ``_id`` stored in MongoDB by Virtool. They are random
    eight-character alphanumeric strings.
    """

    accession: Annotated[
        str,
        StringConstraints(
            min_length=4,
            pattern=r"^\w+\.\d+$",
            strip_whitespace=True,
            to_upper=True,
        ),
    ]
    """The NCBI accession for the sequence.

    Inclusion of the version number was not required.
    """

    definition: Annotated[str, StringConstraints(min_length=10, strip_whitespace=True)]
    """The NCBI definition for the sequence."""

    host: Annotated[str, StringConstraints(strip_whitespace=True)]
    """The host from the sequence source table."""

    segment: str
    """The name of the segment the sequence belongs to.

    For example, "RNA A". This segment name would be found in the OTU schema field.
    """

    sequence: Annotated[
        str,
        StringConstraints(
            min_length=10,
            pattern="^[ATCGRYKMSWBDHVNatcgrykmswbdhvn]+$",
            strip_whitespace=True,
            to_upper=True,
        ),
    ]
    """The nucleotide sequence."""


class LegacyIsolate(BaseModel):
    """An isolate as stored in a legacy OTU."""

    id: str
    """The unique identifier for the isolate.

    These are random eight-character alphanumeric strings.
    """

    default: bool
    """Whether the isolate is the default isolate for the OTU."""

    sequences: Annotated[list[LegacySequence], StringConstraints(min_length=1)]
    """The sequences belonging to the isolate."""

    source_type: LegacySourceType
    """The source type for the isolate.

    For example, "isolate" in the isolate name "Isolate A". The source type is always
    lowercase.
    """

    source_name: Annotated[str, StringConstraints(min_length=1, strip_whitespace=True)]
    """The source name for the isolate.

    For example, "A" in the isolate name "Isolate A".
    """


class LegacySchemaSegment(BaseModel):
    """A schema segment as stored in the ``schema`` field of a legacy OTU."""

    molecule: str
    """The type of molecule the segment is (eg. ssRNA)."""

    name: str
    """The name for the segment (eg. RNA A)."""

    required: bool
    """Whether a sequence assigned as this segment is required in all isolates."""


class LegacyOTU(BaseModel):
    """An OTU as stored in the legacy data format."""

    id: Annotated[str, Field(alias="_id")]
    """The unique identifier for the OTU.

    Originated from the ``_id`` stored in MongoDB by Virtool. They are random
    eight-character alphanumeric strings.
    """

    abbreviation: Annotated[str, StringConstraints(strip_whitespace=True)]
    """The abbreviation for the OTU."""

    isolates: Annotated[list[LegacyIsolate], Field(min_length=1)]
    """The isolates belonging to the OTU."""

    name: Annotated[str, StringConstraints(min_length=5, strip_whitespace=True)]
    """The name of the OTU."""

    otu_schema: Annotated[
        list[LegacySchemaSegment],
        Field(min_length=1, alias="schema"),
    ]
    """The schema for the OTU."""

    taxid: int
    """The NCBI taxonomy ID for the OTU."""

    @classmethod
    @field_validator("isolates")
    def check_default_isolate(cls, v: list[LegacyIsolate]) -> list[LegacyIsolate]:
        """Check if there is only one default isolate."""
        default_isolates = [isolate for isolate in v if isolate.default]

        if len(default_isolates) == 0:
            raise ValueError("At least one isolate must be default")

        if len(default_isolates) > 1:
            raise ValueError("Only one isolate can be default")

        return v

    @classmethod
    @field_validator("otu_schema")
    def check_schema_molecule(
        cls,
        v: list[LegacySchemaSegment],
    ) -> list[LegacySchemaSegment]:
        """Check if all segments in the schema have the same molecule."""
        molecules = {segment.molecule for segment in v}

        if len(molecules) > 1:
            raise ValueError("All segments in a schema must have the same molecule")

        return v

    @classmethod
    @field_validator("otu_schema")
    def check_schema_name(
        cls,
        v: list[LegacySchemaSegment],
    ) -> list[LegacySchemaSegment]:
        """Check if there are duplicate schema segment names."""
        if len({segment.name for segment in v}) != len(v):
            raise ValueError("All schema segments must have a unique name")

        return v

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

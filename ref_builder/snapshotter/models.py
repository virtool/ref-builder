from typing import Annotated

from pydantic import (
    UUID4,
    BaseModel,
    Field,
    TypeAdapter,
    field_validator,
)

from ref_builder.schema import OTUSchema
from ref_builder.utils import IsolateName, IsolateNameType


class OTUSnapshotMeta(BaseModel):
    """Structures metadata about the OTU snapshot itself."""

    at_event: int | None = None
    """The event ID of the last change made to this snapshot."""


class OTUSnapshotSequence(BaseModel):
    """Stores and parses sequence data."""

    id: UUID4
    """The sequence id."""

    accession: str
    """The sequence accession."""

    definition: str
    """The sequence definition."""

    legacy_id: str | None
    """A string based ID carried over from a legacy Virtool reference repository.

    It the sequence was not migrated from a legacy repository, this will be `None`.
    """

    sequence: str
    """The sequence."""

    segment: str
    """The sequence segment."""


class OTUSnapshotIsolate(BaseModel):
    """Stores and parses isolate data."""

    id: UUID4
    """The isolate ID."""

    name: IsolateName
    """The isolate's source name metadata."""

    legacy_id: str | None = None
    """A string based ID carried over from a legacy Virtool reference repository."""

    @field_validator("name", mode="before")
    @classmethod
    def convert_isolate_name(cls, raw: dict) -> IsolateName:
        """Takes a dictionary and converts to IsolateName."""
        return IsolateName(type=IsolateNameType(raw["type"]), value=raw["value"])


class OTUSnapshotOTU(BaseModel):
    """Stores and parses OTU data."""

    id: UUID4
    """The OTU id."""

    taxid: int
    """The NCBI Taxonomy id for this OTU."""

    name: str
    """The name of the OTU (eg. TMV for Tobacco mosaic virus)"""

    acronym: str = ""
    """The OTU acronym (eg. TMV for Tobacco mosaic virus)."""

    legacy_id: str | None
    """A string based ID carried over from a legacy Virtool reference repository."""

    otu_schema: Annotated[OTUSchema | None, Field(alias="schema")] = None

    repr_isolate: UUID4 | None = None


class OTUSnapshotToCIsolate(BaseModel):
    """Stores a table of contents for an isolate."""

    id: UUID4
    """The isolate id."""

    accessions: dict[str, UUID4]
    """A mapping of accessions to sequence ids."""


toc_adapter = TypeAdapter(dict[str, OTUSnapshotToCIsolate])

from typing import Annotated

from pydantic import (
    UUID4,
    AliasChoices,
    BaseModel,
    Field,
    TypeAdapter,
    field_serializer,
    field_validator,
)

from ref_builder.resources import OTUModel
from ref_builder.schema import OTUSchema
from ref_builder.utils import Accession, IsolateName, IsolateNameType


class OTUSnapshotToCIsolate(BaseModel):
    """Stores a table of contents for an isolate."""

    id: UUID4
    """The isolate id."""

    name: IsolateName | None
    """The isolate name."""

    accessions: dict[str, UUID4]
    """A mapping of accessions to sequence ids."""


toc_adapter = TypeAdapter(dict[str, OTUSnapshotToCIsolate])

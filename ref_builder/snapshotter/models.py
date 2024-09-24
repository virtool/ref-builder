from pydantic import BaseModel, TypeAdapter, UUID4

from ref_builder.utils import IsolateName


class OTUSnapshotToCIsolate(BaseModel):
    """Stores a table of contents for an isolate."""

    id: UUID4
    """The isolate id."""

    name: IsolateName | None
    """The isolate name."""

    accessions: dict[str, UUID4]
    """A mapping of accessions to sequence ids."""


toc_adapter = TypeAdapter(dict[str, OTUSnapshotToCIsolate])

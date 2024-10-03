from dataclasses import dataclass
from enum import StrEnum
from uuid import UUID

from pydantic import BaseModel, ConfigDict


@dataclass
class OTUMinimal:
    """A minimal representation of an OTU."""

    id: UUID
    acronym: str
    legacy_id: str | None
    name: str
    taxid: int


class MolType(StrEnum):
    """The in vivo molecule type of a sequence.

    Corresponds to Genbank's moltype field
    """

    CRNA = "cRNA"
    DNA = "DNA"
    MRNA = "mRNA"
    RNA = "RNA"
    TRNA = "tRNA"


class Strandedness(StrEnum):
    """Strandedness of a molecule, either single or double."""

    SINGLE = "single"
    DOUBLE = "double"


class Topology(StrEnum):
    """Topology of a molecule, either linear or circular."""

    LINEAR = "linear"
    CIRCULAR = "circular"


class Molecule(BaseModel):
    """The strandedness, molecule type, and topology of an OTU."""

    model_config = ConfigDict(use_enum_values=True)

    strandedness: Strandedness
    topology: Topology
    type: MolType

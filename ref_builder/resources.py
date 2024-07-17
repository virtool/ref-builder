import datetime
from dataclasses import dataclass
from typing import Any
from uuid import UUID

from pydantic import BaseModel

from ref_builder.schema import OTUSchema
from ref_builder.utils import DataType, IsolateName


class RepoMeta(BaseModel):
    """Represents the metadata for a Virtool reference repository."""

    id: UUID
    """The repository id."""

    created_at: datetime.datetime
    """The date and time the repository was created."""

    data_type: DataType
    """The repository data type."""

    name: str
    """The repository name."""

    organism: str
    """The repository organism."""


@dataclass
class RepoSequence:
    """Represents a sequence in a Virtool reference repository."""

    id: UUID
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

    def dict(self) -> dict[str, Any]:
        """Return a dictionary representation of the sequence."""
        return {
            "id": self.id,
            "accession": self.accession,
            "definition": self.definition,
            "legacy_id": self.legacy_id,
            "sequence": self.sequence,
            "segment": self.segment,
        }


class RepoIsolate:
    """Represents an isolate in a Virtool reference repository."""

    def __init__(
        self,
        uuid: UUID,
        name: IsolateName,
        sequences: list[RepoSequence] | None = None,
        legacy_id: str | None = None,
    ) -> None:
        """Initialize a new isolate."""
        self.id = uuid
        """The isolate id."""

        self.name = name
        """The isolate's source name metadata."""

        self._sequences_by_accession = (
            {}
            if sequences is None
            else {sequence.accession: sequence for sequence in sequences}
        )
        """A dictionary of sequences indexed by accession"""

        self.legacy_id = legacy_id
        """A string based ID carried over from a legacy Virtool reference repository.

        It the isolate was not migrated from a legacy repository, this will be `None`.
        """

        self.__hash__ = None

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> "RepoIsolate":
        """Build a new isolate from .dict() output."""
        return RepoIsolate(
            uuid=data["id"],
            name=IsolateName(type=data["name"]["type"], value=data["name"]["value"]),
            sequences=data.get("sequences"),
            legacy_id=data.get("legacy_id"),
        )

    @property
    def sequences(self) -> list[RepoSequence]:
        """A list of sequences in this isolate."""
        return list(self._sequences_by_accession.values())

    @property
    def accessions(self) -> set[str]:
        """A set of accession numbers for sequences in the isolate."""
        return set(self._sequences_by_accession.keys())

    @property
    def sequence_ids(self) -> set[UUID]:
        """A set of UUIDs for sequences in the isolate."""
        return {sequence.id for sequence in self.sequences}

    def __repr__(self) -> str:
        """Return a string representation of the isolate."""
        return (
            "EventSourcedRepoIsolate("
            f"{self.id}, type={self.name.type}, name={self.name.value}, "
            f"accessions={self.accessions})"
        )

    def add_sequence(self, sequence: RepoSequence) -> None:
        """Add a sequence to the isolate."""
        self._sequences_by_accession[sequence.accession] = sequence

    def get_sequence_by_accession(
        self,
        accession: str,
    ) -> RepoSequence | None:
        """Return a sequence with the given accession if it exists in the isolate,
        else None.
        """
        return self._sequences_by_accession.get(accession)

    def dict(self, exclude_contents: bool = False):
        isolate_dict = {
            "id": self.id,
            "legacy_id": self.legacy_id,
            "name": {"type": self.name.type, "value": self.name.value},
        }

        if not exclude_contents:
            isolate_dict["sequences"] = [sequence.dict() for sequence in self.sequences]

        return isolate_dict

    def __eq__(self, other: "RepoIsolate") -> bool:
        if type(other) is not RepoIsolate:
            return False

        if self.id != other.id:
            return False

        if self.name.type != other.name.type:
            return False

        if self.name.value != other.name.value:
            return False

        if self.accessions != other.accessions:
            return False

        for accession in self.accessions:
            if self.get_sequence_by_accession(
                accession,
            ) != other.get_sequence_by_accession(accession):
                return False

        return True


class RepoOTU:
    """Represents an OTU in a Virtool reference repository."""

    def __init__(
        self,
        uuid: UUID,
        taxid: int,
        name: str,
        acronym: str = "",
        legacy_id: str | None = None,
        schema: OTUSchema | None = None,
        excluded_accessions: list[str] | None = None,
        isolates: list[RepoIsolate] | None = None,
        repr_isolate: UUID | None = None,
    ):
        self.id = uuid
        """The OTU id."""

        self.taxid = taxid
        """The NCBI Taxonomy id for this OTU."""

        self.name = name
        """The name of the OTU (eg. TMV for Tobacco mosaic virus)"""

        self.acronym = acronym
        """The OTU acronym (eg. TMV for Tobacco mosaic virus)."""

        self.legacy_id = legacy_id
        """A string based ID carried over from a legacy Virtool reference repository."""

        self.schema = schema
        """The schema of the OTU"""

        self.excluded_accessions = (
            set() if excluded_accessions is None else set(excluded_accessions)
        )
        """A set of accessions that should not be retrieved
        in future fetch operations"""

        self._isolates_by_id = (
            {} if isolates is None else {isolate.id: isolate for isolate in isolates}
        )
        """A dictionary of isolates indexed by isolate UUID"""

        self.repr_isolate = repr_isolate
        """The UUID of the representative isolate of this OTU"""

        self.__hash__ = None

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> "RepoOTU":
        """Build a new OTU from .dict() output."""
        return RepoOTU(
            uuid=data["id"],
            taxid=data["taxid"],
            name=data["name"],
            acronym=data.get("acronym"),
            schema=data.get("schema"),
            isolates=data.get("isolates"),
            repr_isolate=data.get("repr_isolate"),
        )

    @property
    def isolates(self) -> list[RepoIsolate]:
        """Isolates contained in this OTU."""
        return list(self._isolates_by_id.values())

    @property
    def isolate_ids(self) -> set[UUID]:
        """A set of UUIDs for isolates in the OTU."""
        return set(self._isolates_by_id.keys())

    @property
    def accessions(self) -> set[str]:
        """A set of accessions contained in this isolate."""
        accessions = set()
        for isolate in self.isolates:
            accessions.update(isolate.accessions)

        return accessions

    @property
    def sequence_ids(self) -> set[UUID]:
        """A set of UUIDs for sequences in the OTU."""
        sequence_ids = set()
        for isolate in self.isolates:
            sequence_ids.update(isolate.sequence_ids)

        return sequence_ids

    @property
    def blocked_accessions(self) -> set[str]:
        """Accessions that should not be considered for addition to the OTU.

        This includes accessions that already exist in the OTU and accessions that have
        been explicitly excluded.
        """
        return self.accessions | self.excluded_accessions

    def __repr__(self) -> str:
        """Return a shorthand representation of the OTU's contents."""
        return (
            f"EventSourcedRepoOTU(id='{self.id}', taxid='{self.taxid}', "
            f"name={self.name}, accessions={list(self.accessions)})"
        )

    def add_isolate(self, isolate: RepoIsolate) -> None:
        """Add an isolate to the OTU."""
        self._isolates_by_id[isolate.id] = isolate

    def add_sequence(self, sequence: RepoSequence, isolate_id: UUID) -> None:
        """Add a sequence to a given isolate."""
        self._isolates_by_id[isolate_id].add_sequence(sequence)

    def get_isolate(self, isolate_id: UUID) -> RepoIsolate | None:
        """Get isolate associated with a given ID.

        Returns None if no such isolate exists.

        :param isolate_id: The UUID of the isolate to retrieve
        :return: the isolate or ``None``
        """
        return self._isolates_by_id.get(isolate_id)

    def get_sequence_by_accession(
        self,
        accession: str,
    ) -> RepoSequence | None:
        """Return a sequence corresponding to given accession
        if it exists in this OTU.
        """
        if accession not in self.accessions:
            return None

        for isolate in self.isolates:
            if (sequence := isolate.get_sequence_by_accession(accession)) is not None:
                return sequence

        raise ValueError(f"Accession {accession} found in index, but not in data")

    def get_isolate_id_by_name(self, name: IsolateName) -> UUID | None:
        """Get the ID for the isolate with the passed ``name``.

        Returns None if no such isolate exists.

        :param name: The name of the isolate to retrieve
        :return: The isolate ID or ``None``

        """
        for isolate in self.isolates:
            if isolate.name == name:
                return isolate.id

        return None

    def dict(self, exclude_contents: bool = False):
        otu_dict = {
            "id": self.id,
            "acronym": self.acronym,
            "excluded_accessions": list(self.excluded_accessions),
            "legacy_id": self.legacy_id,
            "name": self.name,
            "schema": self.schema.model_dump() if self.schema is not None else None,
            "taxid": self.taxid,
        }

        if not exclude_contents:
            otu_dict["isolates"] = [
                isolate.dict(exclude_contents=False) for isolate in self.isolates
            ]

        return otu_dict

    def __eq__(self, other: "RepoOTU") -> bool:
        """Check if this OTU is equal to ``other``."""
        for isolate_id in self.isolate_ids:
            if self.get_isolate(isolate_id) != other.get_isolate(isolate_id):
                return False

        return (
            isinstance(other, RepoOTU)
            and self.accessions == other.accessions
            and self.acronym == other.acronym
            and self.id == other.id
            and self.isolate_ids == other.isolate_ids
            and self.legacy_id == other.legacy_id
            and self.name == other.name
            and self.taxid == other.taxid
            and self.repr_isolate == other.repr_isolate
            and self.schema == other.schema
        )

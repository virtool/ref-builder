import datetime
from uuid import UUID

from pydantic import UUID4, BaseModel, Field, field_serializer, field_validator

from ref_builder.models import Molecule
from ref_builder.plan import Plan
from ref_builder.utils import Accession, DataType, IsolateName


class RepoMeta(BaseModel):
    """Represents the metadata for a Virtool reference repository."""

    id: UUID4
    """The repository id."""

    created_at: datetime.datetime
    """The date and time the repository was created."""

    data_type: DataType
    """The repository data type."""

    name: str
    """The repository name."""

    organism: str
    """The repository organism."""


class RepoSettings(BaseModel):
    """The default settings of a Virtool reference repository."""

    default_segment_length_tolerance: float = Field(0.03, ge=0.0, le=1.0)
    """The deviation a sequence is allowed from its plan segment's length before it
    fails validation.
    """


class RepoSequence(BaseModel):
    """Represents a sequence in a Virtool reference repository."""

    id: UUID4
    """The sequence id."""

    accession: Accession
    """The sequence accession."""

    definition: str
    """The sequence definition."""

    legacy_id: str | None
    """A string based ID carried over from a legacy Virtool reference repository.

    It the sequence was not migrated from a legacy repository, this will be `None`.
    """

    sequence: str
    """The sequence."""

    segment: UUID4
    """The sequence segment."""

    @field_validator("accession", mode="before")
    @classmethod
    def convert_accession(cls: "RepoSequence", value: Accession | str) -> Accession:
        """Convert the accession to an Accession object."""
        if isinstance(value, Accession):
            return value

        if isinstance(value, str):
            return Accession.from_string(value)

        raise ValueError(f"Invalid type for accession: {type(value)}")

    @field_serializer("accession")
    @classmethod
    def serialize_accession(cls: "RepoSequence", accession: Accession) -> str:
        """Serialize the accession to a string."""
        return str(accession)


class RepoIsolate(BaseModel):
    """Represents an isolate in a Virtool reference repository."""

    id: UUID4
    """The isolate id."""

    legacy_id: str | None
    """A string based ID carried over from a legacy Virtool reference repository.

    It the isolate was not migrated from a legacy repository, this will be `None`.
    """

    name: IsolateName | None
    """The isolate's source name metadata."""

    sequences: list[RepoSequence]
    """A list of sequences contained by this isolate."""

    @property
    def accessions(self) -> set[str]:
        """A set of accession numbers for sequences in the isolate."""
        return {sequence.accession.key for sequence in self.sequences}

    @property
    def sequence_ids(self) -> set[UUID]:
        """A set of UUIDs for sequences in the isolate."""
        return {sequence.id for sequence in self.sequences}

    def add_sequence(self, sequence: RepoSequence) -> None:
        """Add a sequence to the isolate."""
        self.sequences.append(sequence)

    def delete_sequence(self, sequence_id: UUID) -> None:
        """Delete a sequence from a given isolate."""
        sequence = self.get_sequence_by_id(sequence_id)

        if not sequence:
            raise ValueError("This sequence ID does not exist")

        for sequence in self.sequences:
            if sequence.id == sequence_id:
                self.sequences.remove(sequence)

    def get_sequence_by_accession(
        self,
        accession: str,
    ) -> RepoSequence | None:
        """Return a sequence with the given accession if it exists in the isolate,
        else None.
        """
        for sequence in self.sequences:
            if sequence.accession.key == accession:
                return sequence

        return None

    def get_sequence_by_id(self, sequence_id: UUID) -> RepoSequence | None:
        """Return a sequence with the given ID if it exists in the isolate,
        else None.
        """
        for sequence in self.sequences:
            if sequence.id == sequence_id:
                return sequence

        return None

    @field_validator("name", mode="before")
    @classmethod
    def convert_name(
        cls: "RepoIsolate",
        value: dict | IsolateName | None,
    ) -> IsolateName | None:
        """Convert the name to an IsolateName object."""
        if value is None:
            return value

        if isinstance(value, IsolateName):
            return value

        if isinstance(value, dict):
            return IsolateName(**value)

        raise ValueError(f"Invalid type for name: {type(value)}")

    @field_serializer("name")
    def serialize_name(self, name: IsolateName | None) -> dict[str, str] | None:
        """Serialize the isolate name."""
        if name is None:
            return None

        return {
            "type": name.type,
            "value": name.value,
        }


class RepoOTU(BaseModel):
    """Represents an OTU in a Virtool reference repository."""

    id: UUID4
    """The OTU id."""

    acronym: str
    """The OTU acronym (eg. TMV for Tobacco mosaic virus)."""

    excluded_accessions: set[str]
    """A set of accessions that should not be retrieved in future fetch operations."""

    legacy_id: str | None
    """A string based ID carried over from a legacy Virtool reference repository."""

    name: str
    """The name of the OTU (eg. TMV for Tobacco mosaic virus)"""

    representative_isolate: UUID4 | None
    """The UUID of the representative isolate for the OTU."""

    molecule: Molecule
    """The type of molecular information contained in this OTU."""

    plan: Plan
    """The plan."""

    taxid: int
    """The NCBI Taxonomy id for this OTU."""

    isolates: list[RepoIsolate]
    """Isolates contained in this OTU."""

    _isolates_by_id: dict[UUID4:RepoIsolate]
    """A dictionary of isolates indexed by isolate UUID"""

    _sequences_by_id: dict[UUID4:RepoSequence]
    """A dictionary of sequences indexed by sequence UUID"""

    def __init__(self, **data) -> None:
        super().__init__(**data)

        self._sequences_by_id = {}
        self._isolates_by_id = {}
        for isolate in self.isolates:
            self._isolates_by_id[isolate.id] = isolate
            for sequence in isolate.sequences:
                self._sequences_by_id[sequence.id] = sequence

    @property
    def accessions(self) -> set[str]:
        """A set of accessions contained in this isolate."""
        return {sequence.accession.key for sequence in self._sequences_by_id.values()}

    @property
    def blocked_accessions(self) -> set[str]:
        """Accessions that should not be considered for addition to the OTU.

        This includes accessions that already exist in the OTU and accessions that have
        been explicitly excluded.
        """
        return self.accessions | self.excluded_accessions

    @property
    def isolate_ids(self) -> set[UUID4]:
        """A set of UUIDs for isolates in the OTU."""
        return set(self._isolates_by_id.keys())

    def add_isolate(self, isolate: RepoIsolate) -> None:
        """Add an isolate to the OTU."""
        self.isolates.append(isolate)
        self._isolates_by_id[isolate.id] = isolate

        for sequence in isolate.sequences:
            self.add_sequence(sequence)

    def add_sequence(self, sequence: RepoSequence) -> None:
        """Add a sequence to a given isolate."""
        self._sequences_by_id[sequence.id] = sequence

    def get_sequence_by_id(self, sequence_id: UUID) -> RepoSequence | None:
        return self._sequences_by_id.get(sequence_id)

    def delete_isolate(self, isolate_id: UUID4) -> None:
        """Remove an isolate from the OTU."""
        if self._isolates_by_id.get(isolate_id) is None:
            raise ValueError(f"Isolate {isolate_id} does not exist")

        for isolate in self.isolates:
            if isolate.id == isolate_id:
                for sequence_id in isolate.sequence_ids:
                    if self.get_sequence_by_id(sequence_id) is None:
                        raise ValueError(
                            f"Sequence {sequence_id} not found in the sequence list"
                        )
                    self.delete_sequence(sequence_id)

                self.isolates.remove(isolate)

                self._isolates_by_id.pop(isolate_id)

                break

    def delete_sequence(self, sequence_id: UUID4) -> None:
        """Delete a sequence from a given isolate. Used only during rehydration."""
        self._sequences_by_id.pop(sequence_id)

    def get_isolate(self, isolate_id: UUID4) -> RepoIsolate | None:
        """Get isolate associated with a given ID.

        Returns None if no such isolate exists.

        :param isolate_id: The UUID of the isolate to retrieve
        :return: the isolate or ``None``
        """
        return self._isolates_by_id.get(isolate_id)

    def get_isolate_id_by_name(self, name: IsolateName) -> UUID4 | None:
        """Get the ID for the isolate with the passed ``name``.

        Returns None if no such isolate exists.

        :param name: The name of the isolate to retrieve
        :return: The isolate ID or ``None``

        """
        for isolate in self.isolates:
            if isolate.name == name:
                return isolate.id

        return None

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

    def get_isolate_ids_containing_sequence_id(self, sequence_id: UUID4) -> set[UUID4]:
        """Return a set of isolate IDs where the isolate contains the given sequence."""
        containing_isolate_ids = set()

        if sequence_id not in self._sequences_by_id:
            return containing_isolate_ids

        for isolate in self.isolates:
            if sequence_id in isolate.sequence_ids:
                containing_isolate_ids.add(isolate.id)

        if containing_isolate_ids:
            return containing_isolate_ids

        raise ValueError(f"Sequence ID {sequence_id} found in index, but not in data")

    def link_sequence(
        self, isolate_id: UUID4, sequence_id: UUID4
    ) -> RepoSequence | None:
        """Link the given sequence to the given isolate."""
        self.get_isolate(isolate_id).add_sequence(self.get_sequence_by_id(sequence_id))

        return self.get_isolate(isolate_id).get_sequence_by_id(sequence_id)

    def unlink_sequence(self, isolate_id: UUID4, sequence_id: UUID4) -> None:
        """Unlink the given sequence from the given isolate. Used only during rehydration."""
        self.get_isolate(isolate_id).delete_sequence(sequence_id)

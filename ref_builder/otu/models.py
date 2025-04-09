from collections import Counter
from uuid import UUID

from pydantic import (
    UUID4,
    BaseModel,
    ConfigDict,
    Field,
    field_serializer,
    field_validator,
    model_validator,
)
from pydantic_core import PydanticCustomError

from ref_builder.models import Molecule
from ref_builder.plan import Plan
from ref_builder.resources import RepoIsolate, RepoSequence
from ref_builder.utils import Accession, IsolateName, is_refseq


class SequenceBase(BaseModel):
    """A class representing a sequence with basic validation."""

    model_config = ConfigDict(validate_assignment=True)

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

    @property
    def refseq(self) -> bool:
        """Return True if this sequence was sourced from NCBI's RefSeq database."""
        return is_refseq(self.accession.key)

    @field_validator("accession", mode="before")
    @classmethod
    def convert_accession(cls, value: Accession | str) -> Accession:
        """Convert the accession to an Accession object."""
        if isinstance(value, Accession):
            return value

        if isinstance(value, str):
            return Accession.from_string(value)

        raise ValueError(f"Invalid type for accession: {type(value)}")

    @field_serializer("accession")
    @classmethod
    def serialize_accession(cls, accession: Accession) -> str:
        """Serialize the accession to a string."""
        return str(accession)


class Sequence(SequenceBase):
    """A class representing a sequence with full validation."""


class IsolateBase(BaseModel):
    """A class representing an isolate with basic validation."""

    id: UUID4
    """The isolate id."""

    legacy_id: str | None
    """A string based ID carried over from a legacy Virtool reference repository.

    It the isolate was not migrated from a legacy repository, this will be `None`.
    """

    name: IsolateName | None
    """The isolate's name."""

    sequences: list[SequenceBase]
    """The isolates sequences."""

    @property
    def refseq(self) -> bool:
        """Return True if this isolate was sourced from NCBI's RefSeq database."""
        if self.sequences:
            return is_refseq(self.sequences[0].accession.key)

        return False

    def get_sequence_by_accession(
        self,
        accession: Accession,
    ) -> SequenceBase | None:
        """Get a sequence by its accession.

        Return ``None`` if the sequence is not found.

        :param accession: the accession of the sequence to retrieve
        :return: the sequence with the given accession, or None
        """
        for sequence in self.sequences:
            if sequence.accession == accession:
                return sequence

        return None

    def get_sequence_by_id(self, sequence_id: UUID) -> SequenceBase | None:
        """Get a sequence by its ID.

        Return ``None`` if the sequence is not found.

        :param sequence_id: the ID of the sequence to retrieve
        :return: the sequence with the given ID, or None
        """
        for sequence in self.sequences:
            if sequence.id == sequence_id:
                return sequence

        return None

    @field_validator("sequences", mode="before")
    @classmethod
    def convert_sequence_models(
        cls,
        v: list[dict | SequenceBase | Sequence | RepoSequence | BaseModel],
    ) -> list[SequenceBase]:
        """Automatically revalidate sequence if not already validated."""
        if not v or isinstance(v[0], dict | SequenceBase):
            return v

        return [SequenceBase.model_validate(sequence.model_dump()) for sequence in v]


class Isolate(IsolateBase):
    """A class representing an isolate with full validation."""

    model_config = ConfigDict(validate_assignment=True)

    sequences: list[Sequence] = Field(min_length=1)
    """The isolates sequences.

    A valid isolate must have at least one sequence.
    """

    @field_validator("sequences", mode="before")
    @classmethod
    def convert_sequence_models(
        cls,
        v: list[dict | SequenceBase | Sequence | RepoSequence | BaseModel],
    ) -> list[Sequence]:
        if not v or isinstance(v[0], dict | SequenceBase):
            return v

        return [Sequence.model_validate(sequence.model_dump()) for sequence in v]


class OTUBase(BaseModel):
    """A class representing an OTU with basic validation."""

    id: UUID4
    """The OTU id."""

    acronym: str
    """The OTU acronym (eg. TMV for Tobacco mosaic virus)."""

    excluded_accessions: set[str]
    """A set of accessions that should not be retrieved in future fetch operations."""

    isolates: list[IsolateBase]
    """Isolates contained in this OTU."""

    legacy_id: str | None
    """A string based ID carried over from a legacy Virtool reference repository."""

    molecule: Molecule
    """The type of molecular information contained in this OTU."""

    name: str
    """The name of the OTU (eg. TMV for Tobacco mosaic virus)"""

    plan: Plan
    """The plan for the OTU."""

    representative_isolate: UUID4 | None
    """The UUID of the representative isolate of this OTU"""

    taxid: int
    """The NCBI Taxonomy id for this OTU."""

    @property
    def sequences(self) -> list[SequenceBase]:
        """Sequences contained in this OTU."""
        return [sequence for isolate in self.isolates for sequence in isolate.sequences]

    @field_validator("isolates", mode="before")
    @classmethod
    def convert_isolate_models(
        cls,
        v: list[dict | IsolateBase | Isolate | RepoIsolate | BaseModel],
    ) -> list[IsolateBase]:
        """Automatically revalidate isolates if not already validated."""
        if not v or isinstance(v[0], dict | IsolateBase):
            return v

        return [IsolateBase.model_validate(isolate.model_dump()) for isolate in v]


class OTU(OTUBase):
    """A class representing an OTU with full validation."""

    model_config = ConfigDict(validate_assignment=True)

    isolates: list[Isolate] = Field(min_length=1)
    """Isolates contained in this OTU.

    A valid OTU must have at least one isolate.
    """

    representative_isolate: UUID4
    """The UUID of the representative isolate of this OTU.

    A valid OTU must have a representative isolate.
    """

    @property
    def sequences(self) -> list[Sequence]:
        """Sequences contained in this OTU."""
        return [sequence for isolate in self.isolates for sequence in isolate.sequences]

    @field_validator("isolates", mode="before")
    @classmethod
    def convert_isolate_models(
        cls,
        v: list[dict | Isolate | IsolateBase | RepoIsolate | BaseModel],
    ) -> list[Isolate]:
        """Automatically revalidate isolates if not already validated."""
        if not v or isinstance(v[0], dict | Isolate):
            return v

        return [Isolate.model_validate(isolate.model_dump()) for isolate in v]

    @model_validator(mode="after")
    def check_excluded_accessions(self) -> "OTU":
        """Ensure that excluded accessions are not in the OTU."""
        if accessions := self.excluded_accessions & {
            sequence.accession.key for sequence in self.sequences
        }:
            raise ValueError(
                f"Excluded accessions found in the OTU: {', '.join(accessions)}"
            )

        return self

    @model_validator(mode="after")
    def check_representative_isolate(self) -> "OTU":
        """Ensure that the default isolate is in the OTU."""
        if self.representative_isolate not in {isolate.id for isolate in self.isolates}:
            raise ValueError("Representative isolate must be in the OTU")

        return self

    @model_validator(mode="after")
    def check_unique_isolate_names(self) -> "OTU":
        """Ensure there are no duplicate isolate names in the OTU."""
        counts = Counter(isolate.name for isolate in self.isolates)

        duplicates = ", ".join(
            [str(name) for name, count in counts.items() if name and count > 1]
        )

        if duplicates:
            raise ValueError(
                f"Isolate names must be unique. Non-unique names: {duplicates}",
            )

        return self

    @model_validator(mode="after")
    def check_isolates_against_plan(self) -> "OTU":
        """Check that all isolates satisfy the OTU's plan."""
        for isolate in self.isolates:

            for sequence in isolate.sequences:
                segment = self.plan.get_segment_by_id(sequence.segment)
                if segment is None:
                    raise PydanticCustomError(
                        "segment_not_found",
                        "Sequence segment {sequence_segment} was not found in "
                        + "the list of segments: {plan_segments}.",
                        {
                            "isolate_id": isolate.id,
                            "sequence_segment": sequence.segment,
                            "plan_segments": list(self.plan.segment_ids)
                        }
                    )

                min_length = int(segment.length * (1.0 - segment.length_tolerance))
                max_length = int(segment.length * (1.0 + segment.length_tolerance))

                if len(sequence.sequence) < min_length:
                    raise PydanticCustomError(
                        "sequence_too_short",
                        "Sequence based on {sequence_accession} does not pass validation "
                        + "against segment {segment_id} "
                        + "({sequence_length} < {min_sequence_length})",
                        {
                            "isolate_id": isolate.id,
                            "sequence_id": sequence.id,
                            "sequence_accession": sequence.accession,
                            "sequence_length": len(sequence.sequence),
                            "segment_id": segment.id,
                            "segment_reference_length": segment.length,
                            "min_sequence_length": min_length,
                        }
                    )

                if len(sequence.sequence) > max_length:
                    raise PydanticCustomError(
                        "sequence_too_long",
                        "Sequence based on {sequence_accession} does not pass validation "
                        + "against segment {segment_id}"
                        + "({sequence_length} > {max_sequence_length})",
                        {
                            "isolate_id": isolate.id,
                            "sequence_id": sequence.id,
                            "sequence_accession": sequence.accession,
                            "sequence_length": len(sequence.sequence),
                            "segment_id": segment.id,
                            "segment_reference_length": segment.length,
                            "max_sequence_length": max_length,
                        }
                    )

        return self

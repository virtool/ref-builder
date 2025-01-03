from collections import Counter
from uuid import UUID

from pydantic import UUID4, BaseModel, Field, model_validator

from ref_builder.models import Molecule
from ref_builder.plan import Plan
from ref_builder.resources import RepoIsolate, RepoSequence
from ref_builder.utils import Accession, IsolateName


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

    sequences: list[RepoSequence]
    """The isolates sequences."""

    def get_sequence_by_accession(
        self,
        accession: Accession,
    ) -> RepoSequence | None:
        """Get a sequence by its accession.

        Return ``None`` if the sequence is not found.

        :param accession: the accession of the sequence to retrieve
        :return: the sequence with the given accession, or None
        """
        for sequence in self.sequences:
            if sequence.accession == accession:
                return sequence

        return None

    def get_sequence_by_id(self, sequence_id: UUID) -> RepoSequence | None:
        """Get a sequence by its ID.

        Return ``None`` if the sequence is not found.

        :param sequence_id: the ID of the sequence to retrieve
        :return: the sequence with the given ID, or None
        """
        for sequence in self.sequences:
            if sequence.id == sequence_id:
                return sequence

        return None


class Isolate(IsolateBase):
    """A class representing an isolate with full validation."""

    sequences: list[RepoSequence] = Field(min_length=1)
    """The isolates sequences.

    A valid isolate must have at least one sequence.
    """


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
    """The schema of the OTU"""

    representative_isolate: UUID4 | None
    """The UUID of the representative isolate of this OTU"""

    sequences: list[RepoSequence]
    """Sequences contained in this OTU."""

    taxid: int
    """The NCBI Taxonomy id for this OTU."""


class OTU(OTUBase):
    """A class representing an OTU with full validation."""

    isolates: list[RepoIsolate] = Field(min_length=1)
    """Isolates contained in this OTU.

    A valid OTU must have at least one isolate.
    """

    representative_isolate: UUID4
    """The UUID of the representative isolate of this OTU.

    A valid OTU must have a representative isolate.
    """

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
        """Check that all isolates satisfy the OTU's plan.

        TODO: Implement this method.
        """
        return self

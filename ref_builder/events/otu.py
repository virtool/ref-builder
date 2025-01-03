from pydantic import UUID4, field_serializer

from ref_builder.events.base import ApplicableEvent, EventData, Event, OTUQuery
from ref_builder.resources import RepoOTU
from ref_builder.models import Molecule
from ref_builder.plan import Plan


class CreateOTUData(EventData):
    """The data for a :class:`CreateOTU` event."""

    id: UUID4
    acronym: str
    legacy_id: str | None
    molecule: Molecule
    name: str
    taxid: int
    plan: Plan


class CreateOTU(Event):
    """An event that creates a new OTU."""

    data: CreateOTUData
    query: OTUQuery

    def apply(self) -> RepoOTU:
        return RepoOTU(
            id=self.query.otu_id,
            acronym=self.data.acronym,
            excluded_accessions=set(),
            isolates=[],
            legacy_id=self.data.legacy_id,
            molecule=self.data.molecule,
            name=self.data.name,
            repr_isolate=None,
            plan=self.data.plan,
            taxid=self.data.taxid,
        )


class CreatePlanData(EventData):
    """The data for a :class:`CreatePlan` event."""

    plan: Plan


class CreatePlan(ApplicableEvent):
    """An event that sets the isolate plan for an OTU."""

    data: CreatePlanData
    query: OTUQuery

    def apply(self, otu: RepoOTU) -> RepoOTU:
        """Apply changed plan to OTU and return."""
        otu.plan = self.data.plan

        return otu


class SetReprIsolateData(EventData):
    """The data for a :class:`SetReprIsolate` event."""

    isolate_id: UUID4


class SetReprIsolate(ApplicableEvent):
    """An event that sets the representative isolate for an OTU."""

    data: SetReprIsolateData
    query: OTUQuery

    def apply(self, otu: RepoOTU) -> RepoOTU:
        """Update the OTU's representative isolate and return."""
        otu.repr_isolate = self.data.isolate_id

        return otu


class UpdateAllowedAccessionsData(EventData):
    """The data for an UpdateAllowedAccessions event."""

    accessions: set[str]
    allow: bool

    @field_serializer("accessions")
    def serialize_accessions(self, accessions: set[str]) -> list[str]:
        return sorted(accessions)


class UpdateAllowedAccessions(ApplicableEvent):
    """An event that changes the OTU excluded accessions collection.

    This event is emitted when Genbank accessions are either
    allowed or disallowed from inclusion in the reference.
    """

    data: UpdateAllowedAccessionsData
    query: OTUQuery

    def apply(self, otu: RepoOTU) -> RepoOTU:
        """Add accession allowance changes to OTU and return."""

        if self.data.allow:
            for accession in self.data.accessions:
                otu.excluded_accessions.discard(accession)

        else:
            for accession in self.data.accessions:
                otu.excluded_accessions.add(accession)

        return otu

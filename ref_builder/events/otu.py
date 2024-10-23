from pydantic import UUID4

from ref_builder.events.base import ApplicableEvent, EventData, Event, OTUQuery
from ref_builder.resources import RepoOTU
from ref_builder.models import Molecule
from ref_builder.plan import MonopartitePlan, MultipartitePlan


class CreateOTUData(EventData):
    """The data for a :class:`CreateOTU` event."""

    id: UUID4
    acronym: str
    legacy_id: str | None
    molecule: Molecule
    name: str
    taxid: int
    plan: MonopartitePlan | MultipartitePlan


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


class ExcludeAccessionData(EventData):
    """The data for a :class:`ExcludeAccession` event."""

    accession: str


class ExcludeAccession(ApplicableEvent):
    """An accession exclusion event.

    This event is emitted when a Genbank accession is not going to be allowed in the
    reference.
    """

    data: ExcludeAccessionData
    query: OTUQuery

    def apply(self, otu: RepoOTU) -> RepoOTU:
        """Add excluded accession to OTU and return."""
        otu.excluded_accessions.add(self.data.accession)

        return otu


class CreatePlanData(EventData):
    """The data for a :class:`CreatePlan` event."""

    plan: MonopartitePlan | MultipartitePlan


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

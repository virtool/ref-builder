from pydantic import UUID4

from ref_builder.events.base import EventData, Event, OTUQuery
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


class ExcludeAccessionData(EventData):
    """The data for a :class:`ExcludeAccession` event."""

    accession: str


class ExcludeAccession(Event):
    """An accession exclusion event.

    This event is emitted when a Genbank accession is not going to be allowed in the
    reference.
    """

    data: ExcludeAccessionData
    query: OTUQuery


class CreatePlanData(EventData):
    """The data for a :class:`CreatePlan` event."""

    plan: MonopartitePlan | MultipartitePlan


class CreatePlan(Event):
    """An event that sets the isolate plan for an OTU."""

    data: CreatePlanData
    query: OTUQuery


class SetReprIsolateData(EventData):
    """The data for a :class:`SetReprIsolate` event."""

    isolate_id: UUID4


class SetReprIsolate(Event):
    """An event that sets the representative isolate for an OTU."""

    data: SetReprIsolateData
    query: OTUQuery

"""Events that are emitted to the repository event log."""

import datetime

from pydantic import UUID4, AliasChoices, BaseModel, Field, computed_field

from ref_builder.models import Molecule
from ref_builder.plan import MonopartitePlan, MultipartitePlan
from ref_builder.utils import Accession, IsolateName


class EventQuery(BaseModel):
    """A base class for a query that targets an event at a specific resource."""


class RepoQuery(EventQuery):
    """An event query that targets an event at the repository."""

    repository_id: UUID4


class OTUQuery(EventQuery):
    """An event query that targets an event at an OTU."""

    otu_id: UUID4


class IsolateQuery(OTUQuery):
    """An event query that targets an event at an isolate in a specific OTU."""

    isolate_id: UUID4


class SequenceQuery(IsolateQuery):
    """An event query that targets an event at a sequence in a specific isolate and
    OTU.
    """

    sequence_id: UUID4


class EventData(BaseModel):
    """Represents the data for an event."""


class Event(BaseModel):
    """The base event."""

    id: int
    """The unique identifier for the event.

    Event IDs are serially incremented integers.
    """

    data: EventData
    """The data associated with the event."""

    query: EventQuery
    """The query targeting the event at a specific resource."""

    timestamp: datetime.datetime
    """When the event occurred."""

    @computed_field
    @property
    def type(self) -> str:
        """The type of the event as a string."""
        return self.__class__.__name__


class CreateRepoData(EventData):
    """The data for a :class:`CreateRepo` event."""

    id: UUID4
    data_type: str
    name: str
    organism: str


class CreateRepo(Event):
    """An event that creates a new repository.

    This event is always the first event in a repository's event log.
    """

    data: CreateRepoData
    query: RepoQuery


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


class CreateIsolateData(EventData):
    """The data for a :class:`CreateIsolate` event."""

    id: UUID4
    legacy_id: str | None
    name: IsolateName | None


class CreateIsolate(Event):
    """An event that creates an isolate for a specific OTU."""

    data: CreateIsolateData
    query: IsolateQuery


class DeleteIsolateData(EventData):
    """The data for a :class:`DeleteIsolate` event."""

    rationale: str


class DeleteIsolate(Event):
    """An isolate deletion event."""

    data: DeleteIsolateData
    query: IsolateQuery


class CreateSequenceData(EventData):
    """The data for a :class:`CreateSequence` event."""

    id: UUID4
    accession: Accession
    definition: str
    legacy_id: str | None
    segment: str
    sequence: str


class CreateSequence(Event):
    """An event that creates a sequence for a specific isolate and OTU."""

    data: CreateSequenceData
    query: SequenceQuery


class DeleteSequenceData(EventData):
    """The data for a :class:`DeleteSequence` event."""

    replacement: UUID4
    rationale: str


class DeleteSequence(Event):
    """An event that deletes a sequence.

    The second part of a sequence replacement.
    """

    data: DeleteSequenceData
    query: SequenceQuery


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

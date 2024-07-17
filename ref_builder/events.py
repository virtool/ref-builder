import datetime

from pydantic import UUID4, BaseModel, Field, computed_field

from ref_builder.models import Molecule
from ref_builder.schema import OTUSchema, Segment
from ref_builder.utils import IsolateName


class EventQuery(BaseModel):
    """A base class for representing the query targeting an event at a specific
    resource.
    """


class RepoQuery(EventQuery):
    """An event query that targets an event at a repository."""

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
    """Represents the data for an event.

    Different event data classes are used
    """


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
        return self.__class__.__name__


class CreateRepoData(EventData):
    """The data for a repository creation event (`CreateRepo`)."""

    id: UUID4
    data_type: str
    name: str
    organism: str


class CreateRepo(Event):
    """A repo creation event."""

    data: CreateRepoData
    query: RepoQuery


class CreateOTUData(EventData):
    """The data for the creation of a new OTU ('CreateOTU')."""

    id: UUID4
    acronym: str
    legacy_id: str | None
    name: str
    taxid: int
    rep_isolate: UUID4 | None
    otu_schema: OTUSchema | None = Field(None, alias="schema")


class CreateOTU(Event):
    """An OTU creation event."""

    data: CreateOTUData
    query: OTUQuery


class CreateIsolateData(EventData):
    """The data for the creation of a new isolate."""

    id: UUID4
    legacy_id: str | None
    name: IsolateName


class CreateIsolate(Event):
    """An isolate creation event."""

    data: CreateIsolateData
    query: IsolateQuery


class CreateSequenceData(EventData):
    """The data for the creation of a new sequence."""

    id: UUID4
    accession: str
    definition: str
    legacy_id: str | None
    segment: str
    sequence: str


class CreateSequence(Event):
    """A sequence creation event."""

    data: CreateSequenceData
    query: SequenceQuery


class ExcludeAccessionData(EventData):
    """The data for the exclusion of an accession."""

    accession: str


class ExcludeAccession(Event):
    """An accession exclusion event.

    This event is emitted when a Genbank accession is not going to be allowed in the
    reference.
    """

    data: ExcludeAccessionData
    query: OTUQuery


class CreateSchemaData(EventData):
    """The data for the creation of OTU schema data."""

    molecule: Molecule
    segments: list[Segment]


class CreateSchema(Event):
    """An OTU schema creation event."""

    data: CreateSchemaData
    query: OTUQuery

import datetime

from pydantic import BaseModel, computed_field, UUID4


class EventQuery(BaseModel):
    """A base class for a query that targets an event at a specific resource."""


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

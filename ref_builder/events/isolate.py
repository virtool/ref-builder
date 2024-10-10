from pydantic import UUID4

from ref_builder.events.base import EventData, Event, IsolateQuery
from ref_builder.utils import IsolateName


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

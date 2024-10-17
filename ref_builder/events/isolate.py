from typing import Literal

from pydantic import UUID4

from ref_builder.events.base import EventData, Event, IsolateQuery, LinkSequenceQuery
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
    type: Literal["CreateIsolate"] = "CreateIsolate"


class LinkSequenceData(EventData):
    """The data for a :class:`LinkSequence` event."""


class LinkSequence(Event):
    """An event that links an existing sequence to an isolate."""

    data: LinkSequenceData
    query: LinkSequenceQuery
    type: Literal["LinkSequence"] = "LinkSequence"


class DeleteIsolateData(EventData):
    """The data for a :class:`DeleteIsolate` event."""

    rationale: str


class DeleteIsolate(Event):
    """An isolate deletion event."""

    data: DeleteIsolateData
    query: IsolateQuery
    type: Literal["DeleteIsolate"] = "DeleteIsolate"

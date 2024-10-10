from pydantic import UUID4

from ref_builder.events.base import EventData, Event, SequenceQuery
from ref_builder.utils import Accession


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

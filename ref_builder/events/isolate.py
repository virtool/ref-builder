from pydantic import UUID4

from ref_builder.resources import RepoOTU, RepoIsolate
from ref_builder.events.base import (
    ApplicableEvent,
    EventData,
    IsolateQuery,
)
from ref_builder.utils import IsolateName


class CreateIsolateData(EventData):
    """The data for a :class:`CreateIsolate` event."""

    id: UUID4
    legacy_id: str | None
    name: IsolateName | None


class CreateIsolate(ApplicableEvent):
    """An event that creates an isolate for a specific OTU."""

    data: CreateIsolateData
    query: IsolateQuery

    def apply(self, otu: RepoOTU) -> RepoOTU:
        """Add isolate to OTU and return."""
        otu.add_isolate(
            RepoIsolate(
                id=self.data.id,
                legacy_id=self.data.legacy_id,
                name=self.data.name,
                sequences=[],
            ),
        )

        return otu


class LinkSequenceData(EventData):
    """The data for a :class:`LinkSequence` event."""

    sequence_id: UUID4


class LinkSequence(ApplicableEvent):
    """An event that links an existing sequence to an isolate."""

    data: LinkSequenceData
    query: IsolateQuery

    def apply(self, otu: RepoOTU) -> RepoOTU:
        """Add specified sequence to specified isolate and return."""
        otu.link_sequence(
            isolate_id=self.query.isolate_id,
            sequence_id=self.data.sequence_id,
        )

        return otu


class UnlinkSequenceData(EventData):
    """The data for a :class:`LinkSequence` event."""

    sequence_id: UUID4


class UnlinkSequence(ApplicableEvent):
    """An event that unlinks an existing sequence from an isolate."""

    data: UnlinkSequenceData
    query: IsolateQuery

    def apply(self, otu: RepoOTU) -> RepoOTU:
        """Unlink specified sequence from specified isolate and return OTU."""
        otu.unlink_sequence(
            isolate_id=self.query.isolate_id,
            sequence_id=self.data.sequence_id,
        )

        return otu


class DeleteIsolateData(EventData):
    """The data for a :class:`DeleteIsolate` event."""

    rationale: str


class DeleteIsolate(ApplicableEvent):
    """An isolate deletion event."""

    data: DeleteIsolateData
    query: IsolateQuery

    def apply(self, otu: RepoOTU) -> RepoOTU:
        """Delete the specified isolate and return."""
        otu.delete_isolate(self.query.isolate_id)

        return otu

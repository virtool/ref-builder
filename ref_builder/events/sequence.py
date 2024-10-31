from pydantic import UUID4

from ref_builder.resources import RepoOTU, RepoSequence
from ref_builder.events.base import (
    ApplicableEvent,
    EventData,
    SequenceQuery,
)
from ref_builder.utils import Accession


class CreateSequenceData(EventData):
    """The data for a :class:`CreateSequence` event."""

    id: UUID4
    accession: Accession
    definition: str
    legacy_id: str | None
    segment: str
    sequence: str


class CreateSequence(ApplicableEvent):
    """An event that creates a sequence for a specific isolate and OTU."""

    data: CreateSequenceData
    query: SequenceQuery

    def apply(self, otu: RepoOTU) -> RepoOTU:
        """Create sequence in OTU and return."""
        otu.add_sequence(
            RepoSequence(
                id=self.data.id,
                accession=self.data.accession,
                definition=self.data.definition,
                legacy_id=self.data.legacy_id,
                segment=self.data.segment,
                sequence=self.data.sequence,
            ),
        )

        return otu


class DeleteSequenceData(EventData):
    """The data for a :class:`DeleteSequence` event."""

    sequence_id: UUID4
    replacement: UUID4
    rationale: str


class DeleteSequence(ApplicableEvent):
    """An event that deletes a sequence.

    The second part of a sequence replacement.
    """

    data: DeleteSequenceData
    query: SequenceQuery

    def apply(self, otu: RepoOTU) -> RepoOTU:
        """Delete sequence from OTU and return."""
        otu.delete_sequence(self.data.sequence_id)

        return otu

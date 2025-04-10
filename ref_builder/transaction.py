from uuid import UUID


class AbortTransactionError(Exception):
    """Raised when a modification lock is aborted."""


class Transaction:
    """A lock on the repository during modification.

    Locks prevent any other ``ref-builder`` processes from modifying the repository.

    They also make multi-event operations atomic by storing a marker for the last
    validated event. If the lock is released without committing the transaction, the
    marker is not updated and events beyond the marker are pruned.

    """

    def __init__(self) -> None:
        self._otu_ids = set()

    @property
    def affected_otu_ids(self) -> set[UUID]:
        """The IDs of OTUs affected by this transaction."""
        return self._otu_ids

    def add_otu_id(self, otu_id: UUID) -> None:
        """Associate a change in an OTU with this transaction.

        :param otu_id: The OTU ID to track.
        """
        self._otu_ids.add(otu_id)

    def abort(self) -> None:
        """Abort the transaction.

        The pending events will be discarded and no changes will be made to the
        repository.
        """
        raise AbortTransactionError

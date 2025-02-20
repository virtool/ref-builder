class AbortTransactionError(Exception):
    """Raised when a modification lock is aborted."""


class Transaction:
    """A lock on the repository during modification.

    Locks prevent any other ``ref-builder`` processes from modifying the repository.

    They also make multi-event operations atomic by storing a marker for the last
    validated event. If the lock is released without committing the transaction, the
    marker is not updated and events beyond the marker are pruned.

    """

    def abort(self) -> None:
        """Abort the transaction.

        The pending events will be discarded and no changes will be made to the
        repository.
        """
        raise AbortTransactionError

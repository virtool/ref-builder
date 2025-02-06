class LockConflictError(Exception):
    """Raised when a lock is attempted on a locked repository."""

    def __init__(self):
        super().__init__("Repository is already locked by another process.")


class LockRequiredError(Exception):
    """Raised when a transaction is attempted without a lock."""

    def __init__(self) -> None:
        super().__init__("Repository must be locked to use a transaction.")


class TransactionExistsError(Exception):
    """Raised when an error occurs during a transaction."""

    def __init__(self) -> None:
        super().__init__("An active transaction already exists.")


class TransactionRequiredError(Exception):
    """Raised when an error occurs during a transaction."""

    def __init__(self) -> None:
        super().__init__("Writing events requires a transaction and none was found.")

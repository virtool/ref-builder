import os
from pathlib import Path

from structlog import get_logger

from ref_builder.errors import LockConflictError

logger = get_logger()


class Lock:
    """A simple lock file implementation."""

    def __init__(self, path: Path) -> None:
        self.lock_path = path / "lock"
        self.locked = False

        if self.lock_path.exists():
            raise LockConflictError

    def lock(self) -> None:
        """Lock the repository by creating a lock file.

        Locking is idempotent. If the repository is already locked by this process, this
        method will do nothing.
        """
        if self.locked:
            logger.warning("A lock was attempted, but this process already has a lock.")
            return

        with open(self.lock_path, "w") as lock_file:
            lock_file.write(str(os.getpid()))

        self.locked = True

    def unlock(self) -> None:
        """Unlock the repository by removing the lock file."""
        self.lock_path.unlink()
        self.locked = False

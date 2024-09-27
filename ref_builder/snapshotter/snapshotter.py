from collections.abc import Generator, Iterable
from pathlib import Path
from uuid import UUID

from pydantic.dataclasses import dataclass
from structlog import get_logger

from ref_builder.resources import RepoOTU
from ref_builder.snapshotter.index import Index

logger = get_logger()


@dataclass
class Snapshot:
    """A snapshot of an OTU at a specific event."""

    at_event: int
    """The event at which the snapshot was taken."""

    otu: RepoOTU
    """The OTU that was snapshotted."""


class Snapshotter:
    """Load and cache OTU snapshots."""

    def __init__(self, path: Path) -> None:
        """Create a new snapshotter."""
        self.path = path
        """The path to the snapshot root directory."""

        self.path.mkdir(exist_ok=True, parents=True)

        self._index = Index(self.path / "index.db")
        """The index of all OTUs in the snapshotter."""

    def get_id_by_legacy_id(self, legacy_id: str) -> UUID | None:
        """Get an OTU ID by its legacy ID.

        Returns ``None`` if the legacy ID is not found.

        :param legacy_id: The legacy ID to search for.
        :return: The ID of the OTU with the given legacy ID or ``None``.

        """
        return self._index.get_id_by_legacy_id(legacy_id)

    def get_id_by_name(self, name: str) -> UUID | None:
        """Get an OTU ID by its name.

        Returns ``None`` if the name is not found.

        :param name: The name to search for.
        :return: The ID of the OTU with the given name or ``None``.

        """
        return self._index.get_id_by_name(name)

    def get_id_by_taxid(self, taxid: int) -> UUID | None:
        """Get an OTU ID by its taxonomy ID.

        Returns ``None`` if the taxonomy ID is not found.

        :param taxid: The taxonomy ID to search for.
        :return: The ID of the OTU with the given taxonomy ID or ``None``.
        """
        return self._index.get_id_by_taxid(taxid)

    @property
    def otu_ids(self) -> set[UUID]:
        """A list of OTU ids of snapshots."""
        return self._index.otu_ids

    def snapshot(self, otus: Iterable[RepoOTU], at_event: int) -> None:
        """Take a new snapshot."""
        for otu in otus:
            self.cache_otu(otu, at_event)

    def iter_otus(self) -> Generator[RepoOTU, None, None]:
        """Iterate over the OTUs in the snapshot."""
        for otu_id in self.otu_ids:
            yield self.load_by_id(otu_id).otu

    def cache_otu(
        self,
        otu: "RepoOTU",
        at_event: int,
    ) -> None:
        """Create a snapshot for a single OTU."""
        self._index.update(otu, at_event)

    def load_by_id(self, otu_id: UUID) -> Snapshot | None:
        """Load the most recently snapshotted form of an OTU by its ID.

        Returns ``None`` if the OTU is not found.

        :param otu_id: The ID of the OTU to load.
        :return: The OTU or ``None`` if it is not found.

        """
        return self._index.load(otu_id)

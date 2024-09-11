from collections.abc import Generator, Iterable
from dataclasses import dataclass
from pathlib import Path
from uuid import UUID

import orjson
from structlog import get_logger

from ref_builder.resources import RepoOTU
from ref_builder.snapshotter.otu import OTUSnapshot

logger = get_logger()


@dataclass
class Snapshot:
    """A snapshot of an OTU at a specific event."""

    at_event: int
    """The event at which the snapshot was taken."""

    otu: RepoOTU
    """The OTU that was snapshotted."""


@dataclass
class OTUKeys:
    """Stores indexable data about OTUs."""

    id: UUID
    taxid: int
    name: str
    acronym: str = ""
    legacy_id: str | None = None

    @classmethod
    def from_otu(cls, otu: RepoOTU) -> "OTUKeys":
        """Create a new OTUKeys instance from a ``RepoOTU``."""
        return OTUKeys(
            id=otu.id,
            taxid=otu.taxid,
            name=otu.name,
            acronym=otu.acronym,
            legacy_id=otu.legacy_id,
        )

    def dict(self) -> dict[str, str | int | UUID]:
        """Return a dictionary representation of the OTU keys."""
        return {
            "id": self.id,
            "name": self.name,
            "taxid": self.taxid,
            "acronym": self.acronym,
            "legacy_id": self.legacy_id,
        }


class Snapshotter:
    """Load and cache OTU snapshots."""

    def __init__(self, path: Path) -> None:
        """Create a new snapshotter."""
        self.path = path
        """The path to the snapshot root directory."""

        self.path.mkdir(exist_ok=True, parents=True)

        self._index_path = self.path / "index.json"
        """The path to the index data."""

        self._index = self._load_index()
        """The index data of this snapshot index."""

        if self._index is None:
            self._build_index()

    def get_id_by_legacy_id(self, legacy_id: str) -> UUID | None:
        """Get an OTU ID by its legacy ID.

        Returns ``None`` if the legacy ID is not found.

        :param legacy_id: The legacy ID to search for.
        :return: The ID of the OTU with the given legacy ID or ``None``.

        """
        legacy_id_index = {
            self._index[otu_id].legacy_id: otu_id for otu_id in self._index
        }

        legacy_id_index.pop(None, None)

        return legacy_id_index.get(
            legacy_id,
        )

    def get_id_by_name(self, name: str) -> UUID | None:
        """Get an OTU ID by its name.

        Returns ``None`` if the name is not found.

        :param name: The name to search for.
        :return: The ID of the OTU with the given name or ``None``.

        """
        return {self._index[otu_id].name: otu_id for otu_id in self._index}.get(name)

    def get_id_by_taxid(self, taxid: int) -> UUID | None:
        """Get an OTU ID by its taxonomy ID.

        Returns ``None`` if the taxonomy ID is not found.

        :param taxid: The taxonomy ID to search for.
        :return: The ID of the OTU with the given taxonomy ID or ``None``.
        """
        return {self._index[otu_id].taxid: otu_id for otu_id in self._index}.get(taxid)

    @property
    def otu_ids(self) -> set[UUID]:
        """A list of OTU ids of snapshots."""
        return set(self._index.keys())

    def snapshot(self, otus: Iterable[RepoOTU], at_event: int) -> None:
        """Take a new snapshot."""
        _index = {}

        for otu in otus:
            self.cache_otu(otu, at_event)

            _index[otu.id] = OTUKeys(
                id=otu.id,
                taxid=otu.taxid,
                name=otu.name,
                acronym=otu.acronym,
                legacy_id=otu.legacy_id,
            )

        self._index = _index

        with open(self._index_path, "wb") as f:
            f.write(
                orjson.dumps(
                    {str(otu_id): self._index[otu_id].dict() for otu_id in self._index},
                ),
            )

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
        logger.debug("Writing a snapshot", otu_id=otu.id)

        otu_snap = OTUSnapshot(self.path / f"{otu.id}")
        otu_snap.cache(otu, at_event)

        self._index[otu.id] = OTUKeys.from_otu(otu)

    def load_by_id(self, otu_id: UUID) -> Snapshot | None:
        """Load the most recently snapshotted form of an OTU by its ID.

        Returns ``None`` if the OTU is not found.

        :param otu_id: The ID of the OTU to load.
        :return: The OTU or ``None`` if it is not found.

        """
        try:
            otu_snapshot = OTUSnapshot(self.path / f"{otu_id}")
        except FileNotFoundError:
            return None

        otu = otu_snapshot.load()

        return Snapshot(at_event=otu_snapshot.at_event, otu=otu)

    def _build_index(self) -> None:
        """Build a new index from the contents of the snapshot cache directory."""
        index = {}

        for subpath in self.path.iterdir():
            try:
                otu_id = UUID(subpath.stem)
            except ValueError:
                continue

            otu = self.load_by_id(otu_id)
            if otu is None:
                raise FileNotFoundError("OTU not found")
            index[otu.id] = OTUKeys(
                id=otu.id,
                taxid=otu.taxid,
                name=otu.name,
                acronym=otu.acronym,
                legacy_id=otu.legacy_id,
            )

        logger.debug("Snapshot index built", index=index)

        self._index = index

    def _load_index(self) -> dict | None:
        """Load the index from disk."""
        try:
            with open(self._index_path, "rb") as f:
                index_dict = orjson.loads(f.read())
        except FileNotFoundError:
            return None

        _index = {}

        for key in index_dict:
            try:
                otu_id = UUID(key)
            except ValueError:
                logger.warning(
                    f"Corruption in cached index: {key}",
                    cached_metadata=index_dict[key],
                )
                continue

            _index[otu_id] = OTUKeys(**index_dict[key])

        return _index

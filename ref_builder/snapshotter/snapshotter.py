from collections.abc import Generator, Iterable
from dataclasses import dataclass
from pathlib import Path
from uuid import UUID

import orjson
from structlog import get_logger

from ref_builder.resources import RepoMeta, RepoOTU
from ref_builder.snapshotter.otu import OTUSnapshot

logger = get_logger()


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

    def __repr__(self):
        return (
            f"<OTUMetadata {self.id}: taxid={self.taxid} name={self.name} "
            f"acronym={self.acronym} legacy_id={self.legacy_id}>"
        )


class Snapshotter:
    """Load and cache OTU snapshots."""

    def __init__(self, path: Path) -> None:
        """Create a new snapshotter."""
        self.path = path
        """The path to the snapshot root directory."""

        self._meta_path = self.path / "meta.json"
        """The path to the reconstructed Repo metadata."""

        self._index_path = self.path / "index.json"
        """The path to the index data."""

        _index = self._load_index()

        self._index = _index if _index is not None else self._build_index()
        """The index data of this snapshot index."""

    @classmethod
    def new(cls, path: Path, metadata: RepoMeta) -> "Snapshotter":
        """Create a new snapshot index."""
        path.mkdir(exist_ok=True)

        with open(path / "meta.json", "wb") as f:
            f.write(orjson.dumps(metadata.model_dump()))

        return Snapshotter(path)

    def get_id_by_legacy_id(self, legacy_id: str) -> UUID | None:
        """Get an OTU ID by its legacy ID.

        Returns ``None`` if the legacy ID is not found.

        :param legacy_id: The legacy ID to search for.
        :return: The ID of the OTU with the given legacy ID or ``None``.

        """
        self._update_index()

        index_by_legacy_id = {}

        for otu_id in self._index:
            if (legacy_id := self._index[otu_id].legacy_id) is not None:
                index_by_legacy_id[legacy_id] = otu_id

        return index_by_legacy_id.get(legacy_id)

    def get_id_by_name(self, name: str) -> UUID | None:
        """Get an OTU ID by its name.

        Returns ``None`` if the name is not found.

        :param name: The name to search for.
        :return: The ID of the OTU with the given name or ``None``.

        """
        self._update_index()
        return {self._index[otu_id].name: otu_id for otu_id in self._index}.get(name)

    def get_id_by_taxid(self, taxid: int) -> UUID | None:
        """Get an OTU ID by its taxonomy ID.

        Returns ``None`` if the taxonomy ID is not found.

        :param taxid: The taxonomy ID to search for.
        :return: The ID of the OTU with the given taxonomy ID or ``None``.
        """
        self._update_index()
        return {self._index[otu_id].taxid: otu_id for otu_id in self._index}.get(taxid)

    @property
    def otu_ids(self) -> set[UUID]:
        """A list of OTU ids of snapshots."""
        self._update_index()
        return set(self._index.keys())

    def snapshot(
        self,
        otus: Iterable[RepoOTU],
        at_event: int | None = None,
        indent: bool = False,
    ) -> None:
        """Take a new snapshot."""
        options = orjson.OPT_INDENT_2 if indent else None

        _index = {}

        for otu in otus:
            self.cache_otu(otu, at_event=at_event, options=options)
            metadata = OTUKeys(
                id=otu.id,
                taxid=otu.taxid,
                name=otu.name,
                acronym=otu.acronym,
                legacy_id=otu.legacy_id,
            )
            _index[otu.id] = metadata

        self._index = _index
        self._cache_index()

    def iter_otus(self) -> Generator[RepoOTU, None, None]:
        """Iterate over the OTUs in the snapshot."""
        for otu_id in self.otu_ids:
            yield self.load_by_id(otu_id)

    def cache_otu(
        self,
        otu: "RepoOTU",
        at_event: int | None = None,
        options=None,
    ) -> None:
        """Create a snapshot for a single OTU."""
        logger.debug("Writing a snapshot", otu_id=otu.id)

        otu_snap = OTUSnapshot(self.path / f"{otu.id}")
        otu_snap.cache(otu, at_event, options)

        self._index[otu.id] = OTUKeys.from_otu(otu)
        self._cache_index()

    def load_by_id(self, otu_id: UUID) -> RepoOTU | None:
        """Load the most recently snapshotted form of an OTU by its ID.

        Returns ``None`` if the OTU is not found.

        :param otu_id: The ID of the OTU to load.
        :return: The OTU or ``None`` if it is not found.

        """
        try:
            otu_snap = OTUSnapshot(self.path / f"{otu_id}")
        except FileNotFoundError:
            return None

        return otu_snap.load()

    def load_by_name(self, name: str) -> RepoOTU | None:
        """Load the most recently snapshotted form of an OTU by its name.

        Returns ``None`` if the OTU is not found.

        :param name: The name of the OTU to load.
        :return: The OTU or ``None`` if it is not found.

        """
        if otu_id := self.get_id_by_name(name):
            return self.load_by_id(otu_id)

        return None

    def load_by_taxid(self, taxid: int) -> RepoOTU | None:
        """Load the most recently snapshotted form of an OTU by its taxonomy ID.

        Returns ``None`` if the OTU is not found.

        :param taxid: The taxonomy ID of the OTU to load.
        :return: The OTU or ``None`` if it is not found.

        """
        if otu_id := self.get_id_by_taxid(taxid):
            return self.load_by_id(otu_id)

        return None

    def _build_index(self) -> dict[UUID, OTUKeys]:
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

        return index

    def _cache_index(self) -> None:
        """Cache the index to disk."""
        with open(self._index_path, "wb") as f:
            f.write(
                orjson.dumps(
                    {str(otu_id): self._index[otu_id].dict() for otu_id in self._index},
                ),
            )

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

    def _update_index(self) -> None:
        """Update the index in memory."""
        filename_index = {str(otu_id) for otu_id in self._index}

        for path in self.path.iterdir():
            if not path.is_dir() or path.stem in filename_index:
                continue

            try:
                unlisted_otu_id = UUID(path.stem)
            except ValueError:
                continue

            unindexed_otu = self.load_by_id(unlisted_otu_id)

            self._index[unindexed_otu.id] = OTUKeys.from_otu(unindexed_otu)

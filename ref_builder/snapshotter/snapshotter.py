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
    def from_otu(cls, otu: RepoOTU):
        return OTUKeys(
            id=otu.id,
            taxid=otu.taxid,
            name=otu.name,
            acronym=otu.acronym,
            legacy_id=otu.legacy_id,
        )

    def dict(self):
        return {
            "id": self.id,
            "name": self.name,
            "taxid": self.taxid,
            "acronym": self.acronym,
            "legacy_id": self.legacy_id,
        }

    def __repr__(self):
        return (
            f"<OTUMetadata {self.id}: taxid={self.taxid} name={self.name}"
            + " "
            + f"acronym={self.acronym} legacy_id={self.legacy_id}>"
        )


class Snapshotter:
    """Load and cache OTU snapshots."""

    def __init__(self, path: Path):
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
    def new(cls, path: Path, metadata: RepoMeta):
        """Create a new snapshot index."""
        path.mkdir(exist_ok=True)

        with open(path / "meta.json", "wb") as f:
            f.write(orjson.dumps(metadata.model_dump()))

        return Snapshotter(path)

    @property
    def id_to_taxid(self) -> dict[UUID, int]:
        """A mapping of OTU id to Taxonomy IDs."""
        self._update_index()

        return {otu_id: self._index[otu_id].taxid for otu_id in self._index}

    @property
    def index_by_taxid(self) -> dict[int, UUID]:
        """A mapping of Taxonomy ID to OTU id."""
        self._update_index()

        return {self._index[otu_id].taxid: otu_id for otu_id in self._index}

    @property
    def index_by_name(self) -> dict[str, UUID]:
        """A mapping of OTU organism name to OTU Id"""
        self._update_index()

        return {self._index[otu_id].name: otu_id for otu_id in self._index}

    @property
    def index_by_legacy_id(self) -> dict[str, UUID]:
        """A mapping of legacy Virtool id to OTU UUID"""
        self._update_index()

        index_by_legacy_id = {}
        for otu_id in self._index:
            if (legacy_id := self._index[otu_id].legacy_id) is not None:
                index_by_legacy_id[legacy_id] = otu_id

        return index_by_legacy_id

    @property
    def otu_ids(self) -> set[UUID]:
        """A list of OTU ids of snapshots."""
        self._update_index()

        return set(self._index.keys())

    @property
    def accessions(self) -> set[str]:
        return set(self._get_accession_index().keys())

    def snapshot(
        self,
        otus: Iterable[RepoOTU],
        at_event: int | None = None,
        indent: bool = False,
    ):
        """Take a new snapshot"""
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
        """Iterate over the OTUs in the snapshot"""
        for otu_id in self.otu_ids:
            yield self.load_by_id(otu_id)

    def cache_otu(
        self,
        otu: "RepoOTU",
        at_event: int | None = None,
        options=None,
    ):
        """Snapshots a single OTU"""
        logger.debug(f"Writing a snapshot for {otu.taxid}...")
        otu_snap = OTUSnapshot(self.path / f"{otu.id}")
        otu_snap.cache(otu, at_event, options)

        self._index[otu.id] = OTUKeys.from_otu(otu)
        self._cache_index()

    def load_by_id(self, otu_id: UUID) -> RepoOTU | None:
        """Loads an OTU from the most recent repo snapshot"""
        try:
            otu_snap = OTUSnapshot(self.path / f"{otu_id}")
        except FileNotFoundError:
            return None

        return otu_snap.load()

    def load_by_name(self, name: str) -> RepoOTU | None:
        """Takes an OTU name and returns an OTU from the most recent snapshot."""
        otu_id = self.index_by_name.get(name)

        if otu_id:
            return self.load_by_id(otu_id)

        return None

    def load_by_taxid(self, taxid: int) -> RepoOTU | None:
        """Takes a Taxonomy ID and returns an OTU from the most recent snapshot."""
        otu_id = self.index_by_taxid[taxid]

        if otu_id:
            return self.load_by_id(otu_id)

        return None

    def _build_index(self) -> dict[UUID, OTUKeys]:
        """Build a new index from the contents of the snapshot cache directory"""
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

    def _cache_index(self):
        """Cache the index as a dictionary with OTU ID as the key"""
        dict_index = {str(otu_id): self._index[otu_id].dict() for otu_id in self._index}
        with open(self._index_path, "wb") as f:
            f.write(orjson.dumps(dict_index))

    def _load_index(self) -> dict | None:
        """Load the index from file."""
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

    def _update_index(self):
        """Update the index in memory."""
        filename_index = {str(otu_id) for otu_id in self._index}

        for subpath in self.path.iterdir():
            if not subpath.is_dir() or subpath.stem in filename_index:
                continue

            try:
                unlisted_otu_id = UUID(subpath.stem)
            except ValueError:
                continue

            unindexed_otu = self.load_by_id(unlisted_otu_id)

            self._index[unindexed_otu.id] = OTUKeys.from_otu(unindexed_otu)

    def _get_accession_index(self):
        """Return a mapping of all accessions in the snapshot to their parent OTU id."""
        accession_dict = {}
        for otu in self.iter_otus():
            for accession in otu.accessions:
                accession_dict[accession] = otu.id

        return accession_dict

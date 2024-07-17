import shutil
from pathlib import Path
from uuid import UUID

from pydantic import ValidationError

from ref_builder.resources import (
    RepoIsolate,
    RepoOTU,
    RepoSequence,
)
from ref_builder.snapshotter.models import (
    OTUSnapshotIsolate,
    OTUSnapshotMeta,
    OTUSnapshotOTU,
    OTUSnapshotSequence,
    OTUSnapshotToCIsolate,
    toc_adapter,
)


class OTUSnapshotToC:
    """Manages the building and loading of the table of contents."""

    def __init__(self, path: Path) -> None:
        """Initialize the table of contents."""
        self.path: Path = path
        """The path to the table of contents."""

    @staticmethod
    def generate_from_otu(
        otu: "RepoOTU",
    ) -> dict[str, OTUSnapshotToCIsolate]:
        """Return a new table of contents from an OTU."""
        toc = {}
        for isolate in otu.isolates:
            toc[f"{isolate.name}"] = OTUSnapshotToC._generate_table_from_isolate(
                isolate,
            )

        return toc

    def load(self) -> dict[str, OTUSnapshotToCIsolate] | None:
        """Load a table of contents from file."""
        if not self.path.exists():
            return None

        with open(self.path, "rb") as f:
            return toc_adapter.validate_json(f.read())

    def write(
        self,
        data: dict[str, OTUSnapshotToCIsolate],
        indent: int | None = None,
    ) -> None:
        """Write a table of contents to file."""
        with open(self.path, "wb") as f:
            f.write(toc_adapter.dump_json(data, indent=indent))

    def add_sequence(
        self,
        sequence: RepoSequence,
        isolate_id: UUID,
        indent: int | None = None,
    ):
        """Add a new sequence to the table of contents."""
        toc = self.load()

        for key in toc:
            if toc[key].id == isolate_id:
                toc[key].accessions[sequence.accession] = sequence.id
                break

        self.write(toc, indent=indent)

    def add_isolate(self, isolate: RepoIsolate, indent: int | None = None):
        """Add a new isolate to the table of contents."""
        toc = self.load()
        toc[f"{isolate.name}"] = self._generate_table_from_isolate(isolate)
        self.write(data=toc, indent=indent)

    @staticmethod
    def _generate_table_from_isolate(
        isolate: RepoIsolate,
    ) -> OTUSnapshotToCIsolate:
        """Take an isolate and return a table of contents listing for it."""
        return OTUSnapshotToCIsolate(
            id=isolate.id,
            accessions={
                accession: isolate.get_sequence_by_accession(accession).id
                for accession in sorted(isolate.accessions)
            },
        )


class OTUSnapshotDataStore:
    """Stores and retrieves OTU data in snapshot models."""

    def __init__(self, path: Path):
        if not path.exists():
            path.mkdir()

        self.path = path
        """The path to the snapshot's data store directory."""

    @property
    def contents(self):
        """A list of the data store's contents."""
        return list(self.path.glob("*.json"))

    def clean(self):
        """Delete and remake the data store directory."""
        shutil.rmtree(self.path)
        self.path.mkdir()

    def load_isolate(self, isolate_id: UUID) -> OTUSnapshotIsolate:
        """Load and parse an isolate from the data store."""
        with open(self.path / f"{isolate_id}.json", "rb") as f:
            return OTUSnapshotIsolate.model_validate_json(f.read())

    def cache_isolate(
        self,
        isolate: RepoIsolate,
        indent: int | None = None,
    ):
        """Serialize and cache an isolate to the data store."""
        validated_isolate = OTUSnapshotIsolate(**isolate.dict(exclude_contents=True))
        with open(self.path / f"{isolate.id}.json", "w") as f:
            f.write(validated_isolate.model_dump_json(indent=indent))

    def load_sequence(self, sequence_id: UUID) -> OTUSnapshotSequence:
        """Load and parse a sequence from the data store."""
        with open(self.path / f"{sequence_id}.json", "rb") as f:
            return OTUSnapshotSequence.model_validate_json(f.read())

    def cache_sequence(
        self,
        sequence: RepoSequence,
        indent: int | None = None,
    ):
        """Serialize and cache a sequence to the data store."""
        validated_sequence = OTUSnapshotSequence(**sequence.dict())
        with open(self.path / f"{sequence.id}.json", "w") as f:
            f.write(validated_sequence.model_dump_json(indent=indent))


class OTUSnapshot:
    """Manages snapshot data for a single OTU."""

    def __init__(self, path: Path):
        if not path.exists():
            path.mkdir()

        self.path = path
        """The path of this snapshot's directory."""

        self._data = OTUSnapshotDataStore(self.path / "data")
        """The data store of this snapshot. Holds isolate and sequence data."""

        self._toc = OTUSnapshotToC(self.path / "toc.json")
        """The path to this snapshot's table of contents."""

        self._metadata_path = self.path / "metadata.json"
        """The path of this snapshot's metadata file."""

        self._metadata = (
            self._load_metadata() if self._metadata_path.exists() else OTUSnapshotMeta()
        )
        """The metadata for this snapshot."""

    @property
    def at_event(self) -> int | None:
        """The event at which the snapshot was created."""
        return self._metadata.at_event

    @property
    def _otu_path(self):
        """The path to the OTU's taxonomy data."""
        return self.path / "otu.json"

    @property
    def _toc_path(self):
        """The path to the OTU's table of contents."""
        return self.path / "toc.json"

    def clean(self):
        """Delete and remake OTUSnapshot directory structure."""
        shutil.rmtree(self.path)
        self.path.mkdir()
        self._data.clean()

    def cache(
        self,
        otu: "RepoOTU",
        at_event: int | None = None,
        indent: int | None = None,
    ):
        """Cache an OTU at a given event."""
        self._metadata.at_event = at_event

        validated_otu = OTUSnapshotOTU(**otu.dict(exclude_contents=True))
        with open(self._otu_path, "w") as f:
            f.write(validated_otu.model_dump_json(indent=indent))

        for isolate in otu.isolates:
            self._data.cache_isolate(isolate)

            for sequence in isolate.sequences:
                self._data.cache_sequence(sequence)

        self._toc.write(data=OTUSnapshotToC.generate_from_otu(otu), indent=indent)

        self._write_metadata(indent=None)

    def load(self) -> "RepoOTU":
        """Load an OTU from the snapshot."""
        with open(self._otu_path, "rb") as f:
            otu_structure = OTUSnapshotOTU.model_validate_json(f.read())

        toc = self._toc.load()

        isolates = []
        for key in toc:
            isolate_entry = toc[key]

            isolate_structure = self._data.load_isolate(isolate_entry.id)

            sequences = []

            for accession in toc[key].accessions:
                sequence_id = toc[key].accessions[accession]
                sequence_structure = self._data.load_sequence(sequence_id)

                sequence = RepoSequence(**sequence_structure.model_dump())

                sequences.append(sequence)

            isolate_dict = isolate_structure.model_dump()
            isolate_dict["uuid"] = isolate_dict.pop("id")

            isolate = RepoIsolate(**isolate_dict, sequences=sequences)

            isolates.append(isolate)

        otu_dict = otu_structure.model_dump(by_alias=True)
        otu_dict["uuid"] = otu_dict.pop("id")

        return RepoOTU(**otu_dict, isolates=isolates)

    def _write_metadata(self, indent: int | None = None) -> None:
        """Write the snapshot's metadata to file."""
        with open(self._metadata_path, "w") as f:
            f.write(self._metadata.model_dump_json(indent=indent))

    def _load_metadata(self) -> OTUSnapshotMeta | None:
        """Load the snapshot's metadata from file."""
        try:
            with open(self._metadata_path, "rb") as f:
                return OTUSnapshotMeta.model_validate_json(f.read())
        except (FileNotFoundError, ValidationError):
            return None

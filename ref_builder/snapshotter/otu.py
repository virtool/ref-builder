from pathlib import Path
from uuid import UUID

import orjson
from structlog import get_logger

from ref_builder.resources import (
    IsolateModel,
    OTUModel,
    RepoIsolate,
    RepoOTU,
    RepoSequence,
)
from ref_builder.snapshotter.models import OTUSnapshotToCIsolate, toc_adapter

logger = get_logger("snapshotter.otu")


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
            toc[str(isolate.id)] = OTUSnapshotToC._generate_table_from_isolate(
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

    def add_isolate(self, isolate: RepoIsolate, indent: int | None = None) -> None:
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
            name=isolate.name,
            accessions={
                accession: isolate.get_sequence_by_accession(accession).id
                for accession in sorted(isolate.accessions)
            },
        )


class OTUSnapshotDataStore:
    """Stores and retrieves OTU data in snapshot models."""

    def __init__(self, path: Path) -> None:
        """Initialize the data store."""
        if not path.exists():
            path.mkdir()

        self.path = path
        """The path to the snapshot's data store directory."""

    @property
    def contents(self) -> list[Path]:
        """A list of the data store's contents."""
        return list(self.path.glob("*.json"))

    def load_isolate(self, isolate_id: UUID) -> IsolateModel:
        """Load and parse an isolate from the data store."""
        with open(self.path / f"{isolate_id}.json", "rb") as f:
            return IsolateModel.model_validate_json(f.read())

    def cache_isolate(
        self,
        isolate: RepoIsolate,
        indent: int | None = None,
    ) -> None:
        """Serialize and cache an isolate to the data store."""
        isolate_metadata = IsolateModel(
            **isolate.model_dump(exclude={"sequences"}),
        )
        with open(self.path / f"{isolate.id}.json", "w") as f:
            f.write(isolate_metadata.model_dump_json(indent=indent))

    def load_sequence(self, sequence_id: UUID) -> RepoSequence:
        """Load and parse a sequence from the data store."""
        with open(self.path / f"{sequence_id}.json", "rb") as f:
            return RepoSequence.model_validate_json(f.read())

    def cache_sequence(
        self,
        sequence: RepoSequence,
        indent: int | None = None,
    ) -> None:
        """Serialize and cache a sequence to the data store."""
        with open(self.path / f"{sequence.id}.json", "w") as f:
            f.write(sequence.model_dump_json(indent=indent))


class OTUSnapshotter:
    """Manages snapshot data for a single OTU."""

    def __init__(self, path: Path) -> None:
        """Initialize the snapshot."""
        if not path.exists():
            path.mkdir()

        self.path = path
        """The path of this snapshot's directory."""

        self._data = OTUSnapshotDataStore(self.path / "data")
        """The data store of this snapshot. Holds isolate and sequence data."""

        self._toc = OTUSnapshotToC(self.path / "toc.json")
        """The path to this snapshot's table of contents."""

        self.at_event: int | None = None

    @property
    def _otu_path(self) -> Path:
        """The path to the OTU's taxonomy data."""
        return self.path / "otu.json"

    def cache(
        self,
        otu: "RepoOTU",
        at_event: int,
    ) -> None:
        """Cache an OTU at a given event."""
        self.at_event = at_event

        otu_metadata = OTUModel(**otu.model_dump())

        with open(self._otu_path, "wb") as f:
            f.write(
                orjson.dumps(
                    {
                        "at_event": at_event,
                        "data": otu_metadata.model_dump(mode="json")
                    },
                ),
            )

        for isolate in otu.isolates:
            self._data.cache_isolate(isolate)

            for sequence in isolate.sequences:
                self._data.cache_sequence(sequence)

        self._toc.write(data=OTUSnapshotToC.generate_from_otu(otu))

    def load(self) -> "RepoOTU":
        """Load an OTU from the snapshot."""
        with open(self._otu_path, "rb") as f:
            data = orjson.loads(f.read())

            self.at_event = data["at_event"]
            otu_dict = data["data"]

        toc = self._toc.load()

        isolates = []

        for key in toc:
            isolate_entry = toc[key]

            isolate_metadata = self._data.load_isolate(isolate_entry.id)

            sequences = []

            for accession in toc[key].accessions:
                sequence_id = toc[key].accessions[accession]

                sequences.append(self._data.load_sequence(sequence_id))

            isolates.append(
                RepoIsolate(
                    **isolate_metadata.model_dump(),
                    sequences=sequences,
                ),
            )

        return RepoOTU(
            **otu_dict,
            isolates=isolates,
        )

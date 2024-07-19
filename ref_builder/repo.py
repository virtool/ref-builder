"""An access layer for a Virtool event-sourced reference repository.

This is a work in progress.

**Check that OTUs have only one representative (default) isolate.**

The default isolate is set using `rep_isolate` in the `CreateOTU` event. Therefore only
one isolate in an OTU can be default.

TODO: Check if excluded accessions exist in the repo.
TODO: Check for accessions filed in wrong isolates.
TODO: Check for accession conflicts.

"""

import shutil
import uuid
from collections import defaultdict
from collections.abc import Generator
from pathlib import Path

import arrow
from orjson import orjson
from structlog import get_logger

from ref_builder.events import (
    CreateIsolate,
    CreateIsolateData,
    CreateOTU,
    CreateOTUData,
    CreateRepo,
    CreateRepoData,
    CreateSchema,
    CreateSchemaData,
    CreateSequence,
    CreateSequenceData,
    Event,
    EventData,
    EventQuery,
    ExcludeAccession,
    ExcludeAccessionData,
    IsolateQuery,
    OTUQuery,
    RepoQuery,
    SequenceQuery,
)
from ref_builder.index import EventIndex, EventIndexError
from ref_builder.models import Molecule
from ref_builder.resources import (
    RepoIsolate,
    RepoMeta,
    RepoOTU,
    RepoSequence,
)
from ref_builder.schema import OTUSchema, Segment
from ref_builder.snapshotter.snapshotter import Snapshotter
from ref_builder.utils import DataType, IsolateName, IsolateNameType, pad_zeroes

logger = get_logger("repo")


OTU_EVENT_TYPES = (
    CreateOTU,
    CreateIsolate,
    CreateSchema,
    CreateSequence,
    ExcludeAccession,
)


class Repo:
    """An event-sourced repository."""

    def __init__(self, path: Path) -> None:
        self.path = path
        """The path to the repo directory."""

        self.cache_path = self.path / ".cache"
        """The path to the cache subdirectory."""

        self._event_store = EventStore(self.path / "src")
        """The event store of the event sourced repository."""

        self._event_index = EventIndex(self.cache_path / "index")
        """The event index cache of the event sourced repository."""

        snapshot_path = path / ".cache/snapshot"

        self._snapshotter = (
            Snapshotter(path=snapshot_path)
            if snapshot_path.exists()
            else Snapshotter.new(path=snapshot_path, metadata=self.meta)
        )
        """The snapshot index. Maintains and caches the read model of the Repo."""

        # Take a new snapshot if no existing data is found.
        if not self._snapshotter.otu_ids:
            logger.debug("No snapshot data found. Building new snapshot...")
            self.snapshot()

    @classmethod
    def new(cls, data_type: DataType, name: str, path: Path, organism: str):
        """Create a new reference repository."""
        if path.is_file():
            raise ValueError("The target path is a file")

        path.mkdir(parents=True, exist_ok=True)

        if any(path.iterdir()):
            raise ValueError("The target path is not empty")

        with open(path / ".gitignore", "w") as f:
            f.write(".cache\n")

        shutil.copytree(
            Path(__file__).parent.parent / "assets/github",
            path / ".github",
        )

        (path / ".cache").mkdir()

        repo_id = uuid.uuid4()

        _src = EventStore(path / "src")
        _src.write_event(
            CreateRepo,
            CreateRepoData(
                id=repo_id,
                data_type=data_type,
                name=name,
                organism=organism,
            ),
            RepoQuery(repository_id=repo_id),
        )

        return Repo(path)

    @property
    def last_id(self):
        """The id of the most recently added event in the event store."""
        return self._event_store.last_id

    @property
    def meta(self):
        """The metadata for the repository."""
        for event in self._event_store.iter_events():
            if isinstance(event, CreateRepo):
                repo = event.data.model_dump()
                return RepoMeta(**repo, created_at=event.timestamp)

        raise ValueError("No repository creation event found")

    def snapshot(self):
        """Create a snapshot using all the OTUs in the event store."""
        self._snapshotter.snapshot(
            self.iter_otus(ignore_cache=True),
            at_event=self.last_id,
            indent=True,
        )

    def iter_otus(self, ignore_cache: bool = False) -> Generator[RepoOTU, None, None]:
        """Iterate over the OTUs in the snapshot."""
        if ignore_cache:
            event_index = defaultdict(list)

            for event in self._event_store.iter_events():
                otu_id = event.query.model_dump().get("otu_id")

                if otu_id is not None:
                    event_index[otu_id].append(event.id)
        else:
            event_index = self._event_index.load()

        for otu_id in event_index:
            yield self.get_otu(otu_id)

    def create_otu(
        self,
        acronym: str,
        legacy_id: str | None,
        name: str,
        schema: [],
        taxid: int,
    ):
        """Create an OTU."""
        if otu := self.get_otu_by_taxid(taxid):
            raise ValueError(
                f"OTU already exists as {otu}",
            )

        if name in self._snapshotter.index_by_name:
            raise ValueError(f"An OTU with the name '{name}' already exists")

        if legacy_id in self._snapshotter.index_by_legacy_id:
            raise ValueError(f"An OTU with the legacy ID '{legacy_id}' already exists")

        otu_logger = logger.bind(taxid=taxid, name=name, legacy_id=legacy_id)
        otu_logger.info(f"Creating new OTU for Taxonomy ID {taxid}...")

        otu_id = uuid.uuid4()

        event = self._event_store.write_event(
            CreateOTU,
            CreateOTUData(
                id=otu_id,
                acronym=acronym,
                legacy_id=legacy_id,
                name=name,
                schema=schema,
                rep_isolate=None,
                taxid=taxid,
            ),
            OTUQuery(otu_id=otu_id),
        )

        otu_logger.debug("OTU written", event_id=event.id, otu_id=str(otu_id))

        otu = self.get_otu(otu_id)

        self._snapshotter.cache_otu(otu, at_event=self.last_id)

        return otu

    def create_isolate(
        self,
        otu_id: uuid.UUID,
        legacy_id: str | None,
        source_name: str,
        source_type: IsolateNameType,
    ) -> RepoIsolate | None:
        """Create and return a new isolate within the given OTU.
        If the isolate name already exists, return None.
        """
        otu = self.get_otu(otu_id)

        name = IsolateName(type=source_type, value=source_name)
        if otu.get_isolate_id_by_name(name) is not None:
            logger.warning(
                "An isolate by this name already exists",
                isolate_name=str(name),
            )
            return None

        isolate_id = uuid.uuid4()

        event = self._event_store.write_event(
            CreateIsolate,
            CreateIsolateData(id=isolate_id, legacy_id=legacy_id, name=name),
            IsolateQuery(isolate_id=isolate_id, otu_id=otu_id),
        )

        logger.debug(
            "Isolate written",
            event_id=event.id,
            isolate_id=str(isolate_id),
            name=str(name),
        )

        isolate = RepoIsolate(
            uuid=isolate_id,
            legacy_id=legacy_id,
            name=name,
            sequences=[],
        )

        otu.add_isolate(isolate)

        self._snapshotter.cache_otu(otu, at_event=self.last_id)

        return isolate

    def create_sequence(
        self,
        otu_id: uuid.UUID,
        isolate_id: uuid.UUID,
        accession: str,
        definition: str,
        legacy_id: str | None,
        segment: str,
        sequence: str,
    ) -> RepoSequence | None:
        """Create and return a new sequence within the given OTU.
        If the accession already exists in this OTU, return None.
        """
        otu = self.get_otu(otu_id)

        if accession in otu.accessions:
            logger.warning(
                "This accession already exists in the OTU.",
                accession=accession,
                otu_id=str(otu_id),
            )
            return None

        sequence_id = uuid.uuid4()

        event = self._event_store.write_event(
            CreateSequence,
            CreateSequenceData(
                id=sequence_id,
                accession=accession,
                definition=definition,
                legacy_id=legacy_id,
                segment=segment,
                sequence=sequence,
            ),
            SequenceQuery(
                otu_id=otu_id,
                isolate_id=isolate_id,
                sequence_id=sequence_id,
            ),
        )

        logger.debug(
            "Sequence written",
            event_id=event.id,
            sequence_id=str(sequence_id),
            accession=accession,
        )

        sequence = RepoSequence(
            id=sequence_id,
            accession=accession,
            definition=definition,
            legacy_id=legacy_id,
            segment=segment,
            sequence=sequence,
        )

        otu.add_sequence(sequence, isolate_id)

        self._snapshotter.cache_otu(otu, at_event=self.last_id)

        return sequence

    def create_schema(
        self,
        otu_id: uuid.UUID,
        molecule: Molecule,
        segments: list[Segment],
    ):
        otu = self.get_otu(otu_id)

        schema_data = {"molecule": molecule, "segments": segments}

        self._event_store.write_event(
            CreateSchema,
            CreateSchemaData(**schema_data),
            OTUQuery(otu_id=otu.id),
        )

        otu.schema = OTUSchema(**schema_data)

        self._snapshotter.cache_otu(otu, at_event=self.last_id)

        return otu.schema

    def exclude_accession(self, otu_id: uuid.UUID, accession: str):
        """Exclude an accession for an OTU.

        This accession will not be allowed in the repository in the future.

        :param otu_id: the id of the OTU
        :param accession: the accession to exclude

        """
        self._event_store.write_event(
            ExcludeAccession,
            ExcludeAccessionData(accession=accession),
            OTUQuery(otu_id=otu_id),
        )

    def get_otu(self, otu_id: uuid.UUID) -> RepoOTU | None:
        """Return an OTU corresponding with a given OTU Id if it exists, else None."""
        if event_ids := self._get_otu_event_ids(otu_id):
            return self._rehydrate_otu(event_ids)

        return None

    def get_otu_by_taxid(self, taxid: int) -> RepoOTU | None:
        """Return an OTU corresponding with a given OTU Id if it exists, else None"""
        if (otu_id := self._snapshotter.index_by_taxid.get(taxid)) is not None:
            otu = self.get_otu(otu_id)
            return otu

        return None

    def _rehydrate_otu(self, event_ids: list[int]) -> RepoOTU:
        """Rebuild OTU data from a list of event IDs."""
        first_event_id = sorted(event_ids)[0]

        event = self._event_store.read_event(first_event_id)

        if not isinstance(event, CreateOTU):
            raise ValueError(
                f"The first event ({first_event_id}) for an OTU is not a CreateOTU "
                "event",
            )

        otu = RepoOTU(
            uuid=event.data.id,
            acronym=event.data.acronym,
            excluded_accessions=[],
            legacy_id=event.data.legacy_id,
            name=event.data.name,
            taxid=event.data.taxid,
            schema=event.data.otu_schema,
        )

        for event_id in event_ids[1:]:
            event = self._event_store.read_event(event_id)

            if isinstance(event, CreateSchema):
                otu.schema = OTUSchema(
                    molecule=event.data.molecule,
                    segments=event.data.segments,
                )

            elif isinstance(event, CreateIsolate):
                otu.add_isolate(
                    RepoIsolate(
                        uuid=event.data.id,
                        legacy_id=event.data.legacy_id,
                        name=event.data.name,
                    ),
                )

            elif isinstance(event, ExcludeAccession):
                otu.excluded_accessions.add(event.data.accession)

            elif isinstance(event, CreateSequence):
                for isolate in otu.isolates:
                    if isolate.id == event.query.isolate_id:
                        isolate.add_sequence(
                            RepoSequence(
                                id=event.data.id,
                                accession=event.data.accession,
                                definition=event.data.definition,
                                legacy_id=event.data.legacy_id,
                                segment=event.data.segment,
                                sequence=event.data.sequence,
                            ),
                        )

        otu.isolates.sort(key=lambda i: f"{i.name.type} {i.name.value}")

        for isolate in otu.isolates:
            isolate.sequences.sort(key=lambda s: s.accession)

        return otu

    def _get_otu_event_ids(self, otu_id: uuid.UUID) -> list[int]:
        """Gets a list of all event IDs associated with ``otu_id``."""
        event_ids = []

        # Start at zero if no events are indexed.
        at_event = 0

        if indexed := self._event_index.get(otu_id):
            if indexed.at_event == self.last_id:
                return indexed.event_ids

            if indexed.at_event > self.last_id:
                raise EventIndexError("Event index is ahead of the event store.")

            at_event = indexed.at_event
            event_ids = indexed.event_ids

        for event in self._event_store.iter_events(start=at_event + 1):
            if (
                type(event) in OTU_EVENT_TYPES
                and event.query.model_dump().get("otu_id") == otu_id
            ):
                if event.id in event_ids:
                    raise ValueError("Event ID already in event list.")

                event_ids.append(event.id)

        event_ids.sort()

        self._event_index.set(
            otu_id,
            event_ids,
            self.last_id,
        )

        return event_ids


class EventStore:
    """Interface for the event store"""

    def __init__(self, path: Path):
        path.mkdir(exist_ok=True)

        self.path = path
        """The path to the event store directory."""

        self.last_id = 0
        """The id of the latest event."""

        # Check that all events are present and set .last_id to the latest event.
        for event_id in self.event_ids:
            if event_id - self.last_id != 1:
                raise ValueError("Event IDs are not sequential.")

            self.last_id = event_id

    @property
    def event_ids(self) -> list:
        event_ids = []

        for event_path in self.path.iterdir():
            try:
                event_ids.append(int(event_path.stem))
            except ValueError:
                continue

        event_ids.sort()

        return event_ids

    def iter_events(
        self,
        start: int = 1,
    ) -> Generator[Event, None, None]:
        """"""
        if start < 1:
            raise IndexError("Start event ID cannot be less than 1")

        if start not in [int(path.stem) for path in self.path.glob("*.json")]:
            raise IndexError(f"Event {start} not found in event store")

        for path in sorted(self.path.glob("*.json")):
            if path.stem == "meta":
                continue

            if int(path.stem) >= start:
                yield EventStore._read_event_at_path(path)

    def read_event(self, event_id: int) -> Event:
        return EventStore._read_event_at_path(
            self.path / f"{pad_zeroes(event_id)}.json",
        )

    def write_event(
        self,
        cls: type[Event],
        data: EventData,
        query: EventQuery,
    ) -> Event:
        """Write a new event to the repository."""
        event_id = self.last_id + 1

        event = cls(
            id=event_id,
            data=data,
            query=query,
            timestamp=arrow.utcnow().naive,
        )

        with open(self.path / f"{pad_zeroes(event_id)}.json", "wb") as f:
            f.write(
                orjson.dumps(
                    event.model_dump(by_alias=True),
                    f,
                ),
            )

        self.last_id = event_id

        return event

    @staticmethod
    def _read_event_at_path(path: Path) -> Event:
        with open(path, "rb") as f:
            loaded = orjson.loads(f.read())

            match loaded["type"]:
                case "CreateRepo":
                    return CreateRepo(**loaded)
                case "CreateOTU":
                    return CreateOTU(**loaded)
                case "CreateIsolate":
                    return CreateIsolate(**loaded)
                case "CreateSequence":
                    return CreateSequence(**loaded)
                case "CreateSchema":
                    return CreateSchema(**loaded)
                case "ExcludeAccession":
                    return ExcludeAccession(**loaded)

            raise ValueError(f"Unknown event type: {loaded['type']}")

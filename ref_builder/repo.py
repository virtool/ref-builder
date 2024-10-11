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
from collections.abc import Generator, Iterator
from pathlib import Path

import arrow
from orjson import orjson
from structlog import get_logger

from ref_builder.events.base import (
    Event,
    EventData,
    EventQuery,
    RepoQuery,
    OTUQuery,
    IsolateQuery,
    SequenceQuery,
)
from ref_builder.events.repo import (
    CreateRepo,
    CreateRepoData,
)
from ref_builder.events.otu import (
    CreateOTU,
    CreateOTUData,
    CreatePlan,
    CreatePlanData,
    ExcludeAccession,
    ExcludeAccessionData,
    SetReprIsolate,
    SetReprIsolateData,
)
from ref_builder.events.isolate import (
    CreateIsolate,
    CreateIsolateData,
    DeleteIsolate,
    DeleteIsolateData,
)
from ref_builder.events.sequence import (
    CreateSequence,
    CreateSequenceData,
    DeleteSequence,
    DeleteSequenceData,
)
from ref_builder.index import Index
from ref_builder.models import Molecule, OTUMinimal
from ref_builder.plan import MonopartitePlan, MultipartitePlan
from ref_builder.resources import (
    RepoIsolate,
    RepoMeta,
    RepoOTU,
    RepoSequence,
    RepoSettings,
)
from ref_builder.utils import (
    Accession,
    DataType,
    IsolateName,
    pad_zeroes,
)

logger = get_logger("repo")


class Repo:
    """An event-sourced repository."""

    def __init__(self, path: Path) -> None:
        """Create a new instance of the repository."""
        self.path = path
        """The path to the repo directory."""

        self._event_store = EventStore(self.path / "src")
        """The event store of the event sourced repository."""

        self._index = Index(self.path / ".cache" / "index.db")
        """An index for fast lookups.

        It allows fast lookup of OTUs be key fields, fetching of complete OTU state,
        and the events associated with a given OTU ID.
        """

        # Populate the index if it is empty.
        if not self._index.otu_ids:
            logger.info("No index found. Rebuilding...")
            for otu in self.iter_otus_from_events():
                self._index.upsert_otu(otu, self.last_id)

            for event in self._event_store.iter_events():
                try:
                    otu_id = event.query.model_dump()["otu_id"]
                except KeyError:
                    continue

                self._index.add_event_id(event.id, otu_id)

    @classmethod
    def new(
        cls,
        data_type: DataType,
        name: str,
        path: Path,
        organism: str,
        default_segment_length_tolerance: float = 0.03,
    ) -> "Repo":
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
                settings=RepoSettings(
                    default_segment_length_tolerance=default_segment_length_tolerance
                ),
            ),
            RepoQuery(repository_id=repo_id),
        )

        return Repo(path)

    @property
    def last_id(self) -> int:
        """The id of the most recently added event in the event store."""
        return self._event_store.last_id

    @property
    def meta(self) -> RepoMeta:
        """The metadata for the repository."""
        for event in self._event_store.iter_events():
            if isinstance(event, CreateRepo):
                return RepoMeta(**event.data.model_dump(), created_at=event.timestamp)

        raise ValueError("No repository creation event found")

    @property
    def settings(self) -> RepoSettings:
        """The settings for the repository."""
        for event in self._event_store.iter_events():
            if isinstance(event, CreateRepo):
                return event.data.settings

        raise ValueError("No repository creation event found")

    def iter_minimal_otus(self) -> Iterator[OTUMinimal]:
        """Iterate over minimal representations of the OTUs in the repository.

        This is more performant than iterating over full OTUs.

        """
        return self._index.iter_minimal_otus()

    def iter_otus(self) -> Iterator[RepoOTU]:
        """Iterate over the OTUs in the repository."""
        for otu_id in self._index.otu_ids:
            yield self.get_otu(otu_id)

    def iter_otus_from_events(self) -> Iterator[RepoOTU]:
        """Iterate over the OTUs, bypassing the index."""
        events_by_otu = defaultdict(list)

        for event in self._event_store.iter_events():
            otu_id = event.query.model_dump().get("otu_id")

            if otu_id is not None:
                events_by_otu[otu_id].append(event.id)

        for otu_id in events_by_otu:
            yield self._rehydrate_otu(events_by_otu[otu_id])

    def create_otu(
        self,
        acronym: str,
        legacy_id: str | None,
        molecule: Molecule,
        name: str,
        plan: MonopartitePlan | MultipartitePlan,
        taxid: int,
    ) -> RepoOTU | None:
        """Create an OTU."""
        if (otu_id := self.get_otu_id_by_taxid(taxid)) is not None:
            otu = self.get_otu(otu_id)
            raise ValueError(
                f"OTU already exists as {otu}",
            )

        if self._index.get_id_by_name(name):
            raise ValueError(f"An OTU with the name '{name}' already exists")

        if legacy_id and self._index.get_id_by_legacy_id(legacy_id):
            raise ValueError(f"An OTU with the legacy ID '{legacy_id}' already exists")

        otu_logger = logger.bind(taxid=taxid, name=name, legacy_id=legacy_id)
        otu_logger.info(f"Creating new OTU for Taxonomy ID {taxid}...")

        otu_id = uuid.uuid4()

        event = self._write_event(
            CreateOTU,
            CreateOTUData(
                id=otu_id,
                acronym=acronym,
                legacy_id=legacy_id,
                molecule=molecule,
                name=name,
                plan=plan,
                taxid=taxid,
            ),
            OTUQuery(otu_id=otu_id),
        )

        otu_logger.debug("OTU written", event_id=event.id, otu_id=str(otu_id))

        return self.get_otu(otu_id)

    def create_isolate(
        self,
        otu_id: uuid.UUID,
        legacy_id: str | None,
        name: IsolateName | None,
    ) -> RepoIsolate | None:
        """Create and isolate for the OTU with ``otu_id``.

        If the isolate name already exists, return None.
        """
        otu = self.get_otu(otu_id)

        if name is not None and otu.get_isolate_id_by_name(name):
            raise ValueError(f"Isolate name already exists: {name}")

        isolate_id = uuid.uuid4()

        event = self._write_event(
            CreateIsolate,
            CreateIsolateData(id=isolate_id, legacy_id=legacy_id, name=name),
            IsolateQuery(isolate_id=isolate_id, otu_id=otu_id),
        )

        logger.debug(
            "Isolate written",
            event_id=event.id,
            isolate_id=str(isolate_id),
            name=str(name) if name is not None else None,
        )

        return self.get_otu(otu_id).get_isolate(isolate_id)

    def delete_isolate(
        self,
        otu_id: uuid.UUID,
        isolate_id: uuid.UUID,
        rationale: str,
    ) -> None:
        """Delete an existing isolate from a given OTU."""
        self._write_event(
            DeleteIsolate,
            DeleteIsolateData(rationale=rationale),
            IsolateQuery(otu_id=otu_id, isolate_id=isolate_id),
        )

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

        versioned_accession = Accession.from_string(accession)

        if versioned_accession.key in otu.accessions:
            extant_sequence = otu.get_sequence_by_accession(versioned_accession.key)

            if extant_sequence is not None:
                if extant_sequence.accession == versioned_accession:
                    logger.warning(
                        "This accession already exists in the OTU.",
                        accession=str(extant_sequence.accession),
                        otu_id=str(otu_id),
                    )
                    return None

                logger.warning(
                    "New version of accession found.",
                    accession=versioned_accession.key,
                )

        sequence_id = uuid.uuid4()

        event = self._write_event(
            CreateSequence,
            CreateSequenceData(
                id=sequence_id,
                accession=versioned_accession,
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
            accession=str(versioned_accession),
        )

        return (
            self.get_otu(otu_id)
            .get_isolate(isolate_id)
            .get_sequence_by_accession(versioned_accession.key)
        )

    def replace_sequence(
        self,
        otu_id: uuid.UUID,
        isolate_id: uuid.UUID,
        accession: str,
        definition: str,
        legacy_id: str | None,
        segment: str,
        sequence: str,
        replaced_sequence_id: uuid.UUID,
        rationale: str,
    ) -> RepoSequence | None:
        """Create a new sequence and delete an existing sequence,
        replacing the old sequence under the isolate.
        """
        new_sequence = self.create_sequence(
            otu_id=otu_id,
            isolate_id=isolate_id,
            accession=accession,
            definition=definition,
            legacy_id=legacy_id,
            segment=segment,
            sequence=sequence,
        )

        if new_sequence is None:
            return None

        self._write_event(
            DeleteSequence,
            DeleteSequenceData(
                replacement=new_sequence.id,
                rationale=rationale,
            ),
            SequenceQuery(
                otu_id=otu_id,
                isolate_id=isolate_id,
                sequence_id=replaced_sequence_id,
            ),
        )

        return new_sequence

    def set_repr_isolate(self, otu_id: uuid.UUID, isolate_id: uuid.UUID) -> uuid.UUID:
        """Set the representative isolate for an OTU."""
        otu = self.get_otu(otu_id)

        self._write_event(
            SetReprIsolate,
            SetReprIsolateData(isolate_id=isolate_id),
            OTUQuery(otu_id=otu.id),
        )

        return self.get_otu(otu_id).repr_isolate

    def exclude_accession(self, otu_id: uuid.UUID, accession: str) -> set:
        """Exclude an accession for an OTU.

        This accession will not be allowed in the repository in the future.

        :param otu_id: the id of the OTU
        :param accession: the accession to exclude

        """
        otu = self.get_otu(otu_id)

        if accession in otu.excluded_accessions:
            logger.debug("Accession is already excluded.", accession=accession)
        else:
            self._write_event(
                ExcludeAccession,
                ExcludeAccessionData(accession=accession),
                OTUQuery(otu_id=otu_id),
            )

        return self.get_otu(otu_id).excluded_accessions

    def get_otu(self, otu_id: uuid.UUID) -> RepoOTU | None:
        """Get the OTU with the given ``otu_id``.

        If the OTU does not exist, ``None`` is returned.

        :param otu_id: the id of the OTU
        :return: the OTU or ``None``

        """
        event_index_item = self._index.get_event_ids_by_otu_id(otu_id)

        if event_index_item is None:
            return None

        otu = self._rehydrate_otu(event_index_item.event_ids)

        self._index.upsert_otu(otu, self.last_id)

        return otu

    def get_otu_by_taxid(self, taxid: int) -> RepoOTU | None:
        """Return the OTU with the given ``taxid``.

        If no OTU is found, return None.

        :param taxid: the taxonomy ID of the OTU
        :return: the matching OTU or ``None``

        """
        if (otu_id := self.get_otu_id_by_taxid(taxid)) is not None:
            return self.get_otu(otu_id)

        return None

    def get_otu_id_by_taxid(self, taxid: int) -> uuid.UUID | None:
        """Return the UUID of the OTU with the given ``taxid``.

        If no OTU is found, return None.

        :param taxid: the taxonomy ID of the OTU
        :return: the UUID of the OTU or ``None``

        """
        return self._index.get_id_by_taxid(taxid)

    def _rehydrate_otu(self, event_ids: list[int]) -> RepoOTU:
        event = self._event_store.read_event(event_ids[0])

        if not isinstance(event, CreateOTU):
            raise TypeError(
                f"The first event ({event_ids[0]}) for an OTU is not a CreateOTU "
                "event",
            )

        otu = RepoOTU(
            id=event.data.id,
            acronym=event.data.acronym,
            excluded_accessions=set(),
            isolates=[],
            legacy_id=event.data.legacy_id,
            molecule=event.data.molecule,
            name=event.data.name,
            repr_isolate=None,
            plan=event.data.plan,
            taxid=event.data.taxid,
        )

        for event_id in event_ids[1:]:
            event = self._event_store.read_event(event_id)

            if isinstance(event, CreatePlan):
                otu.plan = event.data.plan

            elif isinstance(event, SetReprIsolate):
                otu.repr_isolate = event.data.isolate_id

            elif isinstance(event, CreateIsolate):
                otu.add_isolate(
                    RepoIsolate(
                        id=event.data.id,
                        legacy_id=event.data.legacy_id,
                        name=event.data.name,
                        sequences=[],
                    ),
                )

            elif isinstance(event, DeleteIsolate):
                otu.delete_isolate(event.query.isolate_id)

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

            elif isinstance(event, DeleteSequence):
                otu.delete_sequence(
                    event.query.sequence_id,
                    event.query.isolate_id,
                )

        otu.isolates.sort(
            key=lambda i: f"{i.name.type} {i.name.value}"
            if type(i.name) is IsolateName
            else "",
        )

        for isolate in otu.isolates:
            isolate.sequences.sort(key=lambda s: s.accession)

        return otu

    def _write_event(self, cls, data, query):
        event = self._event_store.write_event(cls, data, query)

        if hasattr(query, "otu_id"):
            self._index.add_event_id(event.id, query.otu_id)

        return event


class EventStore:
    """Interface for the event store."""

    def __init__(self, path: Path) -> None:
        """Create a new instance of the event store.

        If no store exists at ``path``, a new store will be created.

        :param path: the path to the event store directory

        """
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
        """Yield all events in the event store.

        Events are yielded by ascending event ID, which corresponds to the order in
        which they were written.

        Optionally, the starting event ID can be specified using the ``start``
        parameter.

        :param start: the event ID to start from
        :return: a generator of events

        """
        if start < 1:
            raise IndexError("Start event ID cannot be less than 1")

        if not Path(self.path / f"{pad_zeroes(start)}.json").exists():
            raise IndexError(f"Event {start} not found in event store")

        for event_id in range(start, self.last_id + 1):
            try:
                yield self.read_event(event_id)
            except FileNotFoundError:
                break

    def read_event(self, event_id: int) -> Event:
        """Read the event with the given ``event_id``.

        :param event_id: the ID of the event to read
        :return: the event

        """
        with open(self.path / f"{pad_zeroes(event_id)}.json", "rb") as f:
            loaded = orjson.loads(f.read())

            try:
                cls = {
                    "CreateRepo": CreateRepo,
                    "CreateOTU": CreateOTU,
                    "CreateIsolate": CreateIsolate,
                    "CreateSequence": CreateSequence,
                    "DeleteIsolate": DeleteIsolate,
                    "DeleteSequence": DeleteSequence,
                    "CreatePlan": CreatePlan,
                    "SetReprIsolate": SetReprIsolate,
                    "ExcludeAccession": ExcludeAccession,
                }[loaded["type"]]

                return cls(**loaded)

            except KeyError:
                raise ValueError(f"Unknown event type: {loaded['type']}")

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

from collections.abc import Generator
from pathlib import Path

from orjson import orjson

from ref_builder.events.base import Event
from ref_builder.events.isolate import (
    CreateIsolate,
    DeleteIsolate,
    LinkSequence,
    UnlinkSequence,
    RenameIsolate,
)
from ref_builder.events.otu import (
    CreateOTU,
    CreatePlan,
    SetRepresentativeIsolate,
    UpdateExcludedAccessions,
)
from ref_builder.events.repo import CreateRepo
from ref_builder.events.sequence import CreateSequence, DeleteSequence
from ref_builder.utils import pad_zeroes


class EventStore:
    """Interface for the event store."""

    def __init__(self, path: Path) -> None:
        """Create a new instance of the event store.

        If no store exists at ``path``, a new store will be created.

        :param path: the path to the event store directory

        """
        self.events_path = path / "src"
        """The path to the event store directory."""

        self.events_path.mkdir(exist_ok=True)

        self.last_id = 0
        """The id of the latest event."""

        # Check that all events are present and set .last_id to the latest event.
        for event_id in self.event_ids:
            if event_id - self.last_id != 1:
                raise ValueError("Event IDs are not sequential.")

            self.last_id = event_id

    @property
    def event_ids(self) -> list[int]:
        event_ids = []

        for event_path in self.events_path.iterdir():
            try:
                event_ids.append(int(event_path.stem))
            except ValueError as e:
                if "does not include a valid version" in str(e):
                    continue

                raise

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

        if not Path(self.events_path / f"{pad_zeroes(start)}.json").exists():
            # Yield no events if ``start`` is out of range.
            return None

        for event_id in range(start, self.last_id + 1):
            try:
                yield self.read_event(event_id)
            except FileNotFoundError:
                break

    def prune(self, event_id: int) -> None:
        """Remove all events after the given ``event_id``.

        This does not update the SQLite index. The removed events should either be
        uncommitted or the index should be rebuilt after pruning.

        :param event_id: the ID of the event to prune after

        """
        for event_path in self.events_path.iterdir():
            if int(event_path.stem) > event_id:
                event_path.unlink()

        self.last_id = event_id

    def read_event(self, event_id: int) -> Event:
        """Read the event with the given ``event_id``.

        :param event_id: the ID of the event to read
        :return: the event

        """
        with open(self.events_path / f"{pad_zeroes(event_id)}.json", "rb") as f:
            loaded = orjson.loads(f.read())

            try:
                cls = {
                    "CreateRepo": CreateRepo,
                    "CreateOTU": CreateOTU,
                    "CreateIsolate": CreateIsolate,
                    "CreateSequence": CreateSequence,
                    "LinkSequence": LinkSequence,
                    "UnlinkSequence": UnlinkSequence,
                    "DeleteIsolate": DeleteIsolate,
                    "DeleteSequence": DeleteSequence,
                    "CreatePlan": CreatePlan,
                    "RenameIsolate": RenameIsolate,
                    "SetRepresentativeIsolate": SetRepresentativeIsolate,
                    "UpdateExcludedAccessions": UpdateExcludedAccessions,
                }[loaded["type"]]

                return cls(**loaded)

            except KeyError:
                raise ValueError(f"Unknown event type: {loaded['type']}")

    def write_event(self, event: Event) -> Event:
        """Write a new event to the repository."""
        with open(self.events_path / f"{pad_zeroes(event.id)}.json", "wb") as f:
            f.write(
                orjson.dumps(
                    event.model_dump(by_alias=True, mode="json"),
                    f,
                ),
            )

        self.last_id = event.id

        return event

"""An access layer for a Virtool event-sourced reference repository.

This is a work in progress.

**Check that OTUs have only one representative (default) isolate.**

The default isolate is set using `rep_isolate` in the `CreateOTU` event. Therefore only
one isolate in an OTU can be default.

TODO: Check if excluded accessions exist in the repo.
TODO: Check for accessions filed in wrong isolates.
TODO: Check for accession conflicts.

"""

import datetime
import shutil
import uuid
import warnings
from collections import defaultdict
from collections.abc import Collection, Generator, Iterator
from contextlib import contextmanager
from pathlib import Path

import arrow
from structlog import get_logger

from ref_builder.errors import (
    InvalidInputError,
    LockRequiredError,
    TransactionExistsError,
    TransactionRequiredError,
)
from ref_builder.events.base import (
    ApplicableEvent,
    Event,
    EventData,
    EventQuery,
    IsolateQuery,
    OTUQuery,
    RepoQuery,
    SequenceQuery,
)
from ref_builder.events.isolate import (
    CreateIsolate,
    CreateIsolateData,
    DeleteIsolate,
    DeleteIsolateData,
    LinkSequence,
    LinkSequenceData,
    UnlinkSequence,
    UnlinkSequenceData,
)
from ref_builder.events.otu import (
    CreateOTU,
    CreateOTUData,
    CreatePlan,
    CreatePlanData,
    SetRepresentativeIsolate,
    SetRepresentativeIsolateData,
    UpdateExcludedAccessions,
    UpdateExcludedAccessionsData,
)
from ref_builder.events.repo import (
    CreateRepo,
    CreateRepoData,
)
from ref_builder.events.sequence import (
    CreateSequence,
    CreateSequenceData,
    DeleteSequence,
    DeleteSequenceData,
)
from ref_builder.index import Index
from ref_builder.lock import Lock
from ref_builder.models import Molecule, OTUMinimal
from ref_builder.otu.validate import check_otu_is_valid
from ref_builder.plan import Plan
from ref_builder.resources import (
    RepoIsolate,
    RepoMeta,
    RepoOTU,
    RepoSequence,
    RepoSettings,
)
from ref_builder.store import EventStore
from ref_builder.transaction import AbortTransactionError, Transaction
from ref_builder.utils import (
    Accession,
    DataType,
    ExcludedAccessionAction,
    IsolateName,
    get_accession_key,
)

GITIGNORE_CONTENTS = [".cache", "lock"]

logger = get_logger("repo")


class Repo:
    """An event-sourced repository."""

    def __init__(self, path: Path) -> None:
        """Create a new instance of the repository."""
        self.path = path
        """The path to the repo directory."""

        self._event_store = EventStore(self.path)
        """The event store of the event sourced repository."""

        self._index = Index(self.path / ".cache" / "index.db")
        """An index for fast lookups.

        It allows fast lookup of OTUs be key fields, fetching of complete OTU state,
        and the events associated with a given OTU ID.
        """

        self._lock = Lock(self.path)
        """A lock for the repository."""

        self._transaction: Transaction | None = None
        """The current transaction, if one is active."""

        try:
            with open(self.path / "head") as f:
                self._head_id = int(f.read())
        except FileNotFoundError:
            self._head_id = self.last_id

        self._prune()

        # Populate the index if it is empty.
        if not self._index.otu_ids:
            self.rebuild_index()

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
            f.write("\n".join(GITIGNORE_CONTENTS) + "\n")

        shutil.copytree(
            Path(__file__).parent.parent / "assets/github",
            path / ".github",
        )

        (path / ".cache").mkdir()

        repo_id = uuid.uuid4()

        event = EventStore(path).write_event(
            CreateRepo(
                id=1,
                data=CreateRepoData(
                    id=repo_id,
                    data_type=data_type,
                    name=name,
                    organism=organism,
                    settings=RepoSettings(
                        default_segment_length_tolerance=default_segment_length_tolerance
                    ),
                ),
                query=RepoQuery(repository_id=repo_id),
                timestamp=arrow.utcnow().naive,
            )
        )

        with open(path / "head", "w") as f:
            f.write(str(event.id))

        return Repo(path)

    @property
    def head_id(self) -> int:
        """The id of the validated event most recently  added to the repository."""
        return self._head_id

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

    @contextmanager
    def lock(self) -> Iterator[None]:
        """Lock the repository.

        This prevents read and write access from  other ``ref-builder`` processes.
        """
        self._lock.lock()

        try:
            yield
        finally:
            self._lock.unlock()

    @contextmanager
    def use_transaction(self) -> Generator[Transaction, None, None]:
        """Lock the repository during modification.

        This prevents writes from  other ``ref-builder`` processes.
        """
        if self._transaction:
            raise TransactionExistsError

        self.prune()

        if self.last_id != self.head_id:
            logger.error(
                "Head ID and last event ID do not match.",
                head_id=self.head_id,
                last_id=self.last_id,
            )

            raise TransactionExistsError

        if not self._lock.locked:
            raise LockRequiredError

        self._transaction = Transaction()

        try:
            yield self._transaction

            for otu_id in self._transaction.affected_otu_ids:
                if not check_otu_is_valid(self.get_otu(otu_id)):
                    self._transaction.abort()

        except Exception as e:
            logger.debug(
                "Error encountered mid-transaction. Pruning events...",
                head_id=self.head_id,
                last_id=self.last_id,
            )

            self.prune()

            if not isinstance(e, AbortTransactionError):
                raise

        else:
            self._head_id = self.last_id

            with open(self.path / "head", "w") as f:
                f.write(str(self._head_id))

        finally:
            self._transaction = None

    def prune(self) -> None:
        """Prune an events ahead of the head.

        This removes the events from the store and updates the index accordingly.
        """
        self._index.prune(self.head_id)
        self._event_store.prune(self.head_id)

    def clear_index(self) -> bool:
        """Delete and replace the repository read index."""
        index_path = self._index.path

        if index_path.exists():
            index_path.unlink()

            (index_path.parent / f"{index_path.stem}.db-shm").unlink()
            (index_path.parent / f"{index_path.stem}.db-wal").unlink()

            return True

        return False

    def rebuild_index(self) -> None:
        """Rebuild the repository read index."""
        logger.info(
            "No repo index found. Rebuilding...",
            path=str(self.path),
        )

        for otu in self.iter_otus_from_events():
            self._index.upsert_otu(otu, self.last_id)

        for event in self._event_store.iter_events():
            try:
                otu_id = event.query.model_dump()["otu_id"]
            except KeyError:
                continue

            self._index.add_event_id(event.id, otu_id, event.timestamp)

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
        event_ids_by_otu = defaultdict(list)

        for event in self._event_store.iter_events():
            if hasattr(event.query, "otu_id"):
                event_ids_by_otu[event.query.otu_id].append(event.id)

        for otu_id in event_ids_by_otu:
            yield self._rehydrate_otu(
                self._event_store.read_event(event_id)
                for event_id in event_ids_by_otu[otu_id]
            )

    def create_otu(
        self,
        acronym: str,
        legacy_id: str | None,
        molecule: Molecule,
        name: str,
        plan: Plan,
        taxid: int,
    ) -> RepoOTU | None:
        """Create an OTU."""
        if (otu_id := self.get_otu_id_by_taxid(taxid)) is not None:
            otu = self.get_otu(otu_id)
            raise ValueError(f"OTU already exists as {otu.id}")

        if self._index.get_id_by_name(name):
            raise ValueError(f"An OTU with the name '{name}' already exists")

        if legacy_id and self._index.get_id_by_legacy_id(legacy_id):
            raise ValueError(f"An OTU with the legacy ID '{legacy_id}' already exists")

        logger.info("Creating new OTU.", taxid=taxid, name=name)

        otu_id = uuid.uuid4()

        self._write_event(
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
        if (otu_ := self.get_otu(otu_id)) is None:
            raise ValueError(f"OTU does not exist: {otu_id}")

        if isolate_id not in otu_.isolate_ids:
            raise ValueError(f"Isolate does not exist: {isolate_id}")

        if isolate_id == otu_.representative_isolate:
            raise ValueError("Representative isolate cannot be deleted.")

        self._write_event(
            DeleteIsolate,
            DeleteIsolateData(rationale=rationale),
            IsolateQuery(otu_id=otu_id, isolate_id=isolate_id),
        )

    def create_sequence(
        self,
        otu_id: uuid.UUID,
        accession: str,
        definition: str,
        legacy_id: str | None,
        segment: uuid.UUID,
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
                raise ValueError(
                    f"Accession {versioned_accession} already exists in the OTU.",
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
                sequence_id=sequence_id,
            ),
        )

        logger.debug(
            "Sequence written",
            event_id=event.id,
            sequence_id=str(sequence_id),
            accession=str(versioned_accession),
        )

        return self.get_otu(otu_id).get_sequence_by_id(sequence_id)

    def replace_sequence(
        self,
        otu_id: uuid.UUID,
        isolate_id: uuid.UUID,
        sequence_id: uuid.UUID,
        replaced_sequence_id: uuid.UUID,
        rationale: str,
    ) -> RepoSequence | None:
        """Link a new sequence, replacing the old sequence under the isolate.
        Delete the original sequence if it is still in the OTU.
        """
        otu = self.get_otu(otu_id)

        new_sequence = otu.get_sequence_by_id(sequence_id)
        if new_sequence is None:
            raise ValueError(f"New sequence does not exist: {sequence_id}")

        self._write_event(
            UnlinkSequence,
            UnlinkSequenceData(
                sequence_id=replaced_sequence_id,
            ),
            IsolateQuery(otu_id=otu_id, isolate_id=isolate_id),
        )

        self._write_event(
            LinkSequence,
            LinkSequenceData(sequence_id=new_sequence.id),
            IsolateQuery(
                otu_id=otu_id,
                isolate_id=isolate_id,
            ),
        )

        if otu.get_sequence_by_id(replaced_sequence_id) is not None:
            self._write_event(
                DeleteSequence,
                DeleteSequenceData(
                    sequence_id=replaced_sequence_id,
                    replacement=new_sequence.id,
                    rationale=rationale,
                ),
                SequenceQuery(
                    otu_id=otu_id,
                    sequence_id=replaced_sequence_id,
                ),
            )

        return self.get_otu(otu_id).get_sequence_by_id(new_sequence.id)

    def set_plan(self, otu_id: uuid.UUID, plan: Plan) -> Plan:
        """Set the isolate plan for an OTU."""
        self._write_event(
            CreatePlan,
            CreatePlanData(plan=plan),
            OTUQuery(otu_id=otu_id),
        )

        return self.get_otu(otu_id).plan

    def link_sequence(
        self, otu_id: uuid.UUID, isolate_id: uuid.UUID, sequence_id: uuid.UUID
    ) -> RepoSequence | None:
        """Link an existing sequence to an existing isolate."""
        log = logger.bind(otu_id=str(otu_id))

        otu = self.get_otu(otu_id)

        if otu is None:
            log.error("OTU not found.")
            return None

        isolate = otu.get_isolate(isolate_id)

        if isolate is None:
            log.error("Isolate not found.", isolate_id=str(isolate_id))
            return None

        if sequence_id in {s.id for s in isolate.sequences}:
            log.warning(
                "Sequence is already linked to isolate.", sequence_id=str(sequence_id)
            )
            return None

        sequence = otu.get_sequence_by_id(sequence_id)

        if sequence is None:
            log.error("Sequence not found.", sequence_id=str(sequence_id))
            return None

        if sequence.segment in {s.segment for s in isolate.sequences}:
            log.warning(
                "Segment is already linked to isolate.",
                segment_id=str(sequence.segment),
                sequence_id=str(sequence_id),
            )
            return None

        event = self._write_event(
            LinkSequence,
            LinkSequenceData(sequence_id=sequence_id),
            IsolateQuery(
                otu_id=otu_id,
                isolate_id=isolate_id,
            ),
        )

        log.debug(
            "Sequence linked to isolate",
            event_id=event.id,
            sequence_id=str(sequence_id),
            isolate_id=str(isolate_id),
            isolate_name=str(isolate.name),
            accession=str(sequence.accession),
        )

        return (
            self.get_otu(otu_id).get_isolate(isolate_id).get_sequence_by_id(sequence_id)
        )

    def set_representative_isolate(
        self, otu_id: uuid.UUID, isolate_id: uuid.UUID
    ) -> uuid.UUID:
        """Set the representative isolate for an OTU."""
        otu = self.get_otu(otu_id)

        self._write_event(
            SetRepresentativeIsolate,
            SetRepresentativeIsolateData(isolate_id=isolate_id),
            OTUQuery(otu_id=otu.id),
        )

        return self.get_otu(otu_id).representative_isolate

    def exclude_accession(self, otu_id: uuid.UUID, accession: str) -> set:
        """Exclude an accession for an OTU.

        This accession will not be allowed in the repository in the future.

        :param otu_id: the id of the OTU
        :param accession: the accession to exclude

        """
        otu = self.get_otu(otu_id)

        try:
            accession_key = get_accession_key(accession)

        except ValueError as e:
            if "Invalid accession key" in str(e):
                logger.warning(
                    "Invalid accession included in set. "
                    "No changes were made to excluded accessions.",
                    accession=accession,
                )
            return otu.excluded_accessions

        if accession_key in otu.excluded_accessions:
            logger.debug("Accession is already excluded.", accession=accession_key)
        else:
            self._write_event(
                UpdateExcludedAccessions,
                UpdateExcludedAccessionsData(
                    accessions={accession_key},
                    action=ExcludedAccessionAction.EXCLUDE,
                ),
                OTUQuery(otu_id=otu_id),
            )

        return self.get_otu(otu_id).excluded_accessions

    def exclude_accessions(
        self,
        otu_id: uuid.UUID,
        accessions: Collection[str],
    ) -> set[str]:
        """Add accessions to OTU's excluded accessions."""
        otu = self.get_otu(otu_id)

        try:
            excludable_accessions = {
                get_accession_key(raw_accession) for raw_accession in accessions
            }
        except ValueError as e:
            if "Invalid accession key" in str(e):
                logger.warning(
                    "Invalid accession included in set. "
                    "No changes were made to excluded accessions.",
                    accessions=sorted(accessions),
                )

                return otu.excluded_accessions

            raise

        if unremovable_accessions := excludable_accessions.intersection(otu.accessions):
            logger.warning(
                "Accessions currently in OTU cannot be removed.",
                unremovable_accessions=sorted(unremovable_accessions),
            )
            excludable_accessions -= unremovable_accessions

        if (
            extant_requested_accessions := excludable_accessions
            & otu.excluded_accessions
        ):
            logger.info(
                "Ignoring already excluded accessions",
                requested_exclusions=sorted(extant_requested_accessions),
                old_excluded_accessions=sorted(otu.excluded_accessions),
            )

            excludable_accessions -= otu.excluded_accessions

        if excludable_accessions:
            self._write_event(
                UpdateExcludedAccessions,
                UpdateExcludedAccessionsData(
                    accessions=excludable_accessions,
                    action=ExcludedAccessionAction.EXCLUDE,
                ),
                OTUQuery(otu_id=otu_id),
            )

            logger.info(
                "Added accessions to excluded accession list.",
                taxid=otu.taxid,
                otu_id=str(otu.id),
                new_excluded_accessions=sorted(excludable_accessions),
                old_excluded_accessions=sorted(otu.excluded_accessions),
            )
        else:
            logger.warning("No excludable accessions were given.")

        return self.get_otu(otu_id).excluded_accessions

    def allow_accessions(
        self,
        otu_id: uuid.UUID,
        accessions: Collection[str],
    ) -> set[str]:
        """Remove accessions from OTU's excluded accessions."""
        otu = self.get_otu(otu_id)

        allowable_accessions = set(accessions)

        if redundant_accessions := allowable_accessions - otu.excluded_accessions:
            logger.debug(
                "Ignoring non-excluded accessions",
                non_excluded_accessions=sorted(redundant_accessions),
            )

            allowable_accessions = allowable_accessions - redundant_accessions

        if allowable_accessions:
            self._write_event(
                UpdateExcludedAccessions,
                UpdateExcludedAccessionsData(
                    accessions=set(allowable_accessions),
                    action=ExcludedAccessionAction.ALLOW,
                ),
                OTUQuery(otu_id=otu_id),
            )

            logger.info(
                "Removed accessions from excluded accession list.",
                taxid=otu.taxid,
                otu_id=str(otu.id),
                new_excluded_accessions=sorted(allowable_accessions),
            )

        return self.get_otu(otu_id).excluded_accessions

    def get_otu_id_by_isolate_id(self, isolate_id: uuid.UUID) -> uuid.UUID | None:
        """Get an OTU ID from an isolate ID that belongs to it."""
        return self._index.get_id_by_isolate_id(isolate_id)

    def get_otu(self, otu_id: uuid.UUID) -> RepoOTU | None:
        """Get the OTU with the given ``otu_id``.

        If the OTU does not exist, ``None`` is returned.

        :param otu_id: the id of the OTU
        :return: the OTU or ``None``

        """
        event_index_item = self._index.get_event_ids_by_otu_id(otu_id)

        if event_index_item is None:
            return None

        try:
            events = (
                self._event_store.read_event(event_id)
                for event_id in event_index_item.event_ids
            )

            otu = self._rehydrate_otu(events)

        except FileNotFoundError:
            logger.error("Event exists in index, but not in source. Deleting index...")

            self.clear_index()

            raise

        self._index.upsert_otu(otu, self.last_id)

        return otu

    def iter_otu_events(
        self, otu_id: uuid.UUID
    ) -> Generator[ApplicableEvent, None, None]:
        """Iterate through event log."""
        event_index_item = self._index.get_event_ids_by_otu_id(otu_id)

        if event_index_item is not None:
            for event_id in event_index_item.event_ids:
                yield self._event_store.read_event(event_id)

    def iter_event_metadata(self):
        """Iterate through the event metadata of all events."""
        yield from self._index.iter_event_metadata()

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

    def get_otu_id_by_partial(self, partial: str) -> uuid.UUID | None:
        """Return the UUID of the OTU starting with the given ``partial`` string.
        Raise a ValueErrror if more than one matching OTU id is found.

        If no OTU is found, return None.

        :param partial: a partial segment of the OTU id with a minimum length of 8
        :return: the UUID of the OTU or ``None``

        """
        if len(partial) < 8:
            raise InvalidInputError(
                "Partial ID segment must be at least 8 characters long."
            )

        return self._index.get_id_by_partial(partial)

    def get_isolate(self, isolate_id: uuid.UUID) -> uuid.UUID | None:
        """Return the isolate with the given id if it exists, else None."""
        if otu_id := self.get_otu_id_by_isolate_id(isolate_id):
            return self.get_otu(otu_id).get_isolate(isolate_id)

        return None

    def get_isolate_id_by_partial(self, partial: str) -> uuid.UUID | None:
        """Return the UUID of the isolate starting with the given ``partial`` string.
        Raise a ValueErrror if more than one matching isolate id is found.
        If no isolate is found, return None.

        :param partial: a partial segment of the isolate id with a minimum length of 8
        :return: the UUID of the isolate or ``None``

        """
        if len(partial) < 8:
            raise InvalidInputError(
                "Partial ID segment must be at least 8 characters long."
            )

        return self._index.get_isolate_id_by_partial(partial)

    def get_otu_first_created(self, otu_id: uuid.UUID) -> datetime.datetime | None:
        """Get the timestamp of the first event associated with an OTU.
        If no events can be found for this OTU, return None.
        """
        return self._index.get_first_timestamp_by_otu_id(otu_id)

    def get_otu_last_modified(self, otu_id: uuid.UUID) -> datetime.datetime | None:
        """Get the timestamp of the last event associated with an OTU.
        If no events can be found for this OTU, return None.
        """
        return self._index.get_latest_timestamp_by_otu_id(otu_id)

    def get_otu_last_updated(self, otu_id: uuid.UUID) -> datetime.datetime | None:
        """Get the timestamp of the last time this OTU was automatically updated.
        If this OTU has not been updated since this repo was initialized, return None.
        """
        return self._index.get_last_otu_update_timestamp(otu_id)

    def write_otu_update_history_entry(self, otu_id: uuid.UUID) -> int:
        """Add a new entry to the otu update history log and return the primary key
        of the entry.
        """
        if (
            update_id := self._index.add_otu_update_history_entry(
                otu_id,
                arrow.utcnow().naive,
            )
            is None
        ):
            raise SystemError(
                "OTU update history entry could not be retrieved after writing."
            )

        return update_id

    def get_event(self, event_id: int) -> Event | None:
        """Return event from event store."""
        try:
            return self._event_store.read_event(event_id)

        except FileNotFoundError:
            return None

    @staticmethod
    def _rehydrate_otu(events: Iterator[Event]) -> RepoOTU:
        """Rehydrate an OTU from an event iterator."""
        event = next(events)

        with warnings.catch_warnings(record=True) as warning_list:
            if isinstance(event, CreateOTU):
                otu = event.apply()

            else:
                raise TypeError(
                    f"The first event ({event}) for an OTU is not a CreateOTU " "event",
                )

            for event in events:
                if not isinstance(event, ApplicableEvent):
                    raise TypeError(
                        f"Event {event.id} {event.type} is not an applicable event."
                    )

                otu = event.apply(otu)

        for warning_msg in warning_list:
            logger.warning(
                warning_msg.message,
                otu_id=str(otu.id),
                warning_category=warning_msg.category.__name__,
            )

        otu.isolates.sort(
            key=lambda i: f"{i.name.type} {i.name.value}"
            if type(i.name) is IsolateName
            else "",
        )

        for isolate in otu.isolates:
            isolate.sequences.sort(key=lambda s: s.accession)

        return otu

    def _prune(self) -> None:
        """Prune events beyond the ID in the head file."""
        head_path = self.path / "head"

        if not head_path.exists():
            return

        with open(head_path) as f:
            head_id = int(f.read())

        self._event_store.prune(head_id)

    def _write_event(self, cls: type[Event], data: EventData, query: EventQuery):
        """Write an event to the repository."""
        if self._transaction is None:
            raise TransactionRequiredError

        event = self._event_store.write_event(
            cls(
                id=self.last_id + 1,
                data=data,
                query=query,
                timestamp=arrow.utcnow().naive,
            )
        )

        if hasattr(event.query, "otu_id"):
            self._index.add_event_id(
                event.id,
                event.query.otu_id,
                event.timestamp,
            )

            self._transaction.add_otu_id(event.query.otu_id)

        return event


@contextmanager
def locked_repo(path: Path) -> Generator[Repo, None, None]:
    """Yield a locked Repo."""
    repo = Repo(path)

    with repo.lock():
        yield repo

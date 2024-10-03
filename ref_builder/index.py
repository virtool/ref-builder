"""An index for improving repository performance."""

import binascii
import sqlite3
from collections.abc import Iterator
from dataclasses import dataclass
from pathlib import Path
from typing import Any
from uuid import UUID

import orjson

from ref_builder.models import OTUMinimal
from ref_builder.resources import RepoOTU


@dataclass
class EventIndexItem:
    """The event IDs associated with an OTU."""

    at_event: int
    """The latest event ID when this cache index was last updated"""

    event_ids: list[int]
    """A list of event IDs"""

    otu_id: UUID
    """The OTU ID this cache index is associated with."""


@dataclass
class Snapshot:
    """A snapshot of an OTU at a specific event."""

    at_event: int
    """The event at which the snapshot was taken."""

    otu: RepoOTU
    """The OTU that was snapshotted."""


class Index:
    """An index for rapid access to repository data.

    1. Get OTU IDs by legacy ID, name, or taxonomy ID.
    2. Save and load complete snapshots of OTUs.
    3. Get the event IDS associated with an OTU.

    """

    def __init__(self, path: Path) -> None:
        """Initialize the index."""
        self.path = path
        """The path to the index file."""

        self.path.parent.mkdir(exist_ok=True, parents=True)

        self.con = sqlite3.connect(path, isolation_level=None)
        self.con.execute("PRAGMA journal_mode = WAL")
        self.con.execute("PRAGMA synchronous = NORMAL")

        self.con.execute(
            """
            CREATE TABLE IF NOT EXISTS at_event (
                id INTEGER PRIMARY KEY,
                at_event INTEGER
            );
            """,
        )

        self.con.execute(
            """
            CREATE TABLE IF NOT EXISTS events (
                at_event INTEGER,
                event_id INTEGER PRIMARY KEY,
                otu_id TEXT
            );
            """,
        )

        self.con.execute(
            """
            CREATE TABLE IF NOT EXISTS otus (
                id TEXT PRIMARY KEY,
                acronym TEXT,
                at_event INTEGER,
                legacy_id TEXT,
                name TEXT,
                otu JSONB,
                taxid INTEGER
            );
            """,
        )

        self.con.execute(
            """
            CREATE TABLE IF NOT EXISTS sequences (
                id TEXT PRIMARY KEY,
                crc TEXT,
                otu_id TEXT,
                sequence TEXT
            );
            """,
        )

        for table, column in [
            ("events", "otu_id"),
            ("otus", "id"),
            ("otus", "legacy_id"),
            ("otus", "name"),
            ("otus", "taxid"),
            ("sequences", "otu_id"),
            ("sequences", "crc"),
        ]:
            self.con.execute(
                f"""
                CREATE INDEX IF NOT EXISTS {table}_{column} ON {table} ({column});
                """,
            )

        self.con.commit()

    @property
    def otu_ids(self) -> set[UUID]:
        """A list of all OTUs tracked in the index."""
        return {UUID(row[0]) for row in self.con.execute("SELECT id FROM otus")}

    def add_event_id(self, event_id: int, otu_id: UUID) -> None:
        """Add a new event ID and associated it with a given OTU ID.

        This allows fast lookup of event IDs for a given OTU ID.

        :param event_id: the event ID to add
        :param otu_id: the OTU ID to add the event ID to
        """
        self.con.execute(
            """
            INSERT INTO events (event_id, otu_id)
            VALUES (?, ?)
            """,
            (
                event_id,
                str(otu_id),
            ),
        )

    def get_event_ids_by_otu_id(self, otu_id: UUID) -> EventIndexItem | None:
        """Get a list of event IDs for an OTU."""
        res = self.con.execute(
            """
            SELECT event_id
            FROM events
            WHERE otu_id = ?
            ORDER BY event_id
            """,
            (str(otu_id),),
        )

        event_ids = [row[0] for row in res]

        if event_ids:
            return EventIndexItem(
                event_ids[-1],
                event_ids,
                otu_id,
            )

        return None

    def get_id_by_legacy_id(self, legacy_id: str) -> UUID | None:
        """Get an OTU ID by its legacy ID."""
        cursor = self.con.execute(
            "SELECT id FROM otus WHERE legacy_id = ?",
            (legacy_id,),
        )

        if result := cursor.fetchone():
            return UUID(result[0])

        return None

    def get_id_by_name(self, name: str) -> UUID | None:
        """Get an OTU ID by its name."""
        cursor = self.con.execute(
            "SELECT id FROM otus WHERE name = ?",
            (name,),
        )

        if result := cursor.fetchone():
            return UUID(result[0])

        return None

    def get_id_by_taxid(self, taxid: int) -> UUID | None:
        """Get an OTU ID by its taxonomy ID."""
        cursor = self.con.execute(
            "SELECT id FROM otus WHERE taxid = ?",
            (taxid,),
        )

        if result := cursor.fetchone():
            return UUID(result[0])

        return None

    def delete_otu(self, otu_id: UUID) -> None:
        """Remove an OTU from the index.

        This happens rarely, so we just rewrite the whole index.

        :param otu_id: the ID of the OTU to remove.
        """
        self.con.execute(
            "DELETE FROM events WHERE otu_id = ?",
            (str(otu_id),),
        )

        self.con.execute(
            "DELETE FROM otus WHERE id = ?",
            (str(otu_id),),
        )

        self.con.execute(
            "DELETE FROM sequences WHERE otu_id = ?",
            (str(otu_id),),
        )

        self.con.commit()

    def iter_minimal_otus(self) -> Iterator[OTUMinimal]:
        """Iterate over minimal representations of all OTUs in the index."""
        cursor = self.con.execute(
            "SELECT acronym, id, legacy_id, name, taxid FROM otus ORDER BY name",
        )

        for row in cursor:
            yield OTUMinimal(
                acronym=row[0],
                id=UUID(row[1]),
                legacy_id=row[2],
                name=row[3],
                taxid=row[4],
            )

    def load_snapshot(self, otu_id: UUID) -> Snapshot | None:
        """Load an OTU snapshot."""
        cursor = self.con.execute(
            "SELECT at_event, otu FROM otus WHERE id = ?",
            (str(otu_id),),
        )

        result = cursor.fetchone()

        if not result:
            return None

        at_event, otu_json = result

        otu = orjson.loads(otu_json)

        sequence_ids = [
            str(sequence["id"])
            for isolate in otu["isolates"]
            for sequence in isolate["sequences"]
        ]

        # Fetch all sequences in a single query
        placeholders = ",".join("?" for _ in sequence_ids)

        sequence_map = dict(
            self.con.execute(
                f"SELECT id, sequence FROM sequences WHERE id IN ({placeholders})",
                sequence_ids,
            ).fetchall(),
        )

        # Update the data structure and check for missing sequences
        missing_sequences = []

        for isolate in otu["isolates"]:
            for sequence in isolate["sequences"]:
                sequence_id = str(sequence["id"])

                if sequence_id in sequence_map:
                    sequence["sequence"] = sequence_map[sequence_id]
                else:
                    missing_sequences.append(sequence_id)

        # Raise an error if any sequences are missing
        if missing_sequences:
            raise ValueError(f"Sequences not found: {', '.join(missing_sequences)}")

        return Snapshot(
            at_event=at_event,
            otu=RepoOTU.model_validate(otu),
        )

    def upsert_otu(self, otu: RepoOTU, at_event: int) -> None:
        """Update the index based on a list of OTUs.

        OTUs that aren't already in the index will be added, and those that are will be
        updated.

        Only sequences that have changed will be updated.

        """
        sequence_ids = {
            str(seq.id) for isolate in otu.isolates for seq in isolate.sequences
        }

        otu_dict = otu.model_dump(mode="json")

        placeholders = ",".join("?" for _ in sequence_ids)

        # Delete any sequences that are no longer in the OTU.
        self.con.execute(
            f"""
            DELETE FROM sequences WHERE otu_id = ? AND NOT id IN ({ placeholders });
            """,
            (str(otu.id), *sequence_ids),
        )

        batch = []

        # Insert or update the sequences.
        for isolate in otu_dict["isolates"]:
            for sequence in isolate["sequences"]:
                crc = _calculate_crc32(sequence["sequence"])
                batch.append(
                    (
                        str(sequence["id"]),
                        crc,
                        str(otu.id),
                        sequence["sequence"],
                    ),
                )

        # Perform batch insert
        self.con.executemany(
            """
            INSERT OR REPLACE INTO sequences (id, crc, otu_id, sequence)
            VALUES (?, ?, ?, ?)
            """,
            batch,
        )

        # Update only if CRC has changed
        self.con.executemany(
            """
            UPDATE sequences
            SET crc = ?, otu_id = ?, sequence = ?
            WHERE id = ? AND crc != ?
            """,
            [(crc, otu_id, seq, seq_id, crc) for seq_id, crc, otu_id, seq in batch],
        )

        self.con.execute(
            """
            INSERT OR REPLACE INTO otus (id, acronym, at_event, legacy_id, name, otu, taxid)
            VALUES (?, ?,?, ?, ?, ?, ?)
            """,
            (
                str(otu.id),
                otu.acronym,
                at_event,
                otu.legacy_id,
                otu.name,
                orjson.dumps(otu_dict),
                otu.taxid,
            ),
        )

        self.con.commit()


def _calculate_crc32(sequence: str) -> str:
    """Calculate the CRC32 checksum of a sequence.

    :param sequence: the sequence to calculate the checksum of.
    :return: the CRC32 checksum as a hexadecimal string.

    """
    crc = binascii.crc32(sequence.encode("utf-8"))

    # Convert CRC as a hexadecimal string.
    return hex(crc & 0xFFFFFFFF)[2:].zfill(8)


def _default_json(obj: Any) -> list:
    """Cast sets to lists.

    This function is used as the default argument to `orjson.dumps` to handle sets.

    :param obj: the object to serialize
    :return: the serialized object
    """
    if isinstance(obj, set):
        return list(obj)

    raise TypeError(f"Object of type {obj.__class__.__name__} is not JSON serializable")

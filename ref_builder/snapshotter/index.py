"""An index for current OTU fields.

This index can be used to quickly find OTU IDs by various fields. It can be easily
updated by calling `update` with a list of OTUs.

"""

import binascii
import sqlite3
from dataclasses import dataclass
from pathlib import Path
from uuid import UUID

import orjson

from ref_builder.resources import RepoOTU


def default_json(obj):
    if isinstance(obj, set):
        return list(obj)

    raise TypeError(f"Object of type {obj.__class__.__name__} is not JSON serializable")


@dataclass
class Snapshot:
    at_event: int
    otu: RepoOTU


class Index:
    def __init__(self, path: Path):
        self.con = sqlite3.connect(path, isolation_level=None)

        self.con.execute("PRAGMA journal_mode = WAL")
        self.con.execute("PRAGMA synchronous = NORMAL")

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

    def delete(self, otu_id: UUID) -> None:
        """Remove an OTU from the index.

        This happens rarely, so we just rewrite the whole index.

        :param otu_id: the ID of the OTU to remove.
        """
        self.con.execute(
            "DELETE FROM otus WHERE id = ?",
            (str(otu_id),),
        )

        self.con.commit()

    def load(self, otu_id: UUID) -> Snapshot | None:
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
        sequence_map = dict(
            self.con.execute(
                "SELECT id, sequence FROM sequences WHERE id IN ({})".format(
                    ",".join("?" * len(sequence_ids)),
                ),
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

    @property
    def otu_ids(self) -> set[UUID]:
        """A list of OTU ids of snapshots."""
        return {UUID(row[0]) for row in self.con.execute("SELECT id FROM otus")}

    def update(self, otu: RepoOTU, at_event: int):
        """Update the index based on a list of OTUs.

        OTUs that aren't already in the index will be added, and those that are will be
        updated.

        Only sequences that have changed will be updated.

        """
        sequence_ids = {
            str(seq.id) for isolate in otu.isolates for seq in isolate.sequences
        }

        otu_dict = otu.model_dump()

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
                crc = calculate_crc32(sequence["sequence"])
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
                orjson.dumps(otu_dict, default=default_json),
                otu.taxid,
            ),
        )

        self.con.commit()


def calculate_crc32(sequence: str) -> str:
    """Calculate the CRC32 checksum of a sequence.

    :param sequence: the sequence to calculate the checksum of.
    :return: the CRC32 checksum as a hexadecimal string.

    """
    crc = binascii.crc32(sequence.encode("utf-8"))

    # Convert CRC as a hexadecimal string.
    return hex(crc & 0xFFFFFFFF)[2:].zfill(8)

"""An index for current OTU fields.

This index can be used to quickly find OTU IDs by various fields. It can be easily
updated by calling `update` with a list of OTUs.

"""

import sqlite3
from dataclasses import dataclass
from pathlib import Path
from uuid import UUID

import orjson

from ref_builder.resources import RepoOTU


@dataclass
class Snapshot:
    at_index: int
    otu: RepoOTU


class Index:
    def __init__(self, path: Path):
        self.con = sqlite3.connect(path, isolation_level=None)

        self.con.execute("PRAGMA journal_mode = WAL")

        self.con.execute(
            """
            CREATE TABLE IF NOT EXISTS otus (
                id TEXT PRIMARY KEY,
                acronym TEXT,
                at_index INTEGER,
                legacy_id TEXT,
                name TEXT,
                taxid INTEGER,
                snapshot BLOB
            );
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

    def load(self, otu_id: UUID):
        """Load an OTU snapshot."""
        cursor = self.con.execute(
            "SELECT at_index, snapshot FROM otus WHERE id = ?",
            (str(otu_id),),
        )

        if result := cursor.fetchone():
            at_index, snapshot = result

            snapshot = orjson.loads(snapshot)

            return Snapshot(
                at_index=at_index,
                otu=RepoOTU.from_dict(snapshot),
            )

        return None

    @property
    def otu_ids(self) -> set[UUID]:
        """A list of OTU ids of snapshots."""
        return {UUID(row[0]) for row in self.con.execute("SELECT id FROM otus")}

    def update(self, otu: RepoOTU, at_index: int):
        """Update the index based on a list of OTUs.

        OTUs that aren't already in the index will be added, and those that are will be
        updated.

        """
        self.con.execute(
            """
            INSERT OR REPLACE INTO otus (id, acronym, at_index, legacy_id, name, taxid, snapshot)
            VALUES (?, ?,?, ?, ?, ?, ?)
            """,
            (
                str(otu.id),
                otu.acronym,
                at_index,
                otu.legacy_id,
                otu.name,
                otu.taxid,
                orjson.dumps(otu.dict()),
            ),
        )

        self.con.commit()


def calculate_crc32(sequence: str) -> str:
    """Calculate the CRC32 checksum of a sequence.

    :param sequence: the sequence to calculate the checksum of.
    :return: the CRC32 checksum as a hexadecimal string.

    """
    crc = zlib.crc32(sequence.encode("utf-8"))

    # Convert CRC as a hexadecimal string.
    return hex(crc & 0xFFFFFFFF)[2:].zfill(8)

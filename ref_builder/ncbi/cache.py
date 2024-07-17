import json
import shutil
from pathlib import Path

from structlog import get_logger

base_logger = get_logger()


class NCBICache:
    """Manages caching functionality for NCBI data"""

    def __init__(self, path: Path) -> None:
        """Initialize the cache with a path to store cached data.

        :param path: A directory that will store cached data
        """
        self.path = path

        self._genbank_path = self.path / "genbank"
        self._taxonomy_path = self.path / "taxonomy"

        self.path.mkdir(exist_ok=True)
        self._genbank_path.mkdir(exist_ok=True)
        self._taxonomy_path.mkdir(exist_ok=True)

    def clear(self) -> None:
        """Clear and reset the cache."""
        shutil.rmtree(self.path)
        self.path.mkdir()

        self._genbank_path.mkdir()
        self._taxonomy_path.mkdir()

    def cache_genbank_record(self, data: dict, accession: str) -> None:
        """Add a Genbank record from NCBI Nucleotide to the cache.

        :param data: A data from a Genbank record corresponding
        :param accession: The NCBI accession of the record
        """
        cached_record_path = self._get_genbank_path(f"{accession}")

        with open(cached_record_path, "w") as f:
            json.dump(data, f)

            if not cached_record_path.exists():
                raise FileNotFoundError

    def load_genbank_record(self, accession: str) -> dict | None:
        """Retrieve a NCBI Nucleotide Genbank record from the cache.

        Returns ``None`` if the record is not found in the cache.

        :param accession: The NCBI accession of the record
        :return: Deserialized Genbank data if file is found in cache, else None
        """
        record_path = self._get_genbank_path(accession)

        try:
            with open(record_path) as f:
                return json.load(f)

        except FileNotFoundError:
            return None

    def cache_taxonomy_record(self, data: dict, taxid):
        """Add a NCBI Taxonomy record to the cache

        :param data: NCBI Taxonomy record data
        :param taxid: A NCBI Taxonomy id
        """
        cached_taxonomy_path = self._get_taxonomy_path(taxid)

        with open(cached_taxonomy_path, "w") as f:
            json.dump(data, f)

        if not cached_taxonomy_path.exists():
            raise FileNotFoundError

    def load_taxonomy(self, taxid: int) -> dict | None:
        """Load data from a cached record fetch

        :param taxid: A NCBI Taxonomy id
        :return: Deserialized Taxonomy data if file is found in cache, else None
        """
        try:
            with open(self._get_taxonomy_path(taxid)) as f:
                return json.load(f)
        except FileNotFoundError:
            return None

    def _get_genbank_path(self, accession: str) -> Path:
        """Returns a standardized path for a set of cached NCBI Nucleotide records

        :param accession: The NCBI accession of a Genbank record
        :return: A properly-formatted path to a cached record
        """
        return self._genbank_path / f"{accession}.json"

    def _get_taxonomy_path(self, taxid: int) -> Path:
        """Returns a standardized path for a cached NCBI Taxonomy record

        :param taxid: A NCBI Taxonomy id
        :return: A properly-formatted path to a cached record
        """
        return self._taxonomy_path / f"{taxid}.json"

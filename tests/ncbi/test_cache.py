import pytest
from syrupy import SnapshotAssertion

from ref_builder.ncbi.cache import NCBICache


def test_clear(scratch_ncbi_cache: NCBICache):
    """Test that the cache is cleared correctly.

    Start with a populated cache and make sure data is in it. Then clear it and make
    sure the cache directories are empty.
    """
    genbank_path = scratch_ncbi_cache.path / "genbank"
    taxonomy_path = scratch_ncbi_cache.path / "taxonomy"

    assert len(list(genbank_path.glob("*.json"))) == 74
    assert len(list(taxonomy_path.glob("*.json"))) == 26

    scratch_ncbi_cache.clear()

    assert len(list(genbank_path.glob("*.json"))) == 0
    assert len(list(taxonomy_path.glob("*.json"))) == 0


class TestGenbank:
    """Test the caching and loading of Genbank records."""

    @pytest.mark.parametrize(
        "accession",
        ["AB017504", "NC_003355"],
    )
    def test_cache(
        self,
        accession: str,
        scratch_ncbi_cache: NCBICache,
        snapshot: SnapshotAssertion,
    ):
        """Test that a record is cached, cleared, and loaded correctly."""
        record = scratch_ncbi_cache.load_genbank_record(accession)

        scratch_ncbi_cache.clear()

        assert scratch_ncbi_cache.load_genbank_record(accession) is None

        scratch_ncbi_cache.cache_genbank_record(record, accession, 1)

        assert scratch_ncbi_cache.load_genbank_record(accession) == record == snapshot

    def test_load_not_found(self, scratch_ncbi_cache: NCBICache):
        """Test that None is returned if a record is not found."""
        assert scratch_ncbi_cache.load_genbank_record("not_found") is None


class TestTaxonomy:
    """Test the caching and loading of taxonomy records."""

    def test_load_and_cache(
        self,
        scratch_ncbi_cache: NCBICache,
        snapshot: SnapshotAssertion,
    ):
        """Test that a taxonomy record can be loaded and cached correctly."""
        record = scratch_ncbi_cache.load_taxonomy(1198450)

        scratch_ncbi_cache.clear()

        assert scratch_ncbi_cache.load_taxonomy(1198450) is None

        scratch_ncbi_cache.cache_taxonomy_record(record, 1198450)

        assert scratch_ncbi_cache.load_taxonomy(1198450) == record == snapshot

    def test_load_not_found(
        self,
        scratch_ncbi_cache: NCBICache,
    ):
        """Test that None is returned if a record is not found."""
        assert scratch_ncbi_cache.load_taxonomy(101010) is None

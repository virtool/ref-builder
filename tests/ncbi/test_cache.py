import json
from pathlib import Path

import pytest
from syrupy import SnapshotAssertion

from ref_builder.ncbi.cache import NCBICache


@pytest.fixture()
def _empty_cache_path(tmp_path) -> Path:
    """A path to an empty cache."""
    return tmp_path / "empty_cache"


def _get_test_genbank(accession: str, cache_example_path: Path) -> dict:
    with open(cache_example_path / "ncbi" / "genbank" / f"{accession}.json") as f:
        return json.load(f)


def _get_test_taxonomy(taxon_id: int, cache_example_path: Path) -> dict:
    with open(cache_example_path / "ncbi" / "taxonomy" / f"{taxon_id}.json") as f:
        return json.load(f)


def test_init(_empty_cache_path: Path):
    """Test that the cache is initialized correctly."""
    assert not _empty_cache_path.exists()

    cache = NCBICache(path=_empty_cache_path)

    assert (cache.path / "genbank").is_dir()
    assert (cache.path / "taxonomy").is_dir()


def test_clear(scratch_user_cache_path: Path):
    """Test that the cache is cleared correctly.

    Start with a populated cache and make sure data is in it. Then clear it and make
    sure the cache directories are empty.
    """
    cache = NCBICache(path=scratch_user_cache_path)
    genbank_path = cache.path / "genbank"
    taxonomy_path = cache.path / "taxonomy"

    assert len(list(genbank_path.glob("*.json"))) == 74
    assert len(list(taxonomy_path.glob("*.json"))) == 26

    cache.clear()

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
        scratch_user_cache_path: Path,
        snapshot: SnapshotAssertion,
    ):
        """Test that a record is cached correctly."""
        cache = NCBICache(scratch_user_cache_path)

        record = _get_test_genbank(accession, scratch_user_cache_path)

        cache.cache_genbank_record(data=record, accession=accession)

        with open(
            scratch_user_cache_path / "ncbi" / "genbank" / f"{accession}.json",
        ) as f:
            assert json.load(f) == snapshot

    @pytest.mark.parametrize(
        "accession",
        ["AB017504", "MH200607", "MT240513", "NC_003355"],
    )
    def test_load(
        self,
        accession: str,
        scratch_user_cache_path: Path,
        snapshot: SnapshotAssertion,
    ):
        """Test that a record is loaded from the cache."""
        cache = NCBICache(scratch_user_cache_path)
        assert cache.load_genbank_record(accession) == snapshot

    def test_load_not_found(self, scratch_user_cache_path: Path):
        """Test that None is returned if a record is not found."""
        cache = NCBICache(scratch_user_cache_path)
        assert cache.load_genbank_record("not_found") is None


class TestTaxonomy:
    """Test the caching and loading of taxonomy records."""

    @pytest.mark.parametrize("taxid", [270478, 438782, 1198450])
    def test_load(
        self,
        taxid: int,
        scratch_user_cache_path: Path,
        snapshot: SnapshotAssertion,
    ):
        """Test that a taxonomy record is loaded correctly."""
        scratch_ncbi_cache = NCBICache(scratch_user_cache_path)

        taxonomy = scratch_ncbi_cache.load_taxonomy(taxid)

        assert taxonomy == snapshot(name=f"{taxid}.json")

    def test_load_not_found(
        self,
        scratch_user_cache_path: Path,
    ):
        """Test that None is returned if a record is not found."""
        scratch_ncbi_cache = NCBICache(scratch_user_cache_path)
        assert scratch_ncbi_cache.load_taxonomy(101010) is None

    def test_cache(
        self,
        _empty_cache_path: Path,
        scratch_user_cache_path: Path,
    ):
        """Test that a taxonomy record is cached correctly."""
        taxonomy = _get_test_taxonomy(1198450, scratch_user_cache_path)

        fresh_cache = NCBICache(_empty_cache_path)

        fresh_cache.cache_taxonomy_record(taxonomy, 1198450)

        assert (fresh_cache.path / "taxonomy" / f"{1198450}.json").exists()

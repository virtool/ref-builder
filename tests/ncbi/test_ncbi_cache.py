import json
from pathlib import Path

import pytest
from syrupy import SnapshotAssertion

from ref_builder.ncbi.cache import NCBICache


@pytest.fixture()
def empty_cache_path(tmp_path):
    return tmp_path / "test_cache"


def get_test_record(accession: str, cache_example_path) -> dict:
    with open(cache_example_path / "genbank" / f"{accession}.json") as f:
        return json.load(f)


def get_test_taxonomy(taxon_id: int, cache_example_path) -> dict:
    with open(cache_example_path / "taxonomy" / f"{taxon_id}.json") as f:
        return json.load(f)


def test_cache_init(empty_cache_path):
    assert not empty_cache_path.exists()

    cache = NCBICache(path=empty_cache_path)

    assert cache._genbank_path.is_dir()

    assert cache._taxonomy_path.is_dir()

    assert cache._genbank_path.exists()

    assert cache._taxonomy_path.exists()


def test_cache_clear(scratch_ncbi_cache_path):
    cache = NCBICache(path=scratch_ncbi_cache_path)

    assert list(cache._genbank_path.glob("*.json")) != []
    assert list(cache._taxonomy_path.glob("*.json")) != []

    cache.clear()

    assert list(cache._genbank_path.glob("*.json")) == []
    assert list(cache._taxonomy_path.glob("*.json")) == []


@pytest.mark.parametrize(
    "accessions",
    (
        ["AB017504", "MH200607", "MK431779", "NC_003355"],
        ["NC_036587", "MT240513", "MT240490"],
    ),
)
class TestCacheGenbankOperations:
    def test_cache_genbank_load_record_batch(
        self,
        accessions,
        scratch_ncbi_cache_path,
        snapshot: SnapshotAssertion,
    ):
        scratch_ncbi_cache = NCBICache(scratch_ncbi_cache_path)

        for accession in accessions:
            record = scratch_ncbi_cache.load_genbank_record(accession)

            assert record == snapshot(name=f"{accession}.json")

    def test_cache_genbank_cache_records(
        self,
        accessions,
        scratch_ncbi_cache_path: Path,
        empty_cache_path,
    ):
        assert not empty_cache_path.exists()

        cache = NCBICache(empty_cache_path)

        for accession in accessions:
            record = get_test_record(accession, scratch_ncbi_cache_path)

            cache.cache_genbank_record(data=record, accession=accession)

            assert (cache._genbank_path / f"{accession}.json").exists()


@pytest.mark.parametrize("fake_accession", ["afjshd", "23222", "wheelhouse"])
def test_cache_genbank_load_fail(fake_accession, scratch_ncbi_cache_path):
    scratch_ncbi_cache = NCBICache(scratch_ncbi_cache_path)

    assert scratch_ncbi_cache.load_genbank_record(fake_accession) is None


@pytest.mark.parametrize("taxid", (270478, 438782, 1198450))
class TestCacheTaxonomyOperations:
    def test_cache_taxonomy_load(
        self,
        taxid,
        scratch_ncbi_cache_path,
        snapshot: SnapshotAssertion,
    ):
        scratch_ncbi_cache = NCBICache(scratch_ncbi_cache_path)

        taxonomy = scratch_ncbi_cache.load_taxonomy(taxid)

        assert taxonomy == snapshot(name=f"{taxid}.json")

    def test_cache_taxonomy_cache(
        self,
        taxid: int,
        scratch_ncbi_cache_path: Path,
        empty_cache_path,
    ):
        taxonomy = get_test_taxonomy(taxid, scratch_ncbi_cache_path)

        fresh_cache = NCBICache(empty_cache_path)

        fresh_cache.cache_taxonomy_record(taxonomy, taxid)

        assert (fresh_cache._taxonomy_path / f"{taxid}.json").exists()

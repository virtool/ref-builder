import shutil
from pathlib import Path

import pytest
from structlog import get_logger
from syrupy import SnapshotAssertion

from ref_builder.ncbi.client import NCBIClient
from ref_builder.ncbi.models import NCBIGenbank, NCBISource, NCBITaxonomy

test_logger = get_logger()

SERVER_ERRORS = ["HTTPError", "IncompleteRead"]


@pytest.fixture()
def empty_client(tmp_path: Path) -> NCBIClient:
    path = tmp_path / "dummy_cache"
    path.mkdir()

    return NCBIClient(path, ignore_cache=True)


class TestClientFetchGenbank:
    @pytest.mark.ncbi()
    @pytest.mark.parametrize(
        "accessions",
        [
            ["AB017504", "MH200607", "MK431779", "NC_003355"],
            ["NC_036587", "MT240513", "MT240490"],
        ],
    )
    def test_fetch_genbank_records_from_ncbi(
        self,
        accessions: list[str],
        empty_client: NCBIClient,
        snapshot: SnapshotAssertion,
    ):
        clean_records = empty_client.fetch_genbank_records(accessions=accessions)

        assert clean_records

        for record in clean_records:
            assert type(record) is NCBIGenbank

            assert type(record.source) is NCBISource

        assert {
            path.name for path in empty_client.cache._genbank_path.glob("*.json")
        } == snapshot

    @pytest.mark.parametrize(
        "accessions",
        [
            ["AB017504", "MH200607", "MK431779", "NC_003355"],
            ["NC_036587", "MT240513", "MT240490"],
        ],
    )
    def test_fetch_genbank_records_from_cache(
        self,
        accessions: list[str],
        scratch_user_cache_path,
        snapshot: SnapshotAssertion,
    ):
        client = NCBIClient(scratch_user_cache_path, ignore_cache=False)

        clean_records = client.fetch_genbank_records(accessions=accessions)

        assert clean_records

        for record in clean_records:
            assert record == snapshot(name=f"{record.accession}_validated.json")

    @pytest.mark.ncbi()
    @pytest.mark.parametrize(
        "accessions",
        [["AB017503", "AB017504", "MH200607", "MK431779", "NC_003355"]],
    )
    def test_fetch_partially_cached_genbank_records(
        self,
        accessions: list[str],
        empty_client,
        files_path,
    ):
        ncbi_cache_files = files_path / "cache_test"
        for accession in ("AB017504", "MH200607"):
            shutil.copy(
                ncbi_cache_files / f"genbank/{accession}.json",
                empty_client.cache._genbank_path / f"{accession}.json",
            )

        assert list(empty_client.cache._genbank_path.glob("*.json")) != []

        clean_records = empty_client.fetch_genbank_records(accessions=accessions)

        assert clean_records

        cache_contents = [
            path.stem for path in empty_client.cache._genbank_path.glob("*.json")
        ]

        for accession in accessions:
            assert accession in cache_contents

    @pytest.mark.ncbi()
    def test_fetch_accessions_fail(self, scratch_ncbi_client: NCBIClient):
        records = scratch_ncbi_client.fetch_genbank_records(["friday", "paella", "111"])

        assert not records


class TestClientFetchRawGenbank:
    @pytest.mark.ncbi()
    @pytest.mark.parametrize(
        "accessions",
        [
            ["AB017504", "MH200607", "MK431779", "NC_003355"],
            ["NC_036587", "MT240513", "MT240490"],
        ],
    )
    def test_fetch_raw_via_accessions(self, accessions: list[str]):
        records = NCBIClient.fetch_unvalidated_genbank_records(accessions)

        for record in records:
            assert record.get("GBSeq_locus")
            assert record.get("GBSeq_sequence")

    @pytest.mark.ncbi()
    @pytest.mark.parametrize(
        "accessions",
        [
            ["paella", "MH200607", "MK431779", "NC_003355"],
            ["friday", "MT240513", "MT240490"],
        ],
    )
    def test_fetch_raw_via_accessions_partial(self, accessions: list[str]):
        records = NCBIClient.fetch_unvalidated_genbank_records(accessions)

        assert len(records) == len(accessions) - 1

        for record in records:
            assert record.get("GBSeq_locus", None)
            assert record.get("GBSeq_sequence", None)

    @pytest.mark.ncbi()
    def test_fetch_raw_via_accessions_fail(self):
        records = NCBIClient.fetch_unvalidated_genbank_records(
            ["friday", "paella", "111"],
        )

        assert not records


@pytest.mark.ncbi()
@pytest.mark.parametrize("taxid", [345184, 1198450, 1016856])
def test_fetch_records_by_taxid(taxid: int, empty_client, snapshot: SnapshotAssertion):
    assert not list(empty_client.cache._genbank_path.glob("*.json"))

    records = empty_client.link_from_taxid_and_fetch(taxid)

    assert records

    for record in records:
        assert type(record) is NCBIGenbank

    assert {
        path.name for path in empty_client.cache._genbank_path.glob("*.json")
    } == snapshot


class TestClientFetchTaxonomy:
    @pytest.mark.ncbi()
    @pytest.mark.parametrize("taxid", [438782, 1198450, 1016856])
    def test_fetch_taxonomy_from_ncbi(
        self,
        taxid: int,
        empty_client,
        snapshot: SnapshotAssertion,
    ):
        taxonomy = empty_client.fetch_taxonomy_record(taxid)

        assert type(taxonomy) is NCBITaxonomy

        assert {
            path.name for path in empty_client.cache._taxonomy_path.glob("*.json")
        } == snapshot

    @pytest.mark.ncbi()
    @pytest.mark.parametrize("taxid", [1000000000000, 99999999])
    def test_fetch_taxonomy_from_ncbi_fail(self, taxid: int, empty_client):
        taxonomy = empty_client.fetch_taxonomy_record(taxid)

        assert taxonomy is None

    @pytest.mark.parametrize("taxid", [438782, 1198450])
    def test_fetch_taxonomy_from_cache(
        self,
        taxid: int,
        scratch_ncbi_client,
        snapshot: SnapshotAssertion,
    ):
        assert scratch_ncbi_client.cache.load_taxonomy(taxid)

        taxonomy = scratch_ncbi_client.fetch_taxonomy_record(taxid)

        assert type(taxonomy) is NCBITaxonomy

        assert taxonomy == snapshot(name=f"{taxonomy}_validated.json")


@pytest.mark.ncbi()
@pytest.mark.parametrize(
    "name,taxid",
    [
        ("Apple rubbery wood virus 1", 2164102),
        ("Rhynchosia golden mosaic virus", 117198),
    ],
)
def test_fetch_taxonomy_by_name(name: str, taxid: int):
    assert NCBIClient.fetch_taxonomy_id_by_name(name) == taxid


# @pytest.mark.skip("ESpell is encountering issues.")
@pytest.mark.ncbi()
@pytest.mark.parametrize(
    "misspelled,expected",
    [
        (
            "Hynchosia yellow mosaic India virus",
            "rhynchosia yellow mosaic india virus",
        ),
        (
            "Angelica bush stunt virus",
            "angelica bushy stunt virus",
        ),
    ],
)
def test_fetch_spelling(misspelled: str, expected: str):
    taxon_name = NCBIClient.fetch_spelling(name=misspelled)

    assert taxon_name == expected

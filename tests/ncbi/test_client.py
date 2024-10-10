import pytest
from syrupy import SnapshotAssertion

from ref_builder.ncbi.client import NCBIClient


class TestFetchGenbank:
    @pytest.mark.ncbi()
    def test_fetch_genbank_records_from_ncbi(
        self,
        snapshot: SnapshotAssertion,
        uncached_ncbi_client: NCBIClient,
    ):
        accessions = ["NC_036587", "MT240513", "AB017504"]

        assert all(
            uncached_ncbi_client.cache.load_genbank_record(accession) is None
            for accession in accessions
        )

        assert uncached_ncbi_client.fetch_genbank_records(accessions) == snapshot

        assert all(
            uncached_ncbi_client.cache.load_genbank_record(accession)
            for accession in accessions
        )

    def test_fetch_genbank_records_from_cache(
        self,
        snapshot: SnapshotAssertion,
        scratch_ncbi_client: NCBIClient,
    ):
        assert (
            scratch_ncbi_client.fetch_genbank_records(
                ["NC_036587", "MT240513", "AB017504"],
            )
            == snapshot
        )

    @pytest.mark.ncbi()
    def test_fetch_partially_cached_genbank_records(
        self,
        snapshot: SnapshotAssertion,
        uncached_ncbi_client: NCBIClient,
    ):
        """Test that the client can fetch records when some are cached and some aren't.

        Start by fetching one accession, then ensure it was cached. Then, fetch is again
        along with two other accessions.
        """
        assert uncached_ncbi_client.cache.load_genbank_record("NC_036587") is None
        assert uncached_ncbi_client.fetch_genbank_records(["NC_036587"])
        assert uncached_ncbi_client.cache.load_genbank_record("NC_036587")

        assert (
            uncached_ncbi_client.fetch_genbank_records(
                ["NC_036587", "MT240513", "AB017504"],
            )
            == snapshot
        )

    @pytest.mark.ncbi()
    def test_fetch_non_existent_accession(self, scratch_ncbi_client: NCBIClient):
        """Test that the client returns an empty list when the fetched accession does
        not exist.
        """
        assert scratch_ncbi_client.fetch_genbank_records(["paella"]) == []


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
def test_fetch_records_by_taxid(
    uncached_ncbi_client: NCBIClient,
    snapshot: SnapshotAssertion,
):
    """Test that the client can fetch records by taxid and that they are cached."""
    records = uncached_ncbi_client.link_from_taxid_and_fetch(1198450)

    assert records == snapshot

    assert all(
        uncached_ncbi_client.cache.load_genbank_record(record.accession)
        for record in records
    )


class TestFetchTaxonomy:
    @pytest.mark.ncbi()
    def test_ok(
        self,
        snapshot: SnapshotAssertion,
        uncached_ncbi_client: NCBIClient,
    ):
        """Test that the client can fetch a taxonomy record and that it is cached."""
        # Make sure the taxid is not cached.
        assert uncached_ncbi_client.cache.load_taxonomy(1198450) is None

        # Fetch the taxonomy record and check its contents.
        assert uncached_ncbi_client.fetch_taxonomy_record(1198450) == snapshot

        # Make sure the taxid is now cached.
        assert uncached_ncbi_client.cache.load_taxonomy(1198450)

    @pytest.mark.ncbi()
    def test_not_found(self, uncached_ncbi_client: NCBIClient):
        """Test that the client returns None when the taxid does not exist."""
        assert uncached_ncbi_client.fetch_taxonomy_record(99999999) is None


@pytest.mark.ncbi()
def test_fetch_taxonomy_by_name():
    """Test that the client can fetch a taxid by name."""
    assert (
        NCBIClient.fetch_taxonomy_id_by_name("Rhynchosia golden mosaic virus") == 117198
    )


@pytest.mark.ncbi()
def test_fetch_spelling():
    """Test that the client can fetch the correct spelling of a virus name."""
    assert (
        NCBIClient.fetch_spelling(name="Angelica bush stunt virus")
        == "angelica bushy stunt virus"
    )
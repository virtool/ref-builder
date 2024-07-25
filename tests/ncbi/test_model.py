import pytest
from pydantic import ValidationError
from syrupy import SnapshotAssertion

from ref_builder.ncbi.cache import NCBICache
from ref_builder.ncbi.models import NCBIGenbank, NCBILineage, NCBIRank, NCBITaxonomy


class TestParseGenbank:
    @pytest.mark.parametrize(
        "accession",
        ["AB017504", "MH200607", "NC_036587", "MT240513", "NC_015504"],
    )
    def test_ok(
        self,
        accession: str,
        scratch_ncbi_cache: NCBICache,
        snapshot: SnapshotAssertion,
    ):
        """Test that multiple valid records can be parsed successfully."""
        record = scratch_ncbi_cache.load_genbank_record(accession)
        assert NCBIGenbank(**record).model_dump() == snapshot

    def test_source(
        self,
        scratch_ncbi_cache: NCBICache,
        snapshot: SnapshotAssertion,
    ):
        """Test that the source table is correctly extracted from the feature table."""
        record = scratch_ncbi_cache.load_genbank_record("AB017504")
        assert NCBIGenbank(**record).source.model_dump() == snapshot

    def test_taxid(
        self,
        scratch_ncbi_cache,
    ):
        """Test that the taxid is correctly extracted from the source field."""
        record = scratch_ncbi_cache.load_genbank_record("AB017504")
        assert NCBIGenbank(**record).source.taxid == 1169032

    def test_sequence_validation_fail(self, scratch_ncbi_cache: NCBICache):
        """Test that validation fails when the sequence contains invalid characters."""
        record = scratch_ncbi_cache.load_genbank_record("AB017504")

        assert NCBIGenbank(**record)

        record["GBSeq_sequence"] = "naa"

        try:
            NCBIGenbank(**record)
        except ValidationError as exc:
            for error in exc.errors():
                assert "GBSeq_sequence" in error["loc"]


class TestParseTaxonomy:
    @pytest.mark.parametrize("taxid", [270478, 438782, 1198450, 1077859])
    def test_ok(
        self,
        taxid: int,
        scratch_ncbi_cache: NCBICache,
        snapshot: SnapshotAssertion,
    ):
        """Test that multiple valid records can be parsed successfully."""
        record = scratch_ncbi_cache.load_taxonomy(taxid)

        taxonomy = NCBITaxonomy(**record)

        assert type(taxonomy) is NCBITaxonomy
        assert taxonomy.model_dump() == snapshot

    def test_with_rank(
        self,
        scratch_ncbi_cache: NCBICache,
        snapshot: SnapshotAssertion,
    ):
        """Test that the rank field can be explicitly set using a kwarg in the case that
        it is not present in the record.
        """
        record = scratch_ncbi_cache.load_taxonomy(1016856)

        with pytest.raises(ValidationError):
            assert NCBITaxonomy(**record)

        taxonomy = NCBITaxonomy(rank="isolate", **record)

        assert taxonomy.rank == NCBIRank.ISOLATE
        assert taxonomy.model_dump() == snapshot


def test_create_lineage_item_alias():
    """Test that the NCBILineage model can be created with or without aliases."""
    lineage_data_with_aliases = {
        "TaxId": "2732397",
        "ScientificName": "Pararnavirae",
        "Rank": "kingdom",
    }

    assert NCBILineage(**lineage_data_with_aliases) == NCBILineage(
        id=2732397,
        name="Pararnavirae",
        rank="kingdom",
    )

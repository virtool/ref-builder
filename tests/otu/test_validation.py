from ref_builder.repo import Repo
from ref_builder.otu.validate import check_otu_is_valid


class TestValidateOTU:
    """Test OTU validation."""

    def test_ok(self, scratch_repo: Repo):
        """Validate scratch_repo contents"""
        for otu_metadata in scratch_repo.iter_minimal_otus():
            assert check_otu_is_valid(scratch_repo.get_otu_by_taxid(otu_metadata.taxid))

    def test_no_isolates(self, scratch_repo: Repo):
        """Test that validating an OTU with no isolates returns False."""
        otu_clean = scratch_repo.get_otu_by_taxid(223262)

        assert check_otu_is_valid(otu_clean)

        otu_invalid = otu_clean.model_copy()
        otu_invalid.delete_isolate(otu_clean.representative_isolate)

        assert not otu_invalid.isolates

        assert not check_otu_is_valid(otu_invalid)

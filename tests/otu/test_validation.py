from ref_builder.repo import Repo
from ref_builder.otu.validate import validate_otu, check_repo_for_invalid_otus


class TestValidateOTU:
    """Test OTU validation."""
    
    def test_ok(self, scratch_repo: Repo):
        """Validate scratch_repo contents"""
        for otu_metadata in scratch_repo.iter_minimal_otus():
            assert validate_otu(scratch_repo.get_otu_by_taxid(otu_metadata.taxid))

    def test_no_isolates(self, scratch_repo: Repo):
        """Test that validating an OTU with no isolates returns False."""
        otu_clean = scratch_repo.get_otu_by_taxid(223262)

        assert validate_otu(otu_clean)

        otu_invalid = otu_clean.model_copy()
        otu_invalid.delete_isolate(otu_clean.representative_isolate)

        assert not otu_invalid.isolates

        assert not validate_otu(otu_invalid)


class TestValidateRepo:
    """Test repo validation."""

    def test_ok(self, scratch_repo: Repo):
        """Check that there are no invalid OTUs in scratch_repo"""
        assert not check_repo_for_invalid_otus(scratch_repo)

from ref_builder.repo import Repo
from ref_builder.otu.validate import validate_otu


class TestValidateOTU:
    """Test OTU validation."""
    
    def test_ok(self, scratch_repo: Repo):
        """Validate scratch_repo contents"""
        for otu_metadata in scratch_repo.iter_minimal_otus():
            assert validate_otu(scratch_repo.get_otu_by_taxid(otu_metadata.taxid))

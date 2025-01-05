import pytest

from ref_builder.otu.utils import assign_records_to_segments
from ref_builder.repo import Repo
from tests.fixtures.factories import NCBIGenbankFactory


@pytest.mark.parametrize("taxid", [223262, 3158377])
class TestAssignRecordsToSegments:
    """Test whether a list of records accurately maps to the segment list of a
    multipartite OTU's plan.
    """

    def test_ok(self, scratch_repo: Repo, taxid: int):
        """Test that a list of records generated to match the OTU plan
        can pass assign_records_to_segments().
        """
        otu = scratch_repo.get_otu_by_taxid(taxid)

        mock_records = NCBIGenbankFactory.build_from_metadata(
            plan=otu.plan,
            moltype=otu.molecule.type,
            name=otu.name,
            taxid=otu.taxid,
        )

        assigned_records = assign_records_to_segments(mock_records, otu.plan)

        assert assigned_records.keys() == {
            segment.id for segment in otu.plan.required_segments
        }

    def test_fail(self, scratch_repo: Repo, taxid: int):
        """Test that a randomly generated list of records raises an appropriate
        ValueError.
        """
        otu = scratch_repo.get_otu_by_taxid(taxid)

        mock_records = NCBIGenbankFactory.batch(len(otu.plan.required_segments))

        with pytest.raises(ValueError, match="Missing one or more required segments"):
            assign_records_to_segments(mock_records, otu.plan)

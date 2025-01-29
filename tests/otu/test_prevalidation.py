import pytest
from syrupy import SnapshotAssertion
from syrupy.filters import props

from ref_builder.ncbi.models import NCBISourceMolType
from ref_builder.otu.utils import assign_records_to_segments
from ref_builder.repo import Repo
from tests.fixtures.factories import NCBIGenbankFactory, NCBISourceFactory


@pytest.mark.parametrize("taxid", [223262, 3158377])
class TestAssignRecordsToSegments:
    """Test whether a list of records can be assigned to plan segments."""

    def test_ok(
        self,
        ncbi_genbank_factory: NCBIGenbankFactory,
        ncbi_source_factory: NCBISourceFactory,
        scratch_repo: Repo,
        snapshot: SnapshotAssertion,
        taxid: int,
    ):
        """Test that a list of records generated to match the OTU plan
        can pass assign_records_to_segments().
        """
        otu = scratch_repo.get_otu_by_taxid(taxid)

        records = [
            ncbi_genbank_factory.build(
                source=ncbi_source_factory.build(
                    mol_type=NCBISourceMolType.from_molecule(otu.molecule),
                    segment=str(segment.name),
                    organism=otu.name,
                    taxid=otu.taxid,
                )
            )
            for segment in otu.plan.required_segments
        ]

        assigned_records = assign_records_to_segments(records, otu.plan)

        assert sorted(assigned_records.values(), key=lambda r: r.accession) == snapshot(
            exclude=props("id")
        )

    def test_names_not_in_plan(
        self,
        ncbi_genbank_factory: NCBIGenbankFactory,
        scratch_repo: Repo,
        snapshot: SnapshotAssertion,
        taxid: int,
    ):
        """Test that a randomly generated list of records raises an appropriate
        ValueError.
        """
        otu = scratch_repo.get_otu_by_taxid(taxid)

        records = ncbi_genbank_factory.batch(len(otu.plan.required_segments))

        with pytest.raises(ValueError, match="Segment names not found in plan:") as e:
            assign_records_to_segments(records, otu.plan)

        assert str(e.value) == snapshot(exclude=props("id"))

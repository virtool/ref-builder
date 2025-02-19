import pytest
from syrupy import SnapshotAssertion

from ref_builder.otu.create import create_otu_with_taxid
from ref_builder.otu.update import auto_update_otu
from ref_builder.repo import Repo


@pytest.mark.ncbi()
class TestUpdateOTU:
    """Test automatic OTU update functionality."""

    def test_ok(
        self,
        precached_repo: Repo,
        snapshot: SnapshotAssertion,
    ):
        """Test automatic update behaviour."""
        with precached_repo.lock():
            otu_before = create_otu_with_taxid(
                precached_repo,
                2164102,
                ["NC_055390", "NC_055391", "NC_055392"],
                "",
            )

        assert otu_before.accessions == {"NC_055390", "NC_055391", "NC_055392"}

        assert otu_before.blocked_accessions == {
            "NC_055390",
            "NC_055391",
            "NC_055392",
            "MF062125",
            "MF062126",
            "MF062127",
        }

        with precached_repo.lock():
            auto_update_otu(precached_repo, otu_before)

        otu_after = precached_repo.get_otu(otu_before.id)

        assert otu_after.excluded_accessions == {"MF062125", "MF062126", "MF062127"}

        assert otu_after.id == otu_before.id

        assert otu_after.isolate_ids.issuperset(otu_before.isolate_ids)

        assert otu_after.accessions == {
            "MF062130",
            "MF062131",
            "MF062132",
            "MF062136",
            "MF062137",
            "MF062138",
            "NC_055390",
            "NC_055391",
            "NC_055392",
            "OR889795",
            "OR889796",
            "OR889797",
        }

        assert {
            str(isolate.name): isolate.accessions for isolate in otu_after.isolates
        } == snapshot()

    def test_with_refseq_replacement_ok(
        self,
        precached_repo: Repo,
        snapshot: SnapshotAssertion,
    ):
        """Test that automatic update replaces superceded accessions with RefSeq versions."""
        with precached_repo.lock():
            otu_before = create_otu_with_taxid(
                precached_repo,
                2164102,
                ["MF062125", "MF062126", "MF062127"],
                "",
            )

        assert (
            otu_before.accessions
            == otu_before.get_isolate(otu_before.representative_isolate).accessions
            == {"MF062125", "MF062126", "MF062127"}
        )

        with precached_repo.lock():
            auto_update_otu(precached_repo, otu_before)

        otu_after = precached_repo.get_otu(otu_before.id)

        assert otu_after.get_isolate(otu_after.representative_isolate).accessions == {
            "NC_055390",
            "NC_055391",
            "NC_055392",
        }
        assert {"MF062125", "MF062126", "MF062127"}.isdisjoint(otu_after.accessions)
        assert otu_after.representative_isolate == otu_before.representative_isolate
        assert (
            otu_after.get_isolate(otu_before.representative_isolate).accessions
            != otu_before.get_isolate(otu_before.representative_isolate).accessions
        )
        assert otu_after.id == otu_before.id
        assert otu_after.isolate_ids.issuperset(otu_before.isolate_ids)
        assert otu_after.excluded_accessions == {"MF062125", "MF062126", "MF062127"}

        assert otu_after.accessions == {
            "MF062136",
            "MF062137",
            "MF062138",
            "MF062130",
            "MF062131",
            "MF062132",
            "NC_055390",
            "NC_055391",
            "NC_055392",
            "OR889795",
            "OR889796",
            "OR889797",
        }

        assert {
            str(isolate.name): isolate.accessions for isolate in otu_after.isolates
        } == snapshot()

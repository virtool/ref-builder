import pytest
from syrupy import SnapshotAssertion

from ref_builder.otu.create import create_otu
from ref_builder.otu.update import auto_update_otu
from ref_builder.repo import Repo


@pytest.mark.ncbi()
class TestUpdateOTU:
    """Test automatic OTU update functionality."""

    def test_without_exclusions_ok(
        self,
        precached_repo: Repo,
        snapshot: SnapshotAssertion,
    ):
        """"""
        otu_before = create_otu(
            precached_repo,
            2164102,
            ["MF062136", "MF062137", "MF062138"],
            "",
        )

        assert otu_before.accessions == {"MF062136", "MF062137", "MF062138"}

        auto_update_otu(precached_repo, otu_before)

        otu_after = precached_repo.get_otu(otu_before.id)

        assert otu_after.id == otu_before.id

        assert otu_after.isolate_ids.issuperset(otu_before.isolate_ids)

        assert {
            isolate.name: isolate.accessions for isolate in otu_after.isolates
        } == snapshot()

        assert otu_after.excluded_accessions == {"MF062125", "MF062126", "MF062127"}

        assert otu_after.accessions == {
            "MF062136",
            "MF062137",
            "MF062138",
            "MF062130",
            "MF062131",
            "MF062132",
            "MK936225",
            "MK936226",
            "MK936227",
            "OQ420743",
            "OQ420744",
            "OQ420745",
            "NC_055390",
            "NC_055391",
            "NC_055392",
            "OR889795",
            "OR889796",
            "OR889797",
        }

    def test_with_replacement_ok(
        self,
        precached_repo: Repo,
        snapshot: SnapshotAssertion,
    ):
        otu_before = create_otu(
            precached_repo,
            2164102,
            ["MF062125", "MF062126", "MF062127"],
            "",
        )

        assert (
            otu_before.accessions
            == otu_before.get_isolate(otu_before.repr_isolate).accessions
            == {"MF062125", "MF062126", "MF062127"}
        )

        auto_update_otu(precached_repo, otu_before)

        otu_after = precached_repo.get_otu(otu_before.id)

        assert otu_after.get_isolate(otu_after.repr_isolate).accessions == {
            "NC_055390",
            "NC_055391",
            "NC_055392",
        }

        assert {"MF062125", "MF062126", "MF062127"}.isdisjoint(otu_after.accessions)

        assert otu_after.repr_isolate == otu_before.repr_isolate

        assert (
            otu_after.get_isolate(otu_before.repr_isolate).accessions
            != otu_before.get_isolate(otu_before.repr_isolate).accessions
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
            "MK936225",
            "MK936226",
            "MK936227",
            "OQ420743",
            "OQ420744",
            "OQ420745",
            "NC_055390",
            "NC_055391",
            "NC_055392",
            "OR889795",
            "OR889796",
            "OR889797",
        }

        assert {
            isolate.name: isolate.accessions for isolate in otu_after.isolates
        } == snapshot()

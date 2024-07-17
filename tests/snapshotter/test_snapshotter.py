import pytest

from ref_builder.repo import Repo
from ref_builder.snapshotter.snapshotter import Snapshotter


@pytest.fixture()
def snapshotter(scratch_repo: Repo) -> Snapshotter:
    return Snapshotter(scratch_repo.path / ".cache/snapshot")


class TestSnapshotIndex:
    def test_otu_ids(self, scratch_repo: Repo, snapshotter: Snapshotter):
        true_otu_ids = [otu.id for otu in scratch_repo.get_all_otus(ignore_cache=True)]

        assert snapshotter.otu_ids

        assert len(true_otu_ids) == len(snapshotter.otu_ids)

        for otu_id in true_otu_ids:
            assert otu_id in snapshotter.id_to_taxid

    def test_load_by_id(self, snapshotter: Snapshotter, scratch_repo: Repo):
        """Test that we can load an OTU by its ID."""
        otu_ids = [otu.id for otu in scratch_repo.get_all_otus(ignore_cache=True)]

        for otu_id in otu_ids:
            assert snapshotter.load_by_id(otu_id).id == otu_id

    def test_load_by_taxid(
        self,
        scratch_repo: Repo,
        snapshotter: Snapshotter,
    ):
        """Test that we can load an OTU by its taxid."""
        taxids = [otu.taxid for otu in scratch_repo.get_all_otus(ignore_cache=True)]

        for taxid in taxids:
            assert snapshotter.load_by_taxid(taxid).taxid == taxid

    def test_load_by_name(
        self,
        scratch_repo: Repo,
        snapshotter: Snapshotter,
    ):
        """Test that we can load an OTU by its name."""
        true_otu_names = [
            otu.name for otu in scratch_repo.get_all_otus(ignore_cache=True)
        ]

        for name in true_otu_names:
            assert snapshotter.load_by_name(name).name == name

    def test_accessions(self, scratch_repo: Repo, snapshotter: Snapshotter):
        true_accessions = set()
        for otu in scratch_repo.get_all_otus(ignore_cache=True):
            true_accessions.update(otu.accessions)

        assert true_accessions

        assert true_accessions == snapshotter.accessions


class TestSnapshotIndexCaching:
    @pytest.mark.parametrize(
        ("taxid", "accessions"),
        [
            (
                3158377,
                [
                    "NC_010314",
                    "NC_010315",
                    "NC_010316",
                    "NC_010317",
                    "NC_010318",
                    "NC_010319",
                ],
            ),
        ],
    )
    def test_load_otu_by_taxid(
        self,
        taxid: int,
        accessions: list[str],
        scratch_repo: Repo,
        snapshotter: Snapshotter,
    ):
        scratch_repo.snapshot()

        rehydrated_otu = scratch_repo.get_otu_by_taxid(taxid)

        assert rehydrated_otu

        snapshot_otu = snapshotter.load_by_taxid(taxid)

        assert snapshot_otu

        assert rehydrated_otu.accessions == snapshot_otu.accessions

import pytest

from ref_builder.repo import Repo
from ref_builder.snapshotter.snapshotter import Snapshotter


@pytest.fixture()
def snapshotter(scratch_repo: Repo) -> Snapshotter:
    return Snapshotter(scratch_repo.path / ".cache/snapshot")


class TestSnapshotIndex:
    def test_otu_ids(self, scratch_repo: Repo, snapshotter: Snapshotter):
        true_otu_ids = [otu.id for otu in scratch_repo.iter_otus(ignore_cache=True)]
        assert len(true_otu_ids) == len(snapshotter.otu_ids)

    def test_load_by_id(self, snapshotter: Snapshotter, scratch_repo: Repo):
        """Test that we can load an OTU by its ID."""
        otu_ids = [otu.id for otu in scratch_repo.iter_otus(ignore_cache=True)]

        for otu_id in otu_ids:
            assert snapshotter.load_by_id(otu_id).id == otu_id

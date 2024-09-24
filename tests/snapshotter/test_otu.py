from pathlib import Path
from uuid import UUID

import orjson
import pytest

from ref_builder.repo import Repo
from ref_builder.snapshotter.otu import OTUSnapshotter


@pytest.fixture()
def empty_snapshot_path(tmp_path) -> Path:
    return tmp_path / "otu_snapshot"


@pytest.mark.parametrize("taxid", [1441799, 430059])
class TestOTUSnapshot:
    def test_snapshot_direct(
        self,
        taxid: int,
        empty_snapshot_path: Path,
        scratch_repo: Repo,
    ):
        """Test that OTUSnapshot can build a snapshot directly from a floating OTU."""
        otu = scratch_repo.get_otu_by_taxid(taxid)

        snapshot = OTUSnapshotter(empty_snapshot_path)

        snapshot.cache(otu, 5)

        assert snapshot.at_event == 5

        for stem in ("otu", "toc"):
            assert (snapshot.path / f"{stem}.json").exists()

        data_path = snapshot.path / "data"

        assert data_path.exists()

        assert {isolate.id for isolate in otu.isolates}.issubset(
            {UUID(path.stem) for path in data_path.glob("*.json")},
        )

    def test_load(self, taxid: int, scratch_repo: Repo):
        """Test OTUSnapshot.load()"""
        otu = scratch_repo.get_otu_by_taxid(taxid)

        snapshot = OTUSnapshotter(
            path=scratch_repo.cache_path / f"snapshot/{otu.id}",
        )

        assert otu.accessions == snapshot.load().accessions

    def test_toc(self, taxid: int, scratch_repo: Repo):
        """Test that the table of contents is written correctly."""
        rehydrated_otu = scratch_repo.get_otu_by_taxid(taxid)

        otu_snapshotter = OTUSnapshotter(
            path=scratch_repo.cache_path / f"snapshot/{rehydrated_otu.id}",
        )

        with open((otu_snapshotter.path / "toc.json"), "rb") as f:
            toc_dict = orjson.loads(f.read())

        for isolate in rehydrated_otu.isolates:
            for accession in isolate.accessions:
                assert accession in toc_dict[str(isolate.id)]["accessions"]

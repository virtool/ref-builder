from pathlib import Path
from uuid import uuid4

import pytest
from polyfactory.factories import DataclassFactory

from ref_builder.repo import Repo
from ref_builder.resources import RepoIsolate, RepoOTU, RepoSequence
from ref_builder.snapshotter.snapshotter import Snapshotter
from ref_builder.utils import Accession


@pytest.fixture()
def snapshotter(scratch_repo: Repo) -> Snapshotter:
    return Snapshotter(scratch_repo.path / ".cache/snapshot")


def test_ok(tmp_path: Path):
    snapshotter = Snapshotter(tmp_path / "snapshotter")

    otu = RepoOTU(
        uuid4(),
        12345,
        "Made Up Virus",
        schema,
    )

    class SequenceFactory(DataclassFactory[RepoSequence]):
        __model__ = RepoSequence

    sequences = [
        SequenceFactory.build(
            accession=Accession("ABCDEFA", 1),
        ),
        SequenceFactory.build(
            accession=Accession("ABCDEFB", 1),
        ),
    ]

    class IsolateFactory(DataclassFactory[RepoIsolate]):
        __model__ = RepoIsolate

    isolate = IsolateFactory.build(sequences=sequences)

    class OTUFactory(DataclassFactory[RepoOTU]):
        __model__ = RepoOTU

    otu = OTUFactory.build(
        isolates=[isolate],
    )


class TestIndex:
    def test_otu_ids(self, scratch_repo: Repo, snapshotter: Snapshotter):
        true_otu_ids = [otu.id for otu in scratch_repo.iter_otus(ignore_cache=True)]
        assert len(true_otu_ids) == len(snapshotter.otu_ids)

    def test_load_by_id(self, snapshotter: Snapshotter, scratch_repo: Repo):
        """Test that we can load an OTU by its ID."""
        otu_ids = [otu.id for otu in scratch_repo.iter_otus(ignore_cache=True)]

        for otu_id in otu_ids:
            assert snapshotter.load_by_id(otu_id).id == otu_id

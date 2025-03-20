import datetime
from pathlib import Path

import pytest
from syrupy import SnapshotAssertion

from ref_builder.otu.create import create_otu_with_taxid
from ref_builder.otu.isolate import add_genbank_isolate
from ref_builder.otu.update import (
    BatchFetchIndex,
    auto_update_otu,
    batch_update_repo,
    iter_fetch_list,
    promote_otu_accessions,
)
from ref_builder.repo import Repo


@pytest.fixture()
def mock_repo(precached_repo: Repo) -> Repo:
    with precached_repo.lock():
        create_otu_with_taxid(
            precached_repo,
            2164102,
            ["MF062125", "MF062126", "MF062127"],
            "",
        )

    yield precached_repo


@pytest.fixture()
def mock_fetch_index() -> dict[int, set[str]]:
    """A mock fetch index for NCBI Taxonomy ID 2164102."""
    return {
        2164102: {
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
    }


@pytest.fixture()
def mock_fetch_index_path(
    mock_fetch_index: dict[int, set[str]], tmp_path: Path
) -> Path:
    """A temporary path to a mock fetch index."""
    file_path = tmp_path / "fetch_index_2165102.json"

    fetch_index = BatchFetchIndex.model_validate(mock_fetch_index)

    with open(file_path, "w") as f:
        f.write(fetch_index.model_dump_json())

    yield file_path

    file_path.unlink()


@pytest.mark.ncbi()
class TestPromoteOTU:
    """Test OTU accession promotion from Genbank to RefSeq."""

    def test_ok(self, empty_repo: Repo):
        """Test that RefSeq accessions can be promoted automatically."""
        with empty_repo.lock():
            otu = create_otu_with_taxid(
                empty_repo, 2164102, ["MF062136", "MF062137", "MF062138"], acronym=""
            )

            isolate = add_genbank_isolate(
                empty_repo, otu, ["MF062125", "MF062126", "MF062127"]
            )

        otu_before = empty_repo.get_otu(otu.id)

        assert otu_before.accessions == {
            "MF062125",
            "MF062126",
            "MF062127",
            "MF062136",
            "MF062137",
            "MF062138",
        }

        assert otu_before.get_isolate(isolate.id).accessions == {
            "MF062125",
            "MF062126",
            "MF062127",
        }

        with empty_repo.lock():
            promoted_accessions = promote_otu_accessions(empty_repo, otu_before)

        assert promoted_accessions == {"NC_055390", "NC_055391", "NC_055392"}

        otu_after = empty_repo.get_otu(otu.id)

        assert otu_after.isolate_ids == otu_before.isolate_ids

        assert otu_after.accessions == {
            "NC_055390",
            "NC_055391",
            "NC_055392",
            "MF062136",
            "MF062137",
            "MF062138",
        }

        assert otu_after.get_isolate(isolate.id).accessions == {
            "NC_055390",
            "NC_055391",
            "NC_055392",
        }

        assert otu_after.excluded_accessions == {"MF062125", "MF062126", "MF062127"}


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

    def test_start_date_limit(self, precached_repo: Repo):
        """Test automatic update with the start date set to ``today``."""
        with precached_repo.lock():
            otu_before = create_otu_with_taxid(
                precached_repo,
                2164102,
                ["NC_055390", "NC_055391", "NC_055392"],
                "",
            )

        assert otu_before.accessions == {"NC_055390", "NC_055391", "NC_055392"}

        with precached_repo.lock():
            otu_after = auto_update_otu(
                precached_repo,
                otu_before,
                start_date=datetime.date.today(),
            )

        assert {
            "MF062130",
            "MF062131",
            "MF062132",
            "MF062136",
            "MF062137",
            "MF062138",
            "OR889795",
            "OR889796",
            "OR889797",
        }.isdisjoint(otu_after.accessions)

    def test_with_refseq_replacement_ok(
        self,
        precached_repo: Repo,
        mock_fetch_index: dict[int, set[str]],
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
            otu_after = auto_update_otu(precached_repo, otu_before)

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

        assert otu_after.accessions == mock_fetch_index[2164102]

        assert {
            str(isolate.name): isolate.accessions for isolate in otu_after.isolates
        } == snapshot()

    def test_with_fetch_index_ok(
        self,
        mock_repo: Repo,
        mock_fetch_index: dict[int, set[str]],
        mock_fetch_index_path: Path,
    ):
        """Test with a path to a pre-made fetch index as input."""
        otu_initial = next(mock_repo.iter_otus())

        assert otu_initial

        with mock_repo.lock():
            otu_after = auto_update_otu(
                mock_repo, otu_initial, fetch_index_path=mock_fetch_index_path
            )

        assert otu_after.accessions == mock_fetch_index[2164102]


class TestBatchUpdate:
    """Test rudimentary batch update operation with a single OTU."""

    def test_ok(
        self,
        mock_repo: Repo,
        mock_fetch_index: dict[int, set[str]],
    ):
        """Test that batch update works as expected."""
        otu_initial = next(mock_repo.iter_otus())

        with mock_repo.lock():
            assert len(batch_update_repo(mock_repo)) == 1

        otu_after = next(mock_repo.iter_otus())

        assert otu_after.accessions == mock_fetch_index[otu_initial.taxid]

        with mock_repo.lock():
            assert len(batch_update_repo(mock_repo)) == 0

    def test_with_fetch_index_ok(
        self,
        mock_repo: Repo,
        mock_fetch_index: dict[int, set[str]],
        mock_fetch_index_path: Path,
    ):
        """Test with a path to a pre-made fetch index as input."""
        otu_initial = next(mock_repo.iter_otus())

        with mock_repo.lock():
            assert len(
                batch_update_repo(mock_repo, fetch_index_path=mock_fetch_index_path)
            ) == 1

        otu_after = next(mock_repo.iter_otus())

        assert otu_after.accessions == mock_fetch_index[otu_initial.taxid]

        with mock_repo.lock():
            assert len(
                batch_update_repo(mock_repo, fetch_index_path=mock_fetch_index_path)
            ) == 0


def test_iter_fetch_list(
    mock_fetch_index: dict[int, set[str]],
):
    """Test fetch list iterator."""
    fetch_list = list(mock_fetch_index[2164102])

    for page_size in [1, 3, 5, 20]:
        for fetch_list_segment in iter_fetch_list(fetch_list, page_size):
            assert len(fetch_list_segment)

        regenerated_fetch_list = []
        for chunk in iter_fetch_list(fetch_list, page_size):
            regenerated_fetch_list.extend(chunk)

        assert fetch_list == regenerated_fetch_list

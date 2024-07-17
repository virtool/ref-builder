import shutil
from pathlib import Path

import pytest

from ref_builder.legacy.utils import (
    build_legacy_otu,
    generate_unique_ids,
    iter_legacy_otus,
    replace_otu,
)


def test_iter_legacy_otu(legacy_repo_path: Path):
    """Test that all OTUs in the repository are iterated over."""
    list(iter_legacy_otus(legacy_repo_path / "src"))


class TestReplaceOTU:
    """Test the replace_otu() function."""

    def test_ok(self, legacy_otu: dict, legacy_repo_path: Path):
        """Test that the OTU is replaced with the new data. Only the changed fields should
        be changed in the repository src path.
        """
        legacy_otu["name"] = "Fake OTU"
        legacy_otu["isolates"][1].update(
            {"source_type": "strain", "source_name": "Fake"},
        )

        otu_path = legacy_repo_path / "src" / "a" / "abaca_bunchy_top_virus"

        replace_otu(otu_path, legacy_otu)

        assert build_legacy_otu(otu_path) == legacy_otu

    def test_bad_otu_path(self, legacy_otu: dict, legacy_repo_path: Path):
        """Test that a FileNotFoundError is raised when the otu.json file is not found."""
        legacy_otu["name"] = "Fake OTU"
        legacy_otu["isolates"][1].update(
            {"source_type": "strain", "source_name": "Fake"},
        )

        otu_path = legacy_repo_path / "src" / "a" / "abaca_bunchy_bottom_virus"

        with pytest.raises(FileNotFoundError) as e:
            replace_otu(
                otu_path,
                legacy_otu,
            )

        assert str(e.value) == f"No such OTU directory: {otu_path}"

    @pytest.mark.parametrize("path", ["", "src", "src/a"])
    def test_bad_otu_parent_path(
        self,
        path: Path,
        legacy_otu: dict,
        legacy_repo_path: Path,
    ):
        """Test that a FileNotFoundError is raised when a parent directory of the
        OTU path does not exist.
        """
        legacy_otu["name"] = "Fake OTU"
        legacy_otu["isolates"][1].update(
            {"source_type": "strain", "source_name": "Fake"},
        )

        otu_path = legacy_repo_path / "src" / "a" / "abaca_bunchy_top_virus"

        shutil.rmtree(legacy_repo_path / path)

        with pytest.raises(FileNotFoundError) as e:
            replace_otu(
                otu_path,
                legacy_otu,
            )

        assert str(e.value) == f"No such OTU parent directory: {otu_path.parent}"


def test_generate_unique_ids():
    """Test that the generate_unique_ids() function generates unique IDs."""
    initial_hashes = generate_unique_ids(n=40, length=8, mixed_case=False, excluded=[])

    additional_hashes = generate_unique_ids(
        n=20,
        length=8,
        mixed_case=False,
        excluded=initial_hashes,
    )

    assert initial_hashes.isdisjoint(additional_hashes)

import json
from io import StringIO
from pathlib import Path

import pytest
from pytest_mock import MockerFixture
from rich.console import Console

from ref_builder.legacy.checks import (
    check_isolate_names_against_ncbi,
    check_unique_otu_abbreviations_and_names,
    find_duplicate_accessions,
)


@pytest.fixture()
def _console(mocker: MockerFixture) -> Console:
    """The console object in the legacy console module.

    Patched to use a StringIO file instead of stdout.
    """
    console = Console(file=StringIO())
    mocker.patch("ref_builder.legacy.repo.console", console)

    return console


def test_check_isolate_names_against_ncbi(legacy_repo_path: Path):
    """Test that check_isolate_names_against_ncbi() reports isolates with names that
    don't match what is found in NCBI for their member accessions.
    """
    assert check_isolate_names_against_ncbi(legacy_repo_path) == {}
    assert 0


def test_find_duplicate_accessions(
    legacy_otu: dict,
    legacy_repo_path: Path,
):
    """Test that find_duplicate_accessions() finds non-unique accessions and returns
    the sequences they are in.
    """
    sequence_path = (
        legacy_repo_path
        / "src"
        / "b"
        / "bamboo_mosaic_virus"
        / "lm3xmhrf"
        / "3d2e72sk.json"
    )

    with open(sequence_path) as f:
        data = json.load(f)

    duplicate_accession = legacy_otu["isolates"][0]["sequences"][0]["accession"]

    with open(sequence_path, "w") as f:
        json.dump(
            {
                **data,
                "accession": duplicate_accession,
            },
            f,
        )

    assert find_duplicate_accessions(legacy_repo_path) == {
        "NC_010317.1": [
            "3d2e72sk",
            "5w3oz9pl",
        ],
    }


class TestCheckUniqueOTUNamesAndAbbreviations:
    def test_ok(
        self,
        _console: Console,
        legacy_repo_path: Path,
    ):
        """Test that validation passes when all OTU abbreviations and names are
        unique.
        """
        check_unique_otu_abbreviations_and_names(legacy_repo_path)
        assert _console.file.getvalue() == ""

    @pytest.mark.parametrize("abbreviation", ["ABTV", "BaMV"])
    @pytest.mark.parametrize(
        "name",
        ["Abaca bunchy top virus", "Bamboo mosaic virus"],
    )
    def test_duplicates(
        self,
        abbreviation: str,
        name: str,
        _console: Console,
        legacy_repo_path: Path,
    ):
        """Test that check_unique_otu_abbreviations_and_names() finds non-unique OTU
        abbreviations.
        """
        json_path = legacy_repo_path / "src" / "b" / "bamboo_mosaic_virus" / "otu.json"

        with open(json_path) as f:
            data = json.load(f)

        with open(json_path, "w") as f:
            json.dump(
                {
                    **data,
                    "abbreviation": abbreviation,
                    "name": name,
                },
                f,
            )

        check_unique_otu_abbreviations_and_names(legacy_repo_path)

        logs = _console.file.getvalue()

        if abbreviation == "ABTV":
            assert "Found non-unique OTU abbreviation: ABTV" in logs

        if name == "Abaca bunchy top virus":
            assert "Found non-unique OTU name: abaca bunchy top virus" in logs

        if abbreviation == "BaMV" and name == "Bamboo mosaic virus":
            assert logs == ""

    def test_empty_abbreviation(
        self,
        _console: Console,
        legacy_repo_path: Path,
    ):
        """Test that check passes for duplicate empty abbreviations."""
        for path in (Path("b/bamboo_mosaic_virus"), Path("o/oat_blue_dwarf_virus")):
            json_path = legacy_repo_path / "src" / path / "otu.json"

            with open(json_path) as f:
                data = json.load(f)

            with open(json_path, "w") as f:
                json.dump(
                    {
                        **data,
                        "abbreviation": "",
                    },
                    f,
                )

        check_unique_otu_abbreviations_and_names(legacy_repo_path)

        assert _console.file.getvalue() == ""

import subprocess
from pathlib import Path

import pytest
from click.testing import CliRunner
from syrupy import SnapshotAssertion
from syrupy.filters import props

from ref_builder.cli.otu import otu as otu_command_group
from ref_builder.console import console, print_otu
from ref_builder.repo import Repo


runner = CliRunner()


class TestCreateOTUCommands:
    """Test the behaviour of ``ref-builder otu create``."""

    @pytest.mark.parametrize(
        ("taxid", "accessions"),
        [(1278205, ["NC_020160"]), (345184, ["DQ178610", "DQ178611"])],
    )
    def test_with_taxid_ok(
        self,
        taxid: int,
        accessions: list[str],
        precached_repo: Repo,
        snapshot: SnapshotAssertion,
    ):
        """Test that an OTU can be created with the --taxid option.

        Also check resulting print_otu() console output.
        """
        result = runner.invoke(
            otu_command_group,
            ["--path", str(precached_repo.path)]
            + ["create", "--taxid", str(taxid)]
            + accessions,
        )

        assert result.exit_code == 0

        otus = list(Repo(precached_repo.path).iter_otus())

        assert len(otus) == 1
        assert otus[0].model_dump() == snapshot(
            exclude=props("id", "isolates", "representative_isolate"),
        )

        with console.capture() as capture:
            print_otu(otus[0])

        assert capture.get() in result.output

    @pytest.mark.parametrize(
        ("taxid", "accessions"),
        [(1278205, ["NC_020160"]), (345184, ["DQ178610", "DQ178611"])],
    )
    def test_without_taxid_ok(
        self,
        taxid: int,
        accessions: list[str],
        precached_repo: Repo,
    ):
        """Test that an OTU can be created without the --taxid option."""
        result = runner.invoke(
            otu_command_group,
            ["--path", str(precached_repo.path)] + ["create"] + accessions,
        )

        assert result.exit_code == 0

        otus = list(Repo(precached_repo.path).iter_otus())

        assert len(otus) == 1
        assert otus[0].taxid == taxid

    def test_acronym(self, precached_repo: Repo):
        """Test that Athe --acronym option works as planned."""
        acronym_ = "CabLCJV"

        result = runner.invoke(
            otu_command_group,
            ["--path", str(precached_repo.path)]
            + [
                "create",
                "--taxid",
                str(345184),
                "--acronym",
                acronym_,
            ]
            + ["DQ178610", "DQ178611"],
        )

        assert result.exit_code == 0

        otus = list(Repo(precached_repo.path).iter_otus())

        assert len(otus) == 1
        assert otus[0].acronym == acronym_

    def test_duplicate_accessions(self, precached_repo: Repo):
        """Test that an error is raised when duplicate accessions are provided."""
        runner = CliRunner()

        result = runner.invoke(
            otu_command_group,
            ["--path", str(precached_repo.path)]
            + ["create", "--taxid", str(1169032)]
            + ["MK431779", "MK431779"],
        )

        assert result.exit_code == 2
        assert "Duplicate accessions are not allowed." in result.output

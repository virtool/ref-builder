from uuid import UUID

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


@pytest.mark.ncbi
class TestPromoteOTUCommand:
    """Test that the ``ref-builder otu promote`` command works as planned."""

    def test_promotable_ok(self, empty_repo: Repo):
        """Test that otu promote adds new accessions to OTU."""
        path_option = ["--path", str(empty_repo.path)]

        taxid = 2164102

        original_rep_isolate_accessions = {"MF062125", "MF062126", "MF062127"}

        result = runner.invoke(
            otu_command_group,
            path_option
            + ["create", "--taxid", str(taxid)]
            + list(original_rep_isolate_accessions),
        )

        assert result.exit_code == 0

        otu_before = empty_repo.get_otu_by_taxid(taxid)

        assert otu_before.accessions == original_rep_isolate_accessions

        result = runner.invoke(
            otu_command_group,
            path_option + ["promote", str(taxid)],
        )

        assert result.exit_code == 0

        assert (
            "Sequences promoted"
            and "['NC_055390', 'NC_055391', 'NC_055392']" in result.output
        )

        repo_after = Repo(empty_repo.path)

        otu_after = repo_after.get_otu(otu_before.id)

        assert otu_after.representative_isolate == otu_before.representative_isolate

        assert otu_after.accessions == {"NC_055390", "NC_055391", "NC_055392"}

        assert otu_after.excluded_accessions == original_rep_isolate_accessions

    def test_redundant_ok(self, empty_repo: Repo):
        """Test that command works correctly when the OTU is already up to date."""
        path_option = ["--path", str(empty_repo.path)]

        taxid = 2164102

        rep_isolate_accessions = {"NC_055390", "NC_055391", "NC_055392"}

        result = runner.invoke(
            otu_command_group,
            path_option
            + ["create", "--taxid", str(taxid)]
            + list(rep_isolate_accessions),
        )

        assert result.exit_code == 0

        otu = empty_repo.get_otu_by_taxid(taxid)

        assert otu.excluded_accessions == {"MF062125", "MF062126", "MF062127"}

        result = runner.invoke(
            otu_command_group,
            path_option + ["promote", str(taxid)],
        )

        assert result.exit_code == 0

        assert "Isolate updated" not in result.output


class TestExcludeAccessionsCommand:
    """Test that ``ref-builder otu exclude-accessions`` behaves as expected."""

    def test_excludable_ok(self, scratch_repo):
        """Test that command lists out new excluded accessions"""
        result = runner.invoke(
            otu_command_group,
            ["--path", str(scratch_repo.path)]
            + ["exclude-accessions", str(345184)]
            + ["DQ178608", "DQ178609"],
        )

        assert result.exit_code == 0

        assert "Added accessions to excluded accession list" in result.output

        assert "['DQ178608', 'DQ178609']" in result.output

    def test_redundant_ok(self, scratch_repo):
        """Test that the command informs the user when excluded accessions are already up to date."""
        result = runner.invoke(
            otu_command_group,
            ["--path", str(scratch_repo.path)]
            + ["exclude-accessions", str(345184)]
            + ["DQ178608", "DQ178609"],
        )

        assert result.exit_code == 0

        assert "['DQ178608', 'DQ178609']" in result.output

        result = runner.invoke(
            otu_command_group,
            ["--path", str(scratch_repo.path)]
            + ["exclude-accessions", str(345184)]
            + ["DQ178608", "DQ178609"],
        )

        assert result.exit_code == 0

        assert "Excluded accession list already up to date" in result.output


class TestAllowAccessionsCommand:
    """Test that ``ref-builder otu allow-accessions`` behaves as expected."""

    def test_allowable_ok(self, scratch_repo: Repo):
        """Test that command lists out newly allowable accessions"""
        taxid = 345184

        result = runner.invoke(
            otu_command_group,
            ["--path", str(scratch_repo.path)]
            + ["exclude-accessions", str(taxid)]
            + ["DQ178608", "DQ178609"],
        )

        assert result.exit_code == 0

        result = runner.invoke(
            otu_command_group,
            ["--path", str(scratch_repo.path)]
            + ["allow-accessions", str(taxid)]
            + ["DQ178608"],
        )

        assert result.exit_code == 0

        assert "Removed accessions from excluded accession list" in result.output

        assert "['DQ178608']" in result.output

        assert "Updated excluded accession list" in result.output

        assert "['DQ178609']" in result.output

    def test_redundant_ok(self, scratch_repo: Repo):
        """Test that the command informs the user when excluded accessions are already up to date."""
        result = runner.invoke(
            otu_command_group,
            ["--path", str(scratch_repo.path)]
            + ["allow-accessions", str(345184)]
            + ["DQ178612", "DQ178613"],
        )

        assert result.exit_code == 0

        assert "Excluded accession list already up to date" in result.output

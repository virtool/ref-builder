from click.testing import CliRunner

from ref_builder.utils import IsolateName, IsolateNameType
from ref_builder.cli.otu import otu as otu_command_group
from ref_builder.cli.isolate import isolate as isolate_command_group

runner = CliRunner()


class TestIsolateCreateCommand:
    """Test `ref-builder isolate create`` works as expected."""

    def test_ok(self, precached_repo):
        """Test basic command functionality."""
        path_option_list = ["--path", str(precached_repo.path)]

        taxid = 1169032

        rep_isolate_accessions = ["MF062136", "MF062137", "MF062138"]

        result = runner.invoke(
            otu_command_group,
            path_option_list
            + ["create", "--taxid", str(taxid)]
            + rep_isolate_accessions,
        )

        assert result.exit_code == 0

        second_isolate_accessions = ["MF062125", "MF062126", "MF062127"]

        result = runner.invoke(
            isolate_command_group,
            path_option_list
            + ["create", "--taxid", str(taxid)]
            + second_isolate_accessions,
        )

        assert result.exit_code == 0

        assert "Isolate created" in result.output

    def test_overwrite_name_option_ok(self, precached_repo):
        """Test that --name option exits smoothly"""
        path_option_list = ["--path", str(precached_repo.path)]

        taxid = 345184

        result = runner.invoke(
            otu_command_group,
            path_option_list
            + ["create", "--taxid", str(taxid)]
            + ["DQ178610", "DQ178611"],
        )

        assert result.exit_code == 0

        result = runner.invoke(
            isolate_command_group,
            path_option_list
            + ["create", "--taxid", str(taxid)]
            + ["--name", "isolate", "dummy"]
            + ["DQ178613", "DQ178614"],
        )

        assert result.exit_code == 0

        assert "Isolate created" and "Isolate dummy" in result.output

    def test_unnamed_option_ok(self, precached_repo):
        """Test that --unnamed option exits smoothly."""
        path_option_list = ["--path", str(precached_repo.path)]

        taxid = 345184

        result = runner.invoke(
            otu_command_group,
            path_option_list
            + ["create", "--taxid", str(taxid)]
            + ["DQ178610", "DQ178611"],
        )

        assert result.exit_code == 0

        result = runner.invoke(
            isolate_command_group,
            path_option_list
            + ["create", "--taxid", str(taxid)]
            + ["--unnamed"]
            + ["DQ178613", "DQ178614"],
        )

        assert result.exit_code == 0

        assert "Isolate created" and "Unnamed" in result.output

    def test_duplicate_accessions_error(self, scratch_repo):
        """Test that an error is raised when duplicate accessions are provided."""

        result = runner.invoke(
            isolate_command_group,
            ["--path", str(scratch_repo.path)]
            + ["create", "--taxid", str(345184)]
            + ["DQ178610", "DQ178610"],
        )

        assert result.exit_code == 2

        assert "Duplicate accessions are not allowed." in result.output


class TestIsolateDeleteCommand:
    """Test `ref-builder isolate delete ISOLATE_ID`` works as expected."""

    def test_ok(self, scratch_repo):
        """Test basic command functionality."""
        otu = scratch_repo.get_otu_by_taxid(1169032)

        isolate_id = otu.get_isolate_id_by_name(
            IsolateName(type=IsolateNameType.ISOLATE, value="WMoV-6.3"),
        )

        assert isolate_id in otu.isolate_ids

        result = runner.invoke(
            isolate_command_group,
            ["--path", str(scratch_repo.path), "delete", str(isolate_id)],
        )

        assert result.exit_code == 0

        assert "Isolate deleted" in result.output

        assert isolate_id not in scratch_repo.get_otu_by_taxid(1169032).isolate_ids

    def test_with_partial_id_ok(self, scratch_repo):
        """Test that a partial isolate ID can also be a valid identifier."""

        otu = scratch_repo.get_otu_by_taxid(1169032)

        isolate_id = otu.get_isolate_id_by_name(
            IsolateName(type=IsolateNameType.ISOLATE, value="WMoV-6.3"),
        )

        assert isolate_id in otu.isolate_ids

        result = runner.invoke(
            isolate_command_group,
            ["--path", str(scratch_repo.path), "delete", str(isolate_id)[:8]],
        )

        assert result.exit_code == 0

        assert "Isolate deleted" in result.output

        assert isolate_id not in scratch_repo.get_otu_by_taxid(1169032).isolate_ids

    def test_empty(self, scratch_repo):
        """Test that an empty isolate identifier string exits with an error."""
        result = runner.invoke(
            isolate_command_group,
            ["--path", str(scratch_repo.path), "delete", ""],
        )

        assert result.exit_code == 1

from click.testing import CliRunner

from ref_builder.utils import IsolateName, IsolateNameType
from ref_builder.cli.isolate import isolate as isolate_command_group

runner = CliRunner()


class TestIsolateGetCommand:
    """Test `ref-builder isolate get ISOLATE_ID`` works as expected."""

    def test_ok(self, scratch_repo):
        """Test basic command functionality."""
        otu_id = scratch_repo.get_otu_id_by_taxid(1169032)

        for isolate_id in scratch_repo.get_otu(otu_id).isolate_ids:
            result = runner.invoke(
                isolate_command_group,
                [
                    "--path",
                    str(scratch_repo.path),
                    "get",
                    str(isolate_id),
                ],
            )

            assert result.exit_code == 0

    def test_partial_ok(self, scratch_repo):
        """Test that a partial isolate ID can also be a valid identifier."""
        otu_id = scratch_repo.get_otu_id_by_taxid(1169032)

        for isolate_id in scratch_repo.get_otu(otu_id).isolate_ids:
            result = runner.invoke(
                isolate_command_group,
                [
                    "--path",
                    str(scratch_repo.path),
                    "get",
                    str(isolate_id)[:8],
                ],
            )

            assert result.exit_code == 0

    def test_json_ok(self, scratch_repo):
        otu_id = scratch_repo.get_otu_id_by_taxid(1169032)

        for isolate_id in scratch_repo.get_otu(otu_id).isolate_ids:
            result = runner.invoke(
                isolate_command_group,
                ["--path", str(scratch_repo.path), "get", str(isolate_id), "--json"],
            )

            assert result.exit_code == 0

    def test_empty(self, scratch_repo):
        """Test that an empty isolate identifier string exits with an error."""
        result = runner.invoke(
            isolate_command_group,
            [
                "--path",
                str(scratch_repo.path),
                "get",
                "",
            ],
        )

        assert result.exit_code == 1

        assert "Isolate ID partial length < 8." in result.output


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

        assert "Isolate ID partial length < 8." in result.output

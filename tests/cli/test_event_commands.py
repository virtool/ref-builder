from click.testing import CliRunner

from ref_builder.cli.event import event as event_command_group
from ref_builder.repo import Repo

runner = CliRunner()


class TestEventListCommand:
    """Test basic ``ref-builder event list`` command functionality."""

    def test_event_list_command(self, scratch_repo: Repo):
        result = runner.invoke(
            event_command_group, ["--path", str(scratch_repo.path), "list"]
        )

        assert result.exit_code == 0

        assert result.output

    def test_event_list_with_otu_identifier(self, scratch_repo: Repo):
        """Test ``event list --otu IDENTIFIER`` functionality."""
        for minimal_otu in scratch_repo.iter_minimal_otus():
            result = runner.invoke(
                event_command_group,
                [
                    "--path",
                    str(scratch_repo.path),
                    "list",
                    "--otu",
                    str(minimal_otu.id),
                ],
            )

            assert result.exit_code == 0

            assert result.output


class TestEventGetCommand:
    """Test basic ``ref-builder event get IDENTIFIER`` command functionality."""

    def test_event_get_command(self, scratch_repo: Repo):
        for event_metadata in scratch_repo.iter_event_metadata():
            result = runner.invoke(
                event_command_group,
                ["--path", str(scratch_repo.path), "get", str(event_metadata.id)],
            )

            assert result.exit_code == 0

            assert result.output

    def test_event_get_nonexistent_command(self, scratch_repo: Repo):
        """Test that a nonexistent event prints error to console."""
        result = runner.invoke(
            event_command_group,
            ["--path", str(scratch_repo.path), "get", str(scratch_repo.last_id + 5)],
        )

        assert result.exit_code == 1

        assert "Event not found" in result.output

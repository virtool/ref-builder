from click.testing import CliRunner

from ref_builder.cli.main import entry as top_command
from ref_builder.repo import Repo

runner = CliRunner()


class TestCheckCommand:
    """Check that ``ref-builder check`` behaves as expected."""

    def test_valid_repo(scratch_repo: Repo):
        """Check that ``ref-builder check`` returns expected log on a valid Repo."""
        result = runner.invoke(
            top_command,
            [
                "--verbose",
                "check",
                "--path",
                str(scratch_repo.path),
            ]
        )

        assert result.exit_code == 0

        assert "Repo is clean" in result.output

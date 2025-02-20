from click.testing import CliRunner

from ref_builder.cli.otu import otu as otu_command_group


runner = CliRunner()


class TestPrintCommands:
    """Test otu console printing commands."""

    def test_otu_list(self, scratch_repo):
        result = runner.invoke(
            otu_command_group,
            ["--path", str(scratch_repo.path)] + ["list"],
        )

        assert result.exit_code == 0

    def test_otu_get_ok(self, scratch_repo):
        for otu_ in scratch_repo.iter_otus():
            result = runner.invoke(
                otu_command_group,
                ["--path", str(scratch_repo.path)] + ["get", str(otu_.id)],
            )

            assert result.exit_code == 0

    def test_otu_get_json_ok(self, scratch_repo):
        for otu_ in scratch_repo.iter_otus():
            result = runner.invoke(
                otu_command_group,
                ["--path", str(scratch_repo.path)]
                + ["get", str(otu_.id), "--json"],
            )

            assert result.exit_code == 0

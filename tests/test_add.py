import subprocess
from pathlib import Path

import pytest
from click.testing import CliRunner
from syrupy import SnapshotAssertion
from syrupy.filters import props

from ref_builder.cli.main import isolate_create, otu_create
from ref_builder.otu.create import create_otu
from ref_builder.otu.update import (
    add_and_name_isolate,
    add_genbank_isolate,
    add_unnamed_isolate,
    update_isolate_from_accessions,
)
from ref_builder.otu.utils import RefSeqConflictError
from ref_builder.repo import Repo
from ref_builder.utils import IsolateName, IsolateNameType


def run_create_otu_command(
    path: Path,
    taxid: int,
    accessions: list,
    acronym: str = "",
):
    subprocess.run(
        [
            "ref-builder",
            "otu",
            "create",
            str(taxid),
            *accessions,
            "--path",
            str(path),
            "--acronym",
            acronym,
        ],
        check=False,
    )


def run_update_otu_command(taxid: int, path: Path):
    subprocess.run(
        ["ref-builder", "otu", "update"] + [str(taxid)] + ["--path", str(path)],
        check=False,
    )


class TestCreateOTU:
    def test_duplicate_accessions(self, precached_repo: Repo):
        """Test that an error is raised when duplicate accessions are provided."""
        runner = CliRunner()
        result = runner.invoke(
            otu_create,
            [
                "--path",
                str(precached_repo.path),
                "1169032",
                "MK431779",
                "MK431779",
            ],
        )

        assert result.exit_code == 2
        assert "Duplicate accessions are not allowed." in result.output

    def test_empty_success(
        self,
        precached_repo: Repo,
        snapshot: SnapshotAssertion,
    ):
        """Test that an OTU can be created in an empty repository."""
        otu = create_otu(
            precached_repo,
            345184,
            ["DQ178610", "DQ178611"],
            "",
        )

        assert otu.model_dump() == snapshot(
            exclude=props("id", "isolates", "repr_isolate"),
        )

        # Ensure only one OTU is present in the repository, and it matches the return
        # value of the creation function.
        assert list(precached_repo.iter_otus()) == [otu]

    def test_empty_fail(self, scratch_repo: Repo):
        with pytest.raises(ValueError):
            create_otu(
                scratch_repo,
                345184,
                ["DQ178610", "DQ178611"],
                "",
            )

    @pytest.mark.parametrize(
        ("taxid", "accessions"),
        [(1278205, ["NC_020160"]), (345184, ["DQ178610", "DQ178611"])],
    )
    def test_otu_create(
        self,
        taxid: int,
        accessions: list[str],
        precached_repo: Repo,
    ):
        assert list(precached_repo.iter_otus()) == []

        otu = create_otu(
            precached_repo,
            taxid,
            accessions,
            "",
        )

        assert list(precached_repo.iter_otus())
        assert otu.plan is not None
        assert otu.repr_isolate is not None

    def test_otu_create_refseq_autoexclude(
        self,
        precached_repo: Repo,
        snapshot: SnapshotAssertion,
    ):
        """Test that the superceded accessions included in RefSeq metadata
        are automatically added to the OTU's excluded accessions list.
        """
        otu = create_otu(
            precached_repo,
            3158377,
            [
                "NC_010314",
                "NC_010316",
                "NC_010315",
                "NC_010317",
                "NC_010318",
                "NC_010319",
            ],
            "",
        )

        assert (
            otu.excluded_accessions
            == {
                "EF546808",
                "EF546809",
                "EF546810",
                "EF546811",
                "EF546812",
                "EF546813",
            }
            == snapshot
        )

    def test_otu_create_with_acronym_auto(self, precached_repo: Repo):
        otu = create_otu(
            precached_repo,
            132477,
            ["NC_013006"],
            "",
        )

        assert otu.acronym == "KLV"

    def test_otu_create_with_acronym_manual(self, precached_repo: Repo):
        otu = create_otu(
            precached_repo,
            1441799,
            ["NC_023881"],
            "FBNSV",
        )

        assert otu.acronym == "FBNSV"


class TestCreateOTUCommands:
    @pytest.mark.parametrize(
        "taxid, accessions",
        [(1278205, ["NC_020160"]), (345184, ["DQ178610", "DQ178611"])],
    )
    def test_ok(
        self,
        taxid: int,
        accessions: list[str],
        precached_repo: Repo,
        snapshot: SnapshotAssertion,
    ):
        run_create_otu_command(
            taxid=taxid,
            path=precached_repo.path,
            accessions=accessions,
        )

        otus = list(Repo(precached_repo.path).iter_otus())

        assert len(otus) == 1
        otu = otus[0]

        assert otu.model_dump() == snapshot(
            exclude=props("id", "isolates", "repr_isolate"),
        )

    def test_acronym(self, precached_repo: Repo):
        """Test if the --acronym option works as planned."""
        run_create_otu_command(
            taxid=345184,
            accessions=["DQ178610", "DQ178611"],
            path=precached_repo.path,
            acronym="CabLCJV",
        )

        otus = list(Repo(precached_repo.path).iter_otus())

        assert len(otus) == 1
        assert otus[0].acronym == "CabLCJV"


class TestAddIsolate:
    def test_duplicate_accessions(self, precached_repo: Repo):
        """Test that an error is raised when duplicate accessions are provided."""
        runner = CliRunner()
        result = runner.invoke(
            isolate_create,
            [
                "345184",
                "DQ178610",
                "DQ178610",
            ],
        )

        assert result.exit_code == 2
        assert "Duplicate accessions are not allowed." in result.output

    def test_genbank_ok(self, precached_repo: Repo):
        """Test that add_genbank_isolate() adds an isolate with a correctly parsed
        name.
        """
        isolate_1_accessions = ["DQ178610", "DQ178611"]
        isolate_2_accessions = ["DQ178613", "DQ178614"]

        otu = create_otu(precached_repo, 345184, isolate_1_accessions, acronym="")

        assert otu.accessions == set(isolate_1_accessions)

        isolate = add_genbank_isolate(precached_repo, otu, isolate_2_accessions)

        otu_after = precached_repo.get_otu_by_taxid(345184)

        assert otu_after.accessions == set(isolate_1_accessions).union(
            set(isolate_2_accessions),
        )

        isolate_after = otu_after.get_isolate(isolate.id)

        assert isolate_after.accessions == set(isolate_2_accessions)

        assert isolate_after.name == IsolateName(
            IsolateNameType.ISOLATE,
            "Douglas Castle",
        )

    def test_ignore_name_ok(self, precached_repo: Repo):
        """Test that add_unnamed_isolate() adds the correct isolate with its name set to None."""
        isolate_1_accessions = ["DQ178610", "DQ178611"]
        isolate_2_accessions = ["DQ178613", "DQ178614"]

        otu = create_otu(precached_repo, 345184, isolate_1_accessions, acronym="")

        assert otu.accessions == set(isolate_1_accessions)

        isolate = add_unnamed_isolate(precached_repo, otu, isolate_2_accessions)

        otu_after = precached_repo.get_otu_by_taxid(345184)

        assert otu_after.isolate_ids == {otu_after.repr_isolate, isolate.id}

        isolate_after = otu_after.get_isolate(isolate.id)

        assert isolate_after.name is None

        assert isolate_after.accessions == {"DQ178613", "DQ178614"}

    def test_add_and_name_ok(self, precached_repo: Repo):
        """Test that add_and_name_isolate() creates an isolate with the correct name."""
        isolate_1_accessions = ["DQ178610", "DQ178611"]
        isolate_2_accessions = ["DQ178613", "DQ178614"]

        otu = create_otu(precached_repo, 345184, isolate_1_accessions, acronym="")

        assert otu.accessions == set(isolate_1_accessions)

        isolate = add_and_name_isolate(
            precached_repo,
            otu,
            isolate_2_accessions,
            isolate_name=IsolateName(type=IsolateNameType.ISOLATE, value="dummy"),
        )

        otu_after = precached_repo.get_otu_by_taxid(345184)

        assert otu_after.isolate_ids == {otu_after.repr_isolate, isolate.id}

        isolate_after = otu_after.get_isolate(isolate.id)

        assert isolate_after.name == IsolateName(
            type=IsolateNameType.ISOLATE,
            value="dummy",
        )

        assert isolate_after.accessions == {"DQ178613", "DQ178614"}

    def test_conflict_fail(self, precached_repo: Repo):
        """Test that an isolate cannot be added to an OTU
        if both its name and its accessions are already contained.
        """
        taxid = 2164102
        accessions = ["MF062136", "MF062137", "MF062138"]

        otu = create_otu(
            precached_repo,
            taxid,
            accessions,
            "",
        )

        assert add_genbank_isolate(precached_repo, otu, accessions) is None


@pytest.mark.parametrize(
    "taxid, original_accessions, refseq_accessions",
    [
        (1169032, ["AB017503"], ["NC_003355"]),
        (345184, ["DQ178608", "DQ178609"], ["NC_038792", "NC_038793"]),
    ],
)
class TestReplaceIsolateSequences:
    def test_manual_replace_ok(
        self,
        empty_repo,
        taxid,
        original_accessions,
        refseq_accessions,
    ):
        """Test that a requested replacement occurs as expected."""
        create_otu(
            empty_repo,
            taxid,
            accessions=original_accessions,
            acronym="",
        )

        otu_before = empty_repo.get_otu_by_taxid(taxid)

        assert otu_before.accessions == set(original_accessions)

        isolate = list(otu_before.isolates)[0]

        update_isolate_from_accessions(
            empty_repo,
            otu_before,
            isolate.name,
            refseq_accessions,
        )

        otu_after = empty_repo.get_otu(otu_before.id)

        assert otu_after.accessions == set(refseq_accessions)

        assert otu_after.excluded_accessions == set(original_accessions)

    def test_raise_refseq_exception(
        self,
        empty_repo,
        taxid: int,
        original_accessions: list[str],
        refseq_accessions: list[str],
    ):
        """Test that attempting to add an isolate with RefSeq accessions
        raises RefSeqConflictError
        """
        create_otu(
            empty_repo,
            taxid,
            accessions=original_accessions,
            acronym="",
        )

        otu_before = empty_repo.get_otu_by_taxid(taxid)

        assert otu_before.accessions == set(original_accessions)

        add_genbank_isolate(empty_repo, otu_before, original_accessions)

        with pytest.raises(RefSeqConflictError):
            add_genbank_isolate(empty_repo, otu_before, refseq_accessions)

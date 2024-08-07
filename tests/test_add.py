import subprocess
from pathlib import Path

import pytest
from syrupy import SnapshotAssertion
from syrupy.filters import props

from ref_builder.otu import create_otu, update_otu, add_isolate
from ref_builder.repo import Repo


def run_create_otu_command(
    path: Path,
    taxid: int,
    accessions: list,
    acronym: str = "",
    autofill: bool = False,
):
    autofill_option = ["--autofill"] if autofill else []

    subprocess.run(
        ["ref-builder", "otu", "create"]
        + [str(taxid)]
        + accessions
        + ["--path", str(path)]
        + ["--acronym", acronym]
        + autofill_option,
        check=False,
    )


def run_update_otu_command(taxid: int, path: Path):
    subprocess.run(
        ["ref-builder", "otu", "update"] + [str(taxid)] + ["--path", str(path)],
        check=False,
    )


class TestCreateOTU:
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

        assert otu.dict() == snapshot(exclude=props("id", "isolates", "repr_isolate"))

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
        assert otu.schema is not None
        assert otu.repr_isolate is not None

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

        assert otu.dict() == snapshot(exclude=props("id", "isolates", "repr_isolate"))

    @pytest.mark.ncbi()
    def test_autofill_ok(
        self,
        precached_repo: Repo,
        snapshot: SnapshotAssertion,
    ):
        run_create_otu_command(
            taxid=345184,
            accessions=["DQ178610", "DQ178611"],
            path=precached_repo.path,
            autofill=True,
        )

        otus = list(Repo(precached_repo.path).iter_otus())

        assert len(otus) == 1
        otu = otus[0]

        assert otu.schema == snapshot
        assert {"DQ178610", "DQ178611"}.intersection(otu.accessions)

    def test_add_acronym_ok(self, precached_repo: Repo):
        """Test if the --acronym option works as planned."""
        run_create_otu_command(
            taxid=345184,
            accessions=["DQ178610", "DQ178611"],
            path=precached_repo.path,
            acronym="CabLCJV",
            autofill=True,
        )

        otus = list(Repo(precached_repo.path).iter_otus())

        assert len(otus) == 1
        otu = otus[0]

        assert otu.acronym == "CabLCJV"


class TestAddIsolate:
    def test_ok(self, precached_repo: Repo):
        isolate_1_accessions = ["DQ178610", "DQ178611"]
        isolate_2_accessions = ["DQ178613", "DQ178614"]

        otu = create_otu(precached_repo, 345184, isolate_1_accessions, acronym="")

        assert otu.accessions == set(isolate_1_accessions)

        isolate = add_isolate(precached_repo, otu, isolate_2_accessions)

        otu = precached_repo.get_otu_by_taxid(345184)

        assert otu.accessions == set(isolate_1_accessions).union(set(isolate_2_accessions))

        assert otu.get_isolate(isolate.id).accessions == set(isolate_2_accessions)


@pytest.mark.ncbi()
class TestUpdateOTU:
    def test_without_exclusions_ok(
        self,
        precached_repo: Repo,
        snapshot: SnapshotAssertion,
    ):
        otu = create_otu(
            precached_repo,
            2164102,
            ["MF062136", "MF062137", "MF062138"],
            "",
        )

        assert otu.accessions == {"MF062136", "MF062137", "MF062138"}

        update_otu(precached_repo, otu)

        otu = precached_repo.get_otu(otu.id)

        assert [otu.dict() for otu in precached_repo.iter_otus()] == snapshot(
            exclude=props("id", "isolates", "repr_isolate"),
        )

        assert otu.accessions == {
            "MF062136",
            "MF062137",
            "MF062138",
            "MF062130",
            "MF062131",
            "MF062132",
            "OQ420743",
            "OQ420744",
            "OQ420745",
            "MK936225",
            "MK936226",
            "MK936227",
        }

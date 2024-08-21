import subprocess
from pathlib import Path
from uuid import UUID

import pytest
from syrupy import SnapshotAssertion
from syrupy.filters import props

from ref_builder.otu.create import create_otu
from ref_builder.otu.update import (
    auto_update_otu,
    add_isolate,
    remove_isolate_from_otu,
    update_isolate_from_accessions,
)
from ref_builder.repo import Repo
from ref_builder.utils import IsolateName, IsolateNameType


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

        assert otu.schema.model_dump() == snapshot(exclude=props("segments"))

        for segment in otu.schema.segments:
            assert segment.model_dump() == snapshot(exclude=props("id"))

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
        otu_before = create_otu(
            precached_repo,
            2164102,
            ["MF062136", "MF062137", "MF062138"],
            "",
        )

        assert otu_before.accessions == {"MF062136", "MF062137", "MF062138"}

        auto_update_otu(precached_repo, otu_before)

        otu_after = precached_repo.get_otu(otu_before.id)

        assert otu_after.id == otu_before.id

        assert otu_after.isolate_ids.issuperset(otu_before.isolate_ids)

        assert {isolate.name: isolate.accessions for isolate in otu_after.isolates} == snapshot()

        assert otu_after.excluded_accessions == {"MF062125", "MF062126", "MF062127"}

        assert otu_after.accessions == {
            "MF062136",
            "MF062137",
            "MF062138",
            "MF062130",
            "MF062131",
            "MF062132",
            "MK936225",
            "MK936226",
            "MK936227",
            "OQ420743",
            "OQ420744",
            "OQ420745",
            "NC_055390",
            "NC_055391",
            "NC_055392",
        }


class TestRemoveIsolate:
    @pytest.mark.parametrize(
        "taxid, isolate_name",
        [
            (1169032, IsolateName(type=IsolateNameType.ISOLATE, value="WMoV-6.3"))
        ]
    )
    def test_ok(self, scratch_repo, taxid: int, isolate_name: IsolateName):
        """Test that a given isolate can be removed from the OTU."""
        otu_before = scratch_repo.get_otu_by_taxid(taxid)

        isolate_id = otu_before.get_isolate_id_by_name(isolate_name)

        assert type(isolate_id) is UUID

        remove_isolate_from_otu(scratch_repo, otu_before, isolate_id)

        otu_after = scratch_repo.get_otu_by_taxid(taxid)

        assert otu_after.get_isolate(isolate_id) is None

        assert otu_before.get_isolate(isolate_id).accessions not in otu_after.accessions

        assert len(otu_after.isolate_ids) == len(otu_before.isolate_ids) - 1


class TestReplaceIsolateSequences:
    @pytest.mark.parametrize(
        "taxid, original_accessions, refseq_accessions",
        [
            (1169032, ["AB017503"], ["NC_003355"]),
            (345184, ["DQ178608", "DQ178609"], ["NC_038792", "NC_038793"])
        ]
    )
    def test_manual_replace_ok(
        self,
            empty_repo,
            taxid: int,
            original_accessions: list[str],
            refseq_accessions: list[str]
    ):
        """Test that a given isolate can receive a new set of accessions
        and update its contents accordingly.
        The OTU should also update its excluded list."""
        otu_before = create_otu(empty_repo, taxid, accessions=original_accessions, acronym="")

        otu_before = empty_repo.get_otu(otu_before.id)

        assert otu_before.accessions == set(original_accessions)

        isolate_before = list(otu_before.isolates)[0]

        isolate_after = update_isolate_from_accessions(
            empty_repo, otu_before, isolate_before.name, refseq_accessions
        )

        otu_after = empty_repo.get_otu(otu_before.id)

        assert otu_after.get_isolate(isolate_after.id).accessions == set(refseq_accessions)

        assert otu_after.excluded_accessions == set(original_accessions)

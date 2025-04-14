import pytest
from click.testing import CliRunner
from syrupy import SnapshotAssertion
from syrupy.filters import props


from ref_builder.console import console, print_otu
from ref_builder.cli.otu import otu as otu_command_group
from ref_builder.otu.create import create_otu_with_taxid, create_otu_without_taxid
from ref_builder.otu.isolate import (
    add_and_name_isolate,
    add_genbank_isolate,
    add_unnamed_isolate,
)
from ref_builder.otu.update import (
    update_isolate_from_accessions,
)
from ref_builder.otu.utils import RefSeqConflictError
from ref_builder.repo import Repo
from ref_builder.utils import IsolateName, IsolateNameType


class TestCreateOTU:
    @pytest.mark.parametrize(
        ("taxid", "accessions"),
        [(1278205, ["NC_020160"]), (345184, ["DQ178610", "DQ178611"])],
    )
    def test_ok(
        self,
        taxid: int,
        accessions: list[str],
        precached_repo: Repo,
    ):
        """Test that an OTU can be created."""
        assert list(precached_repo.iter_otus()) == []

        with precached_repo.lock():
            otu_ = create_otu_with_taxid(
                precached_repo,
                taxid,
                accessions,
                "",
            )

        assert list(precached_repo.iter_otus()) == [otu_]
        assert otu_.plan.monopartite == (len(accessions) == 1)
        assert len(otu_.plan.segments) == len(accessions)
        assert otu_.accessions == set(accessions)
        assert otu_.representative_isolate == otu_.isolates[0].id

    def test_empty_repo(
        self,
        precached_repo: Repo,
        snapshot: SnapshotAssertion,
    ):
        """Test that an OTU can be created in an empty repository."""
        with precached_repo.lock():
            otu_ = create_otu_with_taxid(
                precached_repo,
                345184,
                ["DQ178610", "DQ178611"],
                "",
            )

        assert otu_.model_dump() == snapshot(
            exclude=props("id", "isolates", "representative_isolate"),
        )

        # Ensure only one OTU is present in the repository, and it matches the return
        # value of the creation function.
        assert list(precached_repo.iter_otus()) == [otu_]

    def test_duplicate_taxid(self, precached_repo: Repo):
        """Test that an OTU with the same taxid cannot be created."""
        accessions = ["DQ178610", "DQ178611"]
        taxid = 345184

        with precached_repo.lock():
            create_otu_with_taxid(
                precached_repo,
                taxid,
                accessions,
                "",
            )

        with (
            pytest.raises(
                ValueError,
                match="Taxonomy ID 345184 has already been added to this reference.",
            ),
            precached_repo.lock(),
        ):
            create_otu_with_taxid(
                precached_repo,
                taxid,
                accessions,
                "",
            )

    def test_refseq_autoexclude(self, precached_repo: Repo):
        """Test that the superceded accessions included in RefSeq metadata are
        automatically added to the OTU's excluded accessions list.
        """
        with precached_repo.lock():
            otu_ = create_otu_with_taxid(
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

        assert otu_.excluded_accessions == {
            "EF546808",
            "EF546809",
            "EF546810",
            "EF546811",
            "EF546812",
            "EF546813",
        }

    def test_acronym_auto(self, precached_repo: Repo) -> None:
        """Test that the auto-generated acronym is correct."""
        with precached_repo.lock():
            otu_ = create_otu_with_taxid(
                precached_repo,
                132477,
                ["NC_013006"],
                "",
            )

        assert otu_.acronym == "KLV"

    def test_acronym_manual(self, precached_repo: Repo) -> None:
        """Test the acronym can be set manually."""
        with precached_repo.lock():
            otu_ = create_otu_with_taxid(
                precached_repo,
                1441799,
                ["NC_023881"],
                "FBNSV",
            )

        assert otu_.acronym == "FBNSV"

    @pytest.mark.parametrize(
        ("accessions", "expected_taxid"),
        [
            (["DQ178610", "DQ178611"], 345184),
            (["NC_043170"], 3240630),
        ],
    )
    def test_create_without_taxid_ok(self, accessions, expected_taxid, precached_repo):
        with precached_repo.lock():
            made_otu = create_otu_without_taxid(
                precached_repo, accessions=accessions, acronym=""
            )

        assert made_otu.taxid == expected_taxid

        otu = precached_repo.get_otu_by_taxid(expected_taxid)

        assert otu is not None


class TestAddIsolate:
    def test_multipartite(self, precached_repo: Repo):
        with precached_repo.lock():
            otu_before = create_otu_with_taxid(
                precached_repo,
                2164102,
                ["MF062136", "MF062137", "MF062138"],
                acronym="",
            )

        assert otu_before.accessions == {"MF062136", "MF062137", "MF062138"}
        assert len(otu_before.isolate_ids) == 1

        with precached_repo.lock():
            isolate = add_genbank_isolate(
                precached_repo, otu_before, ["MF062125", "MF062126", "MF062127"]
            )

        otu_after = precached_repo.get_otu(otu_before.id)

        assert otu_after.accessions == {
            "MF062125",
            "MF062126",
            "MF062127",
            "MF062136",
            "MF062137",
            "MF062138",
        }

        assert otu_after.get_isolate(isolate.id).accessions == {
            "MF062125",
            "MF062126",
            "MF062127",
        }

    def test_genbank(self, precached_repo: Repo):
        """Test that add_genbank_isolate() adds an isolate with a correctly parsed
        name.
        """
        isolate_1_accessions = ["DQ178610", "DQ178611"]
        isolate_2_accessions = ["DQ178613", "DQ178614"]

        with precached_repo.lock():
            otu_before = create_otu_with_taxid(
                precached_repo, 345184, isolate_1_accessions, acronym=""
            )

            assert otu_before.accessions == set(isolate_1_accessions)

            isolate = add_genbank_isolate(
                precached_repo, otu_before, isolate_2_accessions
            )

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

    def test_ignore_name(self, precached_repo: Repo):
        """Test that add_unnamed_isolate() adds the isolate with a ``None`` name."""
        isolate_1_accessions = ["DQ178610", "DQ178611"]
        isolate_2_accessions = ["DQ178613", "DQ178614"]

        with precached_repo.lock():
            otu = create_otu_with_taxid(
                precached_repo, 345184, isolate_1_accessions, acronym=""
            )

        assert otu.accessions == set(isolate_1_accessions)

        with precached_repo.lock():
            isolate = add_unnamed_isolate(precached_repo, otu, isolate_2_accessions)

        otu_after = precached_repo.get_otu_by_taxid(345184)

        assert otu_after.isolate_ids == {otu_after.representative_isolate, isolate.id}

        isolate_after = otu_after.get_isolate(isolate.id)

        assert isolate_after.name is None
        assert isolate_after.accessions == {"DQ178613", "DQ178614"}

    def test_add_and_name_isolate(self, precached_repo: Repo):
        """Test that add_and_name_isolate() creates an isolate with the correct name."""
        isolate_1_accessions = ["DQ178610", "DQ178611"]
        isolate_2_accessions = ["DQ178613", "DQ178614"]

        with precached_repo.lock():
            otu = create_otu_with_taxid(
                precached_repo, 345184, isolate_1_accessions, acronym=""
            )

        assert otu.accessions == set(isolate_1_accessions)

        with precached_repo.lock():
            isolate = add_and_name_isolate(
                precached_repo,
                otu,
                isolate_2_accessions,
                isolate_name=IsolateName(type=IsolateNameType.ISOLATE, value="dummy"),
            )

        otu_after = precached_repo.get_otu_by_taxid(345184)

        assert otu_after.isolate_ids == {otu_after.representative_isolate, isolate.id}

        isolate_after = otu_after.get_isolate(isolate.id)

        assert isolate_after.name == IsolateName(
            type=IsolateNameType.ISOLATE,
            value="dummy",
        )

        assert isolate_after.accessions == {"DQ178613", "DQ178614"}

    def test_blocked_fail(self, precached_repo: Repo):
        """Test that an isolate cannot be added to an OTU if both its name and its
        accessions are already contained.
        """
        accessions = ["MF062136", "MF062137", "MF062138"]

        with precached_repo.lock():
            otu = create_otu_with_taxid(
                precached_repo,
                2164102,
                accessions,
                "",
            )

            assert add_genbank_isolate(precached_repo, otu, accessions) is None

    def test_name_conflict_fail(self, precached_repo: Repo):
        """Test that an isolate cannot be added to an OTU if its name is already contained."""
        with precached_repo.lock():
            otu_init = create_otu_with_taxid(
                precached_repo,
                2164102,
                ["MF062136", "MF062137", "MF062138"],
                "",
            )

        last_event_id_before_test = precached_repo.last_id

        with precached_repo.lock():
            bad_isolate = add_and_name_isolate(
                precached_repo,
                otu_init,
                ["NC_055390", "NC_055391", "NC_055392"],
                IsolateName(type=IsolateNameType.ISOLATE, value="4342-5"),
            )

        otu_after = precached_repo.get_otu(otu_init.id)

        assert bad_isolate.id not in otu_after.isolate_ids

        assert precached_repo.last_id == last_event_id_before_test

    def test_plan_mismatch_fail(self, scratch_repo: Repo):
        """Test that an isolate cannot be added to an OTU if it does not match the OTU plan."""
        otu_before = scratch_repo.get_otu_by_taxid(345184)

        with scratch_repo.lock():
            assert (
                add_genbank_isolate(
                    scratch_repo,
                    otu_before,
                    accessions=["AB017503"],
                )
                is None
            )

        otu_after = scratch_repo.get_otu_by_taxid(345184)

        assert "AB017503" not in otu_after.accessions


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
        empty_repo: Repo,
        taxid: int,
        original_accessions: list[str],
        refseq_accessions: list[str],
    ):
        """Test that a requested replacement occurs as expected."""
        with empty_repo.lock():
            create_otu_with_taxid(
                empty_repo,
                taxid,
                accessions=original_accessions,
                acronym="",
            )

        otu_before = empty_repo.get_otu_by_taxid(taxid)

        assert otu_before.accessions == set(original_accessions)

        isolate = next(iter(otu_before.isolates))

        with empty_repo.lock():
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
        empty_repo: Repo,
        taxid: int,
        original_accessions: list[str],
        refseq_accessions: list[str],
    ):
        """Test that attempting to add an isolate with RefSeq accessions
        raises RefSeqConflictError
        """
        with empty_repo.lock():
            create_otu_with_taxid(
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

import subprocess

from ref_builder.repo import Repo
from ref_builder.otu.create import create_otu
from ref_builder.otu.update import (
    add_isolate,
    promote_otu_accessions,
    set_representative_isolate,
)


def test_update_representative_isolate(scratch_repo: Repo):
    taxid = 345184

    otu_before = scratch_repo.get_otu_by_taxid(taxid)

    repr_isolate_after = None
    for isolate_id in otu_before.isolate_ids:
        if isolate_id != otu_before.repr_isolate:
            repr_isolate_after = isolate_id
            break

    set_representative_isolate(scratch_repo, otu_before, repr_isolate_after)

    otu_after = scratch_repo.get_otu_by_taxid(taxid)

    assert otu_after.repr_isolate != otu_before.repr_isolate

    assert otu_after.repr_isolate == repr_isolate_after


class TestUpdateRepresentativeIsolateCommand:
    def test_isolate_id_ok(self, scratch_repo: Repo):
        taxid = 345184

        otu_before = scratch_repo.get_otu_by_taxid(taxid)

        repr_isolate_after = None

        for isolate_id in otu_before.isolate_ids:
            if isolate_id != otu_before.repr_isolate:
                repr_isolate_after = isolate_id
                break

        subprocess.run(
            [
                "ref-builder",
                "otu",
                "update",
            ]
            + ["--path", str(scratch_repo.path)]
            + [str(taxid)]
            + ["default"] +
            [str(repr_isolate_after)],
            check=False,
        )

        scratch_repo = Repo(scratch_repo.path)

        otu_after = scratch_repo.get_otu_by_taxid(taxid)

        assert otu_after.repr_isolate != otu_before.repr_isolate

        assert otu_after.repr_isolate == repr_isolate_after

    def test_isolate_name_ok(self, scratch_repo: Repo):
        taxid = 1169032

        otu_before = scratch_repo.get_otu_by_taxid(taxid)

        repr_isolate_name_before = otu_before.get_isolate(otu_before.repr_isolate).name

        repr_isolate_name_after = None

        for isolate_id in otu_before.isolate_ids:
            if isolate_id != otu_before.repr_isolate:
                repr_isolate_name_after = otu_before.get_isolate(isolate_id).name
                break

        subprocess.run(
            [
                "ref-builder",
                "otu",
                "update",
            ]
            + ["--path", str(scratch_repo.path)]
            + [str(taxid)]
            + ["default"] +
            [str(repr_isolate_name_after)],
            check=False,
        )

        scratch_repo = Repo(scratch_repo.path)

        otu_after = scratch_repo.get_otu_by_taxid(taxid)

        assert otu_after.repr_isolate != otu_before.repr_isolate

        assert otu_after.repr_isolate == otu_before.get_isolate_id_by_name(repr_isolate_name_after)


class TestPromoteAccessions:
    def test_ok(self, empty_repo: Repo):
        """Test that RefSeq accessions can be promoted automatically."""
        otu = create_otu(
            empty_repo,
            2164102,
            ["MF062136", "MF062137", "MF062138"],
            acronym=""
        )
        isolate = add_isolate(empty_repo, otu, ["MF062125", "MF062126", "MF062127"])

        otu_before = empty_repo.get_otu(otu.id)

        assert otu_before.accessions == {
            "MF062125", "MF062126", "MF062127",
            "MF062136", "MF062137", "MF062138"
        }

        assert otu_before.get_isolate(isolate.id).accessions == {"MF062125", "MF062126", "MF062127"}

        promoted_accessions = promote_otu_accessions(empty_repo, otu_before)

        assert promoted_accessions == {"NC_055390", "NC_055391", "NC_055392"}

        otu_after = empty_repo.get_otu(otu.id)

        assert otu_after.isolate_ids == otu_before.isolate_ids

        assert otu_after.accessions == {
            "NC_055390", "NC_055391", "NC_055392",
            "MF062136", "MF062137", "MF062138",
        }

        assert otu_after.get_isolate(isolate.id).accessions == {"NC_055390", "NC_055391", "NC_055392"}

        assert otu_after.excluded_accessions == {"MF062125", "MF062126", "MF062127"}

    def test_command_ok(self, empty_repo: Repo):
        otu = create_otu(
            empty_repo,
            2164102,
            ["MF062125", "MF062126", "MF062127"],
            acronym=""
        )

        otu_before = empty_repo.get_otu(otu.id)

        assert otu_before.accessions == {"MF062125", "MF062126", "MF062127"}

        subprocess.run(
            [
                "ref-builder",
                "otu",
                "update",
            ]
            + ["--path", str(empty_repo.path)]
            + [str(2164102)]
            + ["promote"],
            check=False,
        )

        repo_after = Repo(empty_repo.path)

        otu_after = repo_after.get_otu(otu.id)

        assert otu_after.repr_isolate == otu_before.repr_isolate

        assert otu_after.accessions == {"NC_055390", "NC_055391", "NC_055392",}

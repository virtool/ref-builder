import subprocess

from ref_builder.repo import Repo
from ref_builder.otu.update import set_representative_isolate


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

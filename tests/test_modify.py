import subprocess
from uuid import uuid4

import pytest
from syrupy import SnapshotAssertion
from syrupy.filters import props

from ref_builder.repo import Repo
from ref_builder.otu.create import create_otu
from ref_builder.otu.update import (
    add_genbank_isolate,
    add_segments_to_plan,
    promote_otu_accessions,
    replace_isolate_plan,
    set_isolate_plan,
    set_representative_isolate,
)
from ref_builder.plan import MultipartitePlan, SegmentName, SegmentRule, SegmentPlan


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


class TestSetIsolatePlan:
    def test_set_isolate_plan(self, scratch_repo: Repo):
        otu_before = scratch_repo.get_otu_by_taxid(223262)

        original_plan = otu_before.plan

        assert type(original_plan) is MultipartitePlan

        new_isolate_plan = MultipartitePlan(id=uuid4(), segments=original_plan.segments)

        new_isolate_plan.segments.append(
            SegmentPlan(
                id=uuid4(),
                length=2000,
                name=SegmentName(prefix="DNA", key="C"),
                required=SegmentRule.RECOMMENDED,
            ),
        )

        new_isolate_plan.segments.append(
            SegmentPlan(
                id=uuid4(),
                length=1000,
                name=SegmentName(prefix="DNA", key="Z"),
                required=SegmentRule.OPTIONAL,
            ),
        )

        set_isolate_plan(scratch_repo, otu_before, new_isolate_plan)

        assert type(new_isolate_plan) is MultipartitePlan

        otu_after = scratch_repo.get_otu(otu_before.id)

        assert len(otu_after.plan.segments) == len(otu_before.plan.segments) + 2

        assert otu_after.plan == new_isolate_plan

    def test_replace_isolate_plan(self, scratch_repo: Repo):
        otu_before = scratch_repo.get_otu_by_taxid(345184)

        assert otu_before

        new_plan = replace_isolate_plan(
            scratch_repo, otu_before, ["NC_038792", "NC_038793"]
        )

        otu_after = scratch_repo.get_otu(otu_before.id)

        assert otu_after.plan == new_plan

        assert otu_after.plan != otu_before.plan

    def test_add_segments_to_plan(
        self,
        precached_repo: Repo,
        snapshot: SnapshotAssertion,
    ):
        otu_before = create_otu(
            precached_repo,
            2164102,
            ["MF062136", "MF062137"],
            acronym="",
        )

        original_plan = otu_before.plan

        assert type(original_plan) is MultipartitePlan

        expanded_plan = add_segments_to_plan(
            precached_repo,
            otu_before,
            rule=SegmentRule.OPTIONAL,
            accessions=["MF062138"],
        )

        assert len(expanded_plan.segments) == len(original_plan.segments) + 1

        otu_after = precached_repo.get_otu(otu_before.id)

        assert otu_after.plan != otu_before.plan

        assert otu_after.plan.model_dump() == snapshot(exclude=props("id"))


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
            + ["default"]
            + [str(repr_isolate_after)],
            check=False,
        )

        scratch_repo = Repo(scratch_repo.path)

        otu_after = scratch_repo.get_otu_by_taxid(taxid)

        assert otu_after.repr_isolate != otu_before.repr_isolate

        assert otu_after.repr_isolate == repr_isolate_after

    def test_isolate_name_ok(self, scratch_repo: Repo):
        taxid = 1169032

        otu_before = scratch_repo.get_otu_by_taxid(taxid)

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
            + ["default"]
            + [str(repr_isolate_name_after)],
            check=False,
        )

        scratch_repo = Repo(scratch_repo.path)

        otu_after = scratch_repo.get_otu_by_taxid(taxid)

        assert otu_after.repr_isolate != otu_before.repr_isolate

        assert otu_after.repr_isolate == otu_before.get_isolate_id_by_name(
            repr_isolate_name_after
        )


@pytest.mark.ncbi()
class TestPromoteAccessions:
    def test_ok(self, empty_repo: Repo):
        """Test that RefSeq accessions can be promoted automatically."""
        otu = create_otu(
            empty_repo, 2164102, ["MF062136", "MF062137", "MF062138"], acronym=""
        )
        isolate = add_genbank_isolate(
            empty_repo, otu, ["MF062125", "MF062126", "MF062127"]
        )

        otu_before = empty_repo.get_otu(otu.id)

        assert otu_before.accessions == {
            "MF062125",
            "MF062126",
            "MF062127",
            "MF062136",
            "MF062137",
            "MF062138",
        }

        assert otu_before.get_isolate(isolate.id).accessions == {
            "MF062125",
            "MF062126",
            "MF062127",
        }

        promoted_accessions = promote_otu_accessions(empty_repo, otu_before)

        assert promoted_accessions == {"NC_055390", "NC_055391", "NC_055392"}

        otu_after = empty_repo.get_otu(otu.id)

        assert otu_after.isolate_ids == otu_before.isolate_ids

        assert otu_after.accessions == {
            "NC_055390",
            "NC_055391",
            "NC_055392",
            "MF062136",
            "MF062137",
            "MF062138",
        }

        assert otu_after.get_isolate(isolate.id).accessions == {
            "NC_055390",
            "NC_055391",
            "NC_055392",
        }

        assert otu_after.excluded_accessions == {"MF062125", "MF062126", "MF062127"}

    def test_command_ok(self, empty_repo: Repo):
        otu = create_otu(
            empty_repo, 2164102, ["MF062125", "MF062126", "MF062127"], acronym=""
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

        assert otu_after.accessions == {"NC_055390", "NC_055391", "NC_055392"}

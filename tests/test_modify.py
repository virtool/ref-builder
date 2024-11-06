import subprocess
from uuid import UUID

import pytest
from syrupy import SnapshotAssertion
from syrupy.filters import props

from ref_builder.repo import Repo
from ref_builder.resources import RepoSequence
from ref_builder.otu.create import create_otu
from ref_builder.otu.modify import (
    add_segments_to_plan,
    resize_monopartite_plan,
    set_isolate_plan,
    set_plan_length_tolerances,
)
from ref_builder.otu.update import (
    add_genbank_isolate,
    promote_otu_accessions,
    replace_sequence_in_otu,
    set_representative_isolate,
    delete_isolate_from_otu,
)
from ref_builder.utils import IsolateName, IsolateNameType
from ref_builder.plan import (
    MonopartitePlan,
    MultipartitePlan,
    SegmentName,
    SegmentRule,
    Segment,
)


def test_update_representative_isolate(scratch_repo: Repo):
    """Test representative isolate replacement."""
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

        new_isolate_plan = MultipartitePlan.new(segments=original_plan.segments)

        new_isolate_plan.segments.append(
            Segment.new(
                length=2000,
                length_tolerance=scratch_repo.settings.default_segment_length_tolerance,
                name=SegmentName(prefix="DNA", key="C"),
                required=SegmentRule.RECOMMENDED,
            ),
        )

        new_isolate_plan.segments.append(
            Segment.new(
                length=1000,
                length_tolerance=scratch_repo.settings.default_segment_length_tolerance,
                name=SegmentName(prefix="DNA", key="Z"),
                required=SegmentRule.OPTIONAL,
            ),
        )

        set_isolate_plan(scratch_repo, otu_before, new_isolate_plan)

        assert type(new_isolate_plan) is MultipartitePlan

        otu_after = scratch_repo.get_otu(otu_before.id)

        assert len(otu_after.plan.segments) == len(otu_before.plan.segments) + 2

        assert otu_after.plan == new_isolate_plan

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

    def test_resize_monopartite_plan(
        self,
        precached_repo: Repo,
        snapshot: SnapshotAssertion,
    ):
        otu_before = create_otu(
            precached_repo,
            2164102,
            ["MF062136"],
            acronym="",
        )

        assert type(otu_before.plan) is MonopartitePlan

        resize_monopartite_plan(
            precached_repo,
            otu_before,
            name=SegmentName(prefix="RNA", key="L"),
            rule=SegmentRule.RECOMMENDED,
            accessions=["MF062137", "MF062138"],
        )

        otu_after = precached_repo.get_otu(otu_before.id)

        assert type(otu_after.plan) is MultipartitePlan

        assert otu_after.plan.required_segments[0].length == otu_before.plan.length

        assert otu_after.plan.model_dump() == snapshot(exclude=props("id"))

    def test_extend_plan_monopartite_fail(
        self,
        scratch_repo: Repo,
    ):
        """Test that add_segments_to_plan() fails out
        when the original plan is a MonopartitePlan.
        """
        otu_before = scratch_repo.get_otu_by_taxid(96892)

        assert type(otu_before.plan) is MonopartitePlan

        with pytest.raises(ValueError):
            add_segments_to_plan(
                scratch_repo,
                otu_before,
                rule=SegmentRule.OPTIONAL,
                accessions=["NC_010620"],
            )

    @pytest.mark.parametrize("tolerance", [0.05, 0.5, 1.0])
    def test_set_length_tolerances_ok(self, scratch_repo: Repo, tolerance: float):
        """Check that plan length tolerances can be modified by function."""
        otu_before = scratch_repo.get_otu_by_taxid(96892)

        assert (
            otu_before.plan.length_tolerance
            == scratch_repo.settings.default_segment_length_tolerance
        )

        set_plan_length_tolerances(scratch_repo, otu_before, tolerance)

        otu_after = scratch_repo.get_otu(otu_before.id)

        assert otu_after.plan.length_tolerance == tolerance

    @pytest.mark.parametrize("bad_tolerance", [-1.0, 1.1, 100.0])
    def test_set_length_tolerances_fail(self, scratch_repo: Repo, bad_tolerance: float):
        """Check that plan length tolerances cannot be set to an invalid float value."""
        otu_before = scratch_repo.get_otu_by_taxid(96892)

        assert (
            otu_before.plan.length_tolerance
            == scratch_repo.settings.default_segment_length_tolerance
        )

        set_plan_length_tolerances(scratch_repo, otu_before, bad_tolerance)

        otu_after = scratch_repo.get_otu(otu_before.id)

        assert otu_after.plan.length_tolerance == otu_before.plan.length_tolerance


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


class TestDeleteIsolate:
    def test_ok(self, scratch_repo):
        """Test that a given isolate can be deleted from the OTU."""
        taxid = 1169032

        otu_before = scratch_repo.get_otu_by_taxid(taxid)

        isolate_id = otu_before.get_isolate_id_by_name(
            IsolateName(type=IsolateNameType.ISOLATE, value="WMoV-6.3"),
        )

        assert type(isolate_id) is UUID

        delete_isolate_from_otu(scratch_repo, otu_before, isolate_id)

        otu_after = scratch_repo.get_otu_by_taxid(taxid)

        assert isolate_id not in otu_after.isolate_ids

        assert otu_after.get_isolate(isolate_id) is None

        assert otu_before.get_isolate(isolate_id).accessions not in otu_after.accessions

        assert len(otu_after.isolate_ids) == len(otu_before.isolate_ids) - 1


class TestReplaceSequence:
    def test_ok(self, precached_repo):
        """Test sequence replacement and deletion."""
        otu_before = create_otu(
            precached_repo,
            1169032,
            ["MK431779"],
            acronym="",
        )

        isolate_id, old_sequence_id = (
            otu_before.get_sequence_id_hierarchy_from_accession(
                "MK431779",
            )
        )

        assert type(old_sequence_id) is UUID

        sequence = replace_sequence_in_otu(
            repo=precached_repo,
            otu=otu_before,
            new_accession="NC_003355",
            replaced_accession="MK431779",
        )

        assert type(sequence) is RepoSequence

        otu_after = precached_repo.get_otu_by_taxid(1169032)

        assert (
            otu_after.accessions
            == otu_after.get_isolate(isolate_id).accessions
            == {"NC_003355"}
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

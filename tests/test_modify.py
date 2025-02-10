import subprocess
from uuid import UUID, uuid4

import pytest
from syrupy import SnapshotAssertion
from syrupy.filters import props

from ref_builder.otu.create import create_otu
from ref_builder.otu.isolate import add_genbank_isolate
from ref_builder.otu.modify import (
    add_segments_to_plan,
    allow_accessions_into_otu,
    delete_isolate_from_otu,
    exclude_accessions_from_otu,
    rename_plan_segment,
    replace_sequence_in_otu,
    set_plan,
    set_plan_length_tolerances,
    set_representative_isolate,
)
from ref_builder.otu.update import (
    promote_otu_accessions,
)
from ref_builder.plan import (
    Plan,
    Segment,
    SegmentName,
    SegmentRule,
)
from ref_builder.repo import Repo
from ref_builder.resources import RepoSequence
from ref_builder.utils import IsolateName, IsolateNameType


def test_exclude_accessions(scratch_repo: Repo):
    """Test accession exclusion."""
    taxid = 345184

    otu_before = scratch_repo.get_otu_by_taxid(taxid)

    assert not otu_before.excluded_accessions

    with scratch_repo.lock():
        exclude_accessions_from_otu(
            scratch_repo, otu_before, accessions=["DQ178608", "DQ178609"]
        )

    otu_after = scratch_repo.get_otu_by_taxid(taxid)

    assert otu_after.excluded_accessions == {"DQ178608", "DQ178609"}


def test_allow_accessions(scratch_repo: Repo):
    taxid = 345184

    with scratch_repo.lock():
        exclude_accessions_from_otu(
            scratch_repo,
            otu=scratch_repo.get_otu_by_taxid(taxid),
            accessions={"DQ178608", "DQ178609"},
        )

    otu_before = scratch_repo.get_otu_by_taxid(taxid)

    assert otu_before.excluded_accessions == {"DQ178608", "DQ178609"}

    with scratch_repo.lock():
        allow_accessions_into_otu(
            scratch_repo,
            otu=otu_before,
            accessions={"DQ178608"},
        )

    otu_after = scratch_repo.get_otu_by_taxid(taxid)

    assert otu_after.excluded_accessions == {"DQ178609"}


def test_update_representative_isolate(scratch_repo: Repo):
    """Test representative isolate replacement."""
    taxid = 345184

    otu_before = scratch_repo.get_otu_by_taxid(taxid)

    representative_isolate_after = None

    for isolate_id in otu_before.isolate_ids:
        if isolate_id != otu_before.representative_isolate:
            representative_isolate_after = isolate_id
            break

    with scratch_repo.lock():
        set_representative_isolate(
            scratch_repo, otu_before, representative_isolate_after
        )

    otu_after = scratch_repo.get_otu_by_taxid(taxid)

    assert otu_after.representative_isolate != otu_before.representative_isolate

    assert otu_after.representative_isolate == representative_isolate_after


class TestSetPlan:
    """Test functions that make changes to an OTU plan."""

    def test_ok(self, scratch_repo: Repo):
        """Test that an OTU's plan can be replaced."""
        otu_before = scratch_repo.get_otu_by_taxid(223262)

        original_plan = otu_before.plan

        assert type(original_plan) is Plan

        new_plan = Plan.new(segments=original_plan.segments)

        new_plan.segments.append(
            Segment.new(
                length=2000,
                length_tolerance=scratch_repo.settings.default_segment_length_tolerance,
                name=SegmentName(prefix="DNA", key="C"),
                required=SegmentRule.RECOMMENDED,
            ),
        )

        new_plan.segments.append(
            Segment.new(
                length=1000,
                length_tolerance=scratch_repo.settings.default_segment_length_tolerance,
                name=SegmentName(prefix="DNA", key="Z"),
                required=SegmentRule.OPTIONAL,
            ),
        )

        set_plan(scratch_repo, otu_before, new_plan)

        assert type(new_plan) is Plan

        otu_after = scratch_repo.get_otu(otu_before.id)

        assert len(otu_after.plan.segments) == len(otu_before.plan.segments) + 2

        assert otu_after.plan == new_plan

    def test_rename_segment_ok(self, scratch_repo: Repo):
        """Test that a given plan segment can be renamed."""
        otu_before = scratch_repo.get_otu_by_taxid(223262)

        first_segment_id = otu_before.plan.segments[0].id

        new_name = SegmentName(prefix="RNA", key="TestName")

        assert otu_before.plan.get_segment_by_id(first_segment_id).name != new_name

        rename_plan_segment(
            scratch_repo,
            otu_before,
            segment_id=first_segment_id,
            segment_name=SegmentName(prefix="RNA", key="TestName"),
        )

        otu_after = scratch_repo.get_otu_by_taxid(223262)

        assert otu_after.plan.get_segment_by_id(first_segment_id).name == new_name

    def test_rename_segment_fail(self, scratch_repo: Repo):
        """Test that an attempt to rename a nonexistent segment does not change the OTU."""
        otu_before = scratch_repo.get_otu_by_taxid(223262)

        first_segment = otu_before.plan.segments[0]

        new_name = SegmentName(prefix="RNA", key="TestName")

        assert otu_before.plan.get_segment_by_id(first_segment.id).name != new_name

        assert (
            rename_plan_segment(
                scratch_repo,
                otu_before,
                segment_id=uuid4(),
                segment_name=SegmentName(prefix="RNA", key="TestName"),
            )
            is None
        )

        otu_after = scratch_repo.get_otu_by_taxid(223262)

        assert otu_after.plan.model_dump_json() == otu_before.plan.model_dump_json()

        assert (
            otu_after.plan.get_segment_by_id(first_segment.id).name
            == first_segment.name
        )

    @pytest.mark.parametrize("accessions", [["MF062136", "MF062137"], ["MF062136"]])
    def test_add_segments_to_plan_ok(
        self,
        precached_repo: Repo,
        accessions: list[str],
        snapshot: SnapshotAssertion,
    ):
        """Test the addition of segments to an OTU plan."""
        otu_before = create_otu(
            precached_repo,
            2164102,
            accessions,
            acronym="",
        )

        original_plan = otu_before.plan

        assert type(original_plan) is Plan

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

    def test_add_segments_to_plan_fail(
        self,
        scratch_repo: Repo,
    ):
        """Test that segments cannot be added to a monopartite plan with
        a preexisting unnamed segment.
        """
        otu_before = scratch_repo.get_otu_by_taxid(96892)

        assert otu_before.plan.monopartite

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
            otu_before.plan.segments[0].length_tolerance
            == scratch_repo.settings.default_segment_length_tolerance
        )

        set_plan_length_tolerances(scratch_repo, otu_before, tolerance)

        otu_after = scratch_repo.get_otu(otu_before.id)

        assert otu_after.plan.segments[0].length_tolerance == tolerance

    @pytest.mark.parametrize("bad_tolerance", [-1.0, 1.1, 100.0])
    def test_set_length_tolerances_fail(self, scratch_repo: Repo, bad_tolerance: float):
        """Check that plan length tolerances cannot be set to an invalid float value."""
        otu_before = scratch_repo.get_otu_by_taxid(96892)

        assert (
            otu_before.plan.segments[0].length_tolerance
            == scratch_repo.settings.default_segment_length_tolerance
        )

        set_plan_length_tolerances(scratch_repo, otu_before, bad_tolerance)

        otu_after = scratch_repo.get_otu(otu_before.id)

        assert (
            otu_after.plan.segments[0].length_tolerance
            == otu_before.plan.segments[0].length_tolerance
        )


class TestUpdateRepresentativeIsolateCommand:
    def test_isolate_id_ok(self, scratch_repo: Repo):
        taxid = 345184

        otu_before = scratch_repo.get_otu_by_taxid(taxid)

        other_isolate_id = None
        for isolate_id in otu_before.isolate_ids:
            if isolate_id != otu_before.representative_isolate:
                other_isolate_id = isolate_id
                break

        assert type(other_isolate_id) is UUID

        subprocess.run(
            ["ref-builder", "otu"]
            + ["--path", str(scratch_repo.path)]
            + ["set-default-isolate"]
            + [str(taxid), str(other_isolate_id)],
            check=False,
        )

        scratch_repo = Repo(scratch_repo.path)

        otu_after = scratch_repo.get_otu_by_taxid(taxid)

        assert otu_after.representative_isolate != otu_before.representative_isolate

        assert otu_after.representative_isolate == other_isolate_id

    def test_isolate_name_ok(self, scratch_repo: Repo):
        taxid = 1169032

        otu_before = scratch_repo.get_otu_by_taxid(taxid)

        representative_isolate_after = None

        for isolate_id in otu_before.isolate_ids:
            if isolate_id != otu_before.representative_isolate:
                representative_isolate_after = otu_before.get_isolate(isolate_id).name
                break

        subprocess.run(
            ["ref-builder", "otu"]
            + ["--path", str(scratch_repo.path)]
            + ["set-default-isolate"]
            + [str(taxid)]
            + [str(representative_isolate_after)],
            check=False,
        )

        scratch_repo = Repo(scratch_repo.path)

        otu_after = scratch_repo.get_otu_by_taxid(taxid)

        assert otu_after.representative_isolate != otu_before.representative_isolate

        assert otu_after.representative_isolate == otu_before.get_isolate_id_by_name(
            representative_isolate_after
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
            ["ref-builder", "otu", "--path", str(empty_repo.path)]
            + ["promote", str(2164102)],
            check=False,
        )

        repo_after = Repo(empty_repo.path)

        otu_after = repo_after.get_otu(otu.id)

        assert otu_after.representative_isolate == otu_before.representative_isolate

        assert otu_after.accessions == {"NC_055390", "NC_055391", "NC_055392"}

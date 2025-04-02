from uuid import UUID, uuid4

import pytest
from syrupy import SnapshotAssertion
from syrupy.filters import props

from ref_builder.otu.create import create_otu_with_taxid
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
from ref_builder.plan import (
    Plan,
    Segment,
    SegmentName,
    SegmentRule,
)
from ref_builder.repo import Repo
from ref_builder.resources import RepoSequence
from ref_builder.utils import IsolateName, IsolateNameType
from tests.fixtures.factories import IsolateFactory


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
                rule=SegmentRule.RECOMMENDED,
            ),
        )

        new_plan.segments.append(
            Segment.new(
                length=1000,
                length_tolerance=scratch_repo.settings.default_segment_length_tolerance,
                name=SegmentName(prefix="DNA", key="Z"),
                rule=SegmentRule.OPTIONAL,
            ),
        )
        with scratch_repo.lock():
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

        with scratch_repo.lock():
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

    @pytest.mark.parametrize(
        ("initial_accessions", "new_accessions"),
        [
            (["MF062136", "MF062137"], ["MF062138"]),
            (["MF062136"], ["MF062137", "MF062138"]),
        ],
    )
    def test_add_segments_to_plan_ok(
        self,
        precached_repo: Repo,
        initial_accessions: list[str],
        new_accessions: list[str],
        snapshot: SnapshotAssertion,
    ):
        """Test the addition of segments to an OTU plan."""
        with precached_repo.lock():
            otu_before = create_otu_with_taxid(
                precached_repo,
                2164102,
                initial_accessions,
                acronym="",
            )

        original_plan = otu_before.plan

        assert isinstance(original_plan, Plan)

        with precached_repo.lock():
            new_segment_ids = add_segments_to_plan(
                precached_repo,
                otu_before,
                rule=SegmentRule.OPTIONAL,
                accessions=new_accessions,
            )

        assert len(new_segment_ids) == len(new_accessions)

        otu_after = precached_repo.get_otu(otu_before.id)

        assert new_segment_ids.issubset(otu_before.plan.segment_ids)

        assert otu_after.plan.segment_ids.issubset(otu_before.plan.segment_ids)

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

        with scratch_repo.lock():
            assert not add_segments_to_plan(
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

        with scratch_repo.lock():
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

        with scratch_repo.lock():
            set_plan_length_tolerances(scratch_repo, otu_before, bad_tolerance)

        otu_after = scratch_repo.get_otu(otu_before.id)

        assert (
            otu_after.plan.segments[0].length_tolerance
            == otu_before.plan.segments[0].length_tolerance
        )


class TestDeleteIsolate:
    """Test isolate deletion behaviour."""

    def test_ok(self, scratch_repo):
        """Test that a given isolate can be deleted from the OTU."""
        taxid = 1169032

        otu_before = scratch_repo.get_otu_by_taxid(taxid)

        isolate_id = otu_before.get_isolate_id_by_name(
            IsolateName(type=IsolateNameType.ISOLATE, value="WMoV-6.3"),
        )

        assert type(isolate_id) is UUID

        with scratch_repo.lock():
            assert delete_isolate_from_otu(scratch_repo, otu_before, isolate_id)

        otu_after = scratch_repo.get_otu_by_taxid(taxid)

        assert isolate_id not in otu_after.isolate_ids

        assert otu_after.get_isolate(isolate_id) is None

        assert otu_before.get_isolate(isolate_id).accessions not in otu_after.accessions

        assert len(otu_after.isolate_ids) == len(otu_before.isolate_ids) - 1

    def test_representative_isolate_fail(self, scratch_repo: Repo):
        """Test that the representative isolate cannot be deleted."""
        taxid = 1169032

        otu_before = scratch_repo.get_otu_by_taxid(taxid)

        with scratch_repo.lock():
            assert not delete_isolate_from_otu(
                scratch_repo,
                otu_before,
                isolate_id=otu_before.representative_isolate,
            )

        assert scratch_repo.get_otu(otu_before.id).isolate_ids == otu_before.isolate_ids


class TestReplaceSequence:
    def test_ok(self, precached_repo):
        """Test sequence replacement and deletion."""
        with precached_repo.lock():
            otu_init = create_otu_with_taxid(
                precached_repo,
                1169032,
                ["MK431779"],
                acronym="",
            )

        old_sequence = otu_init.get_sequence_by_accession("MK431779")

        isolate_id = next(
            iter(otu_init.get_isolate_ids_containing_sequence_id(old_sequence.id))
        )

        assert isolate_id == otu_init.representative_isolate

        with precached_repo.lock():
            new_sequence = replace_sequence_in_otu(
                repo=precached_repo,
                otu=otu_init,
                new_accession="NC_003355",
                replaced_accession="MK431779",
            )

        assert isinstance(new_sequence, RepoSequence)

        otu_after = precached_repo.get_otu_by_taxid(1169032)

        assert (
            otu_after.accessions
            == otu_after.get_isolate(isolate_id).accessions
            == {"NC_003355"}
        )

    def test_multiple_link_ok(self, precached_repo):
        with precached_repo.lock():
            otu_init = create_otu_with_taxid(
                precached_repo,
                345184,
                ["DQ178608", "DQ178609"],
                acronym="",
            )

        otu_id = otu_init.id

        post_init_otu = precached_repo.get_otu(otu_init.id)

        assert post_init_otu.accessions == {"DQ178608", "DQ178609"}

        rep_isolate = post_init_otu.get_isolate(post_init_otu.representative_isolate)

        mock_isolate = IsolateFactory.build_on_plan(otu_init.plan)
        mock_sequence = mock_isolate.sequences[1]

        print(mock_sequence)

        print(post_init_otu.plan.get_segment_by_id(mock_sequence.segment))

        with precached_repo.lock(), precached_repo.use_transaction():
            sequence_seg2 = precached_repo.create_sequence(
                otu_id=otu_id,
                accession=str(mock_sequence.accession),
                definition=mock_sequence.definition,
                legacy_id=None,
                segment=mock_sequence.segment,
                sequence=mock_sequence.sequence,
            )

            isolate_init = precached_repo.create_isolate(
                otu_id=otu_id,
                legacy_id=None,
                name=mock_isolate.name,
            )

            precached_repo.link_sequence(
                otu_id=otu_id,
                isolate_id=isolate_init.id,
                sequence_id=rep_isolate.sequences[0].id,
            )

            precached_repo.link_sequence(
                otu_id=otu_id,
                isolate_id=isolate_init.id,
                sequence_id=sequence_seg2.id,
            )

        otu_second_isolate = precached_repo.get_otu(otu_id)

        assert otu_second_isolate.accessions == {
            "DQ178608",
            "DQ178609",
            sequence_seg2.accession.key,
        }

        with precached_repo.lock():
            replace_sequence_in_otu(
                precached_repo,
                otu_second_isolate,
                new_accession="NC_038792",
                replaced_accession="DQ178608",
            )

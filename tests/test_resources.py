from uuid import uuid4

import pytest

from ref_builder.models import Molecule, MolType, Strandedness, Topology
from ref_builder.plan import Plan, PlanWarning, Segment, SegmentRule
from ref_builder.repo import Repo
from ref_builder.resources import RepoIsolate, RepoOTU
from ref_builder.utils import IsolateName, IsolateNameType
from tests.fixtures.factories import IsolateFactory, OTUFactory


class TestSequence:
    """Test properties of RepoSequence."""

    @pytest.mark.parametrize(
        ("taxid", "accessions"),
        [
            (
                345184,
                ["DQ178614", "DQ178613", "DQ178610", "DQ178611"],
            ),
        ],
    )
    def test_equivalence(self, taxid: int, accessions: list["str"], scratch_repo: Repo):
        """Test that the == operator works correctly."""
        otu = scratch_repo.get_otu_by_taxid(taxid)

        for accession in accessions:
            assert otu.get_sequence_by_accession(
                accession,
            ) == otu.get_sequence_by_accession(accession)


class TestIsolate:
    """Test properties of RepoIsolate."""

    def test_no_sequences(self):
        """Test that an isolate intializes correctly with no sequences."""
        isolate = RepoIsolate(
            id=uuid4(),
            legacy_id=None,
            name=IsolateName(type=IsolateNameType.ISOLATE, value="A"),
            sequences=[],
        )

        assert isolate.sequences == []

    @pytest.mark.parametrize("taxid", [345184])
    def test_equivalence(self, taxid: int, scratch_repo: Repo):
        """Test that the == operator works correctly."""
        otu = scratch_repo.get_otu_by_taxid(taxid)

        for isolate in otu.isolates:
            assert isolate == otu.get_isolate(isolate.id)


class TestOTU:
    """Test properties of RepoOTU."""

    def test_no_isolates(self):
        """Test that an isolate initializes correctly with no isolates."""
        otu = RepoOTU(
            id=uuid4(),
            acronym="TMV",
            excluded_accessions=set(),
            isolates=[],
            legacy_id=None,
            molecule=Molecule(
                strandedness=Strandedness.SINGLE,
                type=MolType.RNA,
                topology=Topology.LINEAR,
            ),
            name="Tobacco mosaic virus",
            representative_isolate=None,
            plan=Plan.new(
                segments=[
                    Segment.new(
                        length=6395,
                        length_tolerance=0.03,
                        name=None,
                    )
                ]
            ),
            taxid=12242,
        )

        assert otu.isolates == []

    def test_equivalence(self, scratch_repo: Repo):
        """Test that the == operator works correctly."""
        taxid = 345184

        assert scratch_repo.get_otu_by_taxid(taxid) == scratch_repo.get_otu_by_taxid(
            taxid,
        )

    def test_get_sequence_id_hierarchy(self, scratch_repo: Repo):
        """Test that the isolate ID can be found from a sequence ID."""
        otu = scratch_repo.get_otu_by_taxid(345184)

        test_sequence = otu.get_sequence_by_accession("DQ178610")

        containing_isolate_ids = otu.get_isolate_ids_containing_sequence_id(
            test_sequence.id
        )

        assert len(containing_isolate_ids) == 1

        isolate_id = next(iter(containing_isolate_ids))

        assert otu.get_isolate(isolate_id) is not None

    def test_check_get_sequence_by_id_integrity(self, scratch_repo: Repo):
        """Test that RepoOTU.get_sequence() can retrieve every sequence ID
        within each constituent isolate from its own private index.
        """

        otu_ = RepoOTU(**OTUFactory.build().model_dump())

        for i in range(10):
            isolate = RepoIsolate.model_validate(
                IsolateFactory.build_on_plan(otu_.plan).model_dump()
            )
            otu_.add_isolate(isolate)

        sequence_ids_in_isolates = {
            sequence.id for isolate in otu_.isolates for sequence in isolate.sequences
        }

        for sequence_id in sequence_ids_in_isolates:
            assert otu_.get_sequence_by_id(sequence_id) is not None

    def test_check_get_sequence_by_accession_integrity(self, scratch_repo: Repo):
        """Test that RepoOTU.get_sequence() can retrieve every accession
        within each constituent isolate.
        """

        otu_ = RepoOTU.model_validate(OTUFactory.build().model_dump())

        for i in range(10):
            isolate = RepoIsolate.model_validate(
                IsolateFactory.build_on_plan(otu_.plan).model_dump()
            )
            otu_.add_isolate(isolate)

        for isolate in otu_.isolates:
            for sequence in isolate.sequences:
                assert (
                    otu_.get_sequence_by_accession(sequence.accession.key) == sequence
                )

    def test_otu_delete_isolate(self):
        """Test that RepoOTU.delete_isolate() does not affect other isolates."""

        otu_before = RepoOTU.model_validate(OTUFactory.build().model_dump())

        for i in range(10):
            isolate = RepoIsolate.model_validate(
                IsolateFactory.build_on_plan(otu_before.plan).model_dump()
            )
            otu_before.add_isolate(isolate)

        initial_isolate_ids = [isolate.id for isolate in otu_before.isolates]

        deleted_id = initial_isolate_ids[3]

        otu_after = otu_before.model_copy()

        otu_after.delete_isolate(deleted_id)

        assert otu_after.get_isolate(deleted_id) is None

        for i in range(10):
            extant_id = initial_isolate_ids[i + 3]
            assert otu_after.get_isolate(extant_id) == otu_before.get_isolate(extant_id)

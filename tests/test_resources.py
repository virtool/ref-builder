from uuid import uuid4

import pytest

from ref_builder.models import Molecule, MolType, Strandedness, Topology
from ref_builder.plan import MonopartitePlan
from ref_builder.repo import Repo
from ref_builder.resources import RepoIsolate, RepoOTU
from ref_builder.utils import IsolateName, IsolateNameType


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
            repr_isolate=None,
            plan=MonopartitePlan.new(length=6395, length_tolerance=0.03),
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

        isolate_id, sequence_id = otu.get_sequence_id_hierarchy_from_accession(
            "DQ178610",
        )

        assert otu.get_isolate(isolate_id) is not None

        assert sequence_id in otu.sequence_ids

from uuid import uuid4

import pytest

from ref_builder.models import Molecule, MolType, Strandedness, Topology
from ref_builder.resources import RepoIsolate, RepoOTU
from ref_builder.schema import OTUSchema, Segment
from ref_builder.utils import IsolateName, IsolateNameType


class TestSequence:
    @pytest.mark.parametrize(
        "taxid,accessions",
        [
            (
                345184,
                ["DQ178614", "DQ178613", "DQ178610", "DQ178611"],
            ),
        ],
    )
    def test_equivalence(self, taxid, accessions, scratch_repo):
        """Test that the == operator works correctly."""
        otu = scratch_repo.get_otu_by_taxid(taxid)

        for accession in accessions:
            assert otu.get_sequence_by_accession(
                accession,
            ) == otu.get_sequence_by_accession(accession)


class TestIsolate:
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
    def test_equivalence(self, taxid, scratch_repo):
        """Test that the == operator works correctly."""
        otu = scratch_repo.get_otu_by_taxid(taxid)

        for isolate in otu.isolates:
            assert isolate == otu.get_isolate(isolate.id)


class TestOTU:
    def test_no_isolates(self):
        """Test that an isolate initializes correctly with no isolates."""
        otu = RepoOTU(
            uuid=uuid4(),
            taxid=12242,
            name="Tobacco mosaic virus",
            schema=OTUSchema(
                molecule=Molecule(
                    strandedness=Strandedness.SINGLE,
                    type=MolType.RNA,
                    topology=Topology.LINEAR,
                ),
                segments=[
                    Segment(id=uuid4(), name="genomic RNA", length=6395, required=True),
                ],
            ),
        )

        assert otu.isolates == []

    @pytest.mark.parametrize("taxid", [345184])
    def test_equivalence(self, taxid, scratch_repo):
        """Test that the == operator works correctly."""
        otu = scratch_repo.get_otu_by_taxid(taxid)

        otu_copy = RepoOTU(
            uuid=otu.id,
            taxid=otu.taxid,
            name=otu.name,
            acronym=otu.acronym,
            legacy_id=otu.legacy_id,
            schema=otu.schema,
            excluded_accessions=otu.excluded_accessions,
            isolates=otu.isolates,
            repr_isolate=otu.repr_isolate,
        )

        assert otu == otu_copy

    @pytest.mark.parametrize("taxid, accession", [(345184, "DQ178610")])
    def test_get_sequence_id_hierarchy(self, taxid, accession, scratch_repo):
        otu = scratch_repo.get_otu_by_taxid(taxid)

        isolate_id, sequence_id = otu.get_sequence_id_hierarchy_from_accession(
            accession,
        )

        assert otu.get_isolate(isolate_id) is not None

        assert sequence_id in otu.sequence_ids

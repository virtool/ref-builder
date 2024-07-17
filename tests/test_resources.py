from uuid import uuid4

import pytest

from ref_builder.resources import (
    RepoIsolate,
    RepoOTU,
    RepoSequence,
)
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
            sequence = otu.get_sequence_by_accession(accession)

            assert type(sequence) is RepoSequence

            sequence_copy = RepoSequence(**sequence.dict())

            assert sequence == sequence_copy


class TestIsolate:
    def test_no_sequences(self):
        """Test that an isolate intializes correctly with no sequences."""
        isolate = RepoIsolate(
            uuid=uuid4(),
            name=IsolateName(type=IsolateNameType.ISOLATE, value="A"),
        )

        assert isolate.sequences == []

    @pytest.mark.parametrize("taxid", [345184])
    def test_equivalence(self, taxid, scratch_repo):
        """Test that the == operator works correctly."""
        otu = scratch_repo.get_otu_by_taxid(taxid)

        for isolate in otu.isolates:
            isolate_copy = RepoIsolate(
                uuid=isolate.id,
                name=isolate.name,
                legacy_id=isolate.legacy_id,
                sequences=isolate.sequences,
            )

            assert isolate == isolate_copy


class TestOTU:
    def test_no_isolates(self):
        """Test that an isolate intializes correctly with no isolates."""
        otu = RepoOTU(
            uuid=uuid4(),
            taxid=12242,
            name="Tobacco mosaic virus",
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
        )

        assert otu == otu_copy

"""Tests for OTU models."""

import uuid
import warnings

import pytest
from pydantic import ValidationError

from ref_builder.otu.models import OTU, Isolate, OTUBase, Sequence
from ref_builder.plan import PlanWarning, SegmentRule
from ref_builder.utils import Accession, is_refseq
from tests.fixtures.factories import (
    IsolateFactory,
    OTUFactory,
    PlanFactory,
    SegmentFactory,
    SequenceFactory,
)


class TestSequence:
    """Test the ``Sequence`` model which is used for complete validation of sequences."""

    def test_ok(self, sequence_factory: SequenceFactory):
        sequence = sequence_factory.build()

        assert Sequence.model_validate(sequence.model_dump())

    def test_bad_sequence_fail(self, sequence_factory: SequenceFactory):
        bad_sequence = sequence_factory.build(accession=Accession("BADSEQUENCE", 1))

        try:
            assert Sequence.model_validate(bad_sequence.model_dump())

        except ValidationError as e:
            for error in e.errors():
                assert "accession" in error["loc"]
                assert (
                    "Accession BADSEQUENCE.1 does not match a valid accession pattern"
                    in error["msg"]
                )


class TestIsolate:
    """Test the ``Isolate`` model which is used for complete validation of isolates."""

    otu: OTUBase

    @pytest.fixture(autouse=True)
    def _build_otu(self, otu_factory: OTUFactory):
        self.otu = otu_factory.build(
            plan=PlanFactory.build(segments=SegmentFactory.build_series(3))
        )

    def test_accession_consistency_warning(self, isolate_factory: IsolateFactory):
        """Test that a warning is raised if accession provenances are mixed."""
        genbank_isolate = isolate_factory.build_on_plan(self.otu.plan)
        refseq_isolate = isolate_factory.build_on_plan(self.otu.plan, refseq=True)

        for isolate in (genbank_isolate, refseq_isolate):
            bad_isolate = isolate.model_copy()

            if is_refseq(bad_isolate.sequences[0].accession.key):
                bad_isolate.sequences[0].accession = Accession("BD000001", 1)
            else:
                bad_isolate.sequences[0].accession = Accession("NC_000001", 1)

            assert not all(
                [
                    is_refseq(sequence.accession.key)
                    for sequence in bad_isolate.sequences
                ]
            )

            with warnings.catch_warnings(record=True) as warning_list:
                Isolate.model_validate(bad_isolate.model_dump())

            for warning_msg in warning_list:
                assert warning_msg.category.__name__ == "IsolateInconsistencyWarning"

                assert (
                    "Combination of RefSeq and non-RefSeq sequences found in multipartite isolate"
                    in str(warning_msg.message)
                )

                assert (
                    f"{[sequence.accession.key for sequence in bad_isolate.sequences]}"
                    in str(warning_msg)
                )


class TestOTU:
    """Test the ``OTU`` model which is used for complete validation of OTUs."""

    otu: OTUBase

    @pytest.fixture(autouse=True)
    def _build_otu(self, otu_factory: OTUFactory):
        self.otu = otu_factory.build()

    def test_ok(self):
        """Test that a valid OTU passes validation."""
        assert OTU.model_validate(self.otu.model_dump())

    def test_no_required_segments(self):
        """Test that OTU raises a warning if initialized without required segments."""
        mock_otu = OTUFactory.build(
            plan=PlanFactory.build(
                segments=[SegmentFactory.build(rule=SegmentRule.RECOMMENDED)]
            )
        )

        assert not mock_otu.plan.required_segments

        with pytest.warns(PlanWarning):
            OTU.model_validate(mock_otu.model_dump())

    def test_excluded_accessions(self):
        """Test that validation fails if the OTU includes accessions that are included
        in ``excluded_accessions``.
        """
        accession = self.otu.isolates[0].sequences[0].accession.key

        self.otu.excluded_accessions.add(accession)

        with pytest.raises(
            ValueError, match=f"Excluded accessions found in the OTU: {accession}"
        ):
            assert OTU.model_validate(self.otu.model_dump())

    def test_no_isolates(self):
        """Test that validation fails if the OTU has no isolates."""
        self.otu.isolates = []

        with pytest.raises(
            ValueError, match="List should have at least 1 item after validation, not 0"
        ):
            assert OTU.model_validate(self.otu.model_dump())

    def test_no_representative_isolate(self):
        """Test that validation fails if the OTU has no representative isolate."""
        self.otu.representative_isolate = None

        with pytest.raises(
            ValueError,
            match="UUID input should be a string, bytes or UUID object.*input_value="
            "None",
        ):
            assert OTU.model_validate(self.otu.model_dump())

    def test_missing_representative_isolate(self):
        """Test that validation fails if the representative isolate is not in the
        OTU.
        """
        self.otu.representative_isolate = uuid.uuid4()

        with pytest.raises(
            ValueError, match="Representative isolate must be in the OTU"
        ):
            assert OTU.model_validate(self.otu.model_dump())

    def test_unique_isolate_names(self):
        """Test that validation fails if the OTU has duplicate isolate names."""
        duplicate_name = self.otu.isolates[1].name = self.otu.isolates[0].name

        with pytest.raises(
            ValueError,
            match=f"Isolate names must be unique. Non-unique names: {duplicate_name}",
        ):
            assert OTU.model_validate(self.otu.model_dump())

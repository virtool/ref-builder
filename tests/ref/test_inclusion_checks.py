from uuid import UUID

import pytest
from faker import Faker

from ref_builder.resources import RepoOTU
from ref_builder.ncbi.models import NCBIGenbank
from ref_builder.otu.isolate import (
    create_monopartite_isolate,
    create_multipartite_isolate,
)
from ref_builder.otu.utils import check_sequence_length, IsolateName, IsolateNameType
from tests.fixtures.factories import NCBIGenbankFactory, NCBISourceFactory
from tests.fixtures.providers import AccessionProvider, SequenceProvider


faker = Faker()
faker.add_provider(AccessionProvider)
faker.add_provider(SequenceProvider)

REF_STRING = "ABCDEFGH"


def create_mock_record(
    otu: RepoOTU,
    sequence_length: int,
) -> NCBIGenbank:
    """Return a mock records based on a monopartite OTU."""
    mock_record_source = NCBISourceFactory.build(
        taxid=otu.taxid,
        organism=otu.name,
    )

    mock_record = NCBIGenbankFactory.build(
        accession=faker.genbank_accession(),
        source=mock_record_source,
        sequence=faker.sequence(min=sequence_length, max=sequence_length),
    )

    return mock_record


def create_mock_isolate_records(
    otu: RepoOTU,
    sequence_length_multiplier: float,
) -> dict[UUID, NCBIGenbank]:
    """Return a dictionary of mock records based on a multipartite OTU."""
    mock_assigned_records = {}

    accession_starter = 100000

    for i in range(len(otu.plan.required_segments)):
        segment = otu.plan.required_segments[i]

        mock_record_source = NCBISourceFactory.build(
            taxid=otu.taxid,
            organism=otu.name,
            segment=str(segment.name),
        )

        sequence_length = int(segment.length * sequence_length_multiplier)

        mock_assigned_records[segment.id] = NCBIGenbankFactory.build(
            accession=f"AB{accession_starter + i}",
            source=mock_record_source,
            sequence=faker.sequence(
                min=sequence_length, max=sequence_length
            ),
        )

    return mock_assigned_records


class TestSequenceLengthCheck:
    """Test sequence length check functions."""

    @pytest.mark.parametrize(
        ("sequence_lengths", "tolerance"),
        [
            ([1000, 970, 1030], 0.03),
            ([1000, 980, 1020], 0.02),
            ([1000, 960, 1040], 0.04),
        ],
    )
    def test_check_sequence_length_function_success(
        self, sequence_lengths: list[int], tolerance: float
    ):
        """Test that check_sequence_length() returns True
        given different sequences and tolerance levels.
        """
        for sequence_length in sequence_lengths:
            fake_sequence = faker.sequence(min=sequence_length, max=sequence_length)

            assert check_sequence_length(fake_sequence, 1000, tolerance)

    @pytest.mark.parametrize(
        ("sequence_lengths", "tolerance"),
        [
            ([970, 1030], 0.02),
            ([969, 1031], 0.03),
            ([979, 1021], 0.02),
            ([959, 1041], 0.04),
        ],
    )
    def test_check_sequence_length_function_fail(
        self, sequence_lengths: list[int], tolerance: float
    ):
        """Test that check_sequence_length() returns False
        given different sequences and tolerance levels.
        """
        for sequence_length in sequence_lengths:
            fake_sequence = faker.sequence(min=sequence_length, max=sequence_length)

            assert not check_sequence_length(fake_sequence, 1000, tolerance)


class TestAddMonopartiteIsolate:
    """Test that new monopartite isolates follow plan length limits."""

    @pytest.mark.parametrize("sequence_length_multiplier", [1.0, 1.03, 0.98])
    def test_ok(self, scratch_repo, sequence_length_multiplier: float):
        """Test that sequences within recommended length variance
        are added without issue."""
        otu_before = scratch_repo.get_otu_by_taxid(270478)

        mock_record = create_mock_record(
            otu_before,
            int(otu_before.plan.segments[0].length * sequence_length_multiplier),
        )

        assert mock_record

        mock_isolate = create_monopartite_isolate(
            scratch_repo,
            otu_before,
            isolate_name=IsolateName(type=IsolateNameType.ISOLATE, value="mock"),
            record=mock_record,
        )

        otu_after = scratch_repo.get_otu_by_taxid(270478)

        assert mock_isolate.id in otu_after.isolate_ids

        assert mock_record.accession in otu_after.accessions

    @pytest.mark.parametrize("sequence_length_multiplier", [0.5, 1.035, 20.0])
    def test_fail(self, scratch_repo, sequence_length_multiplier: float):
        """Test that sequences that exceed recommended length variance
        are automatically rejected.
        """
        otu_before = scratch_repo.get_otu_by_taxid(270478)

        mock_record = create_mock_record(
            otu_before,
            int(otu_before.plan.segments[0].length * sequence_length_multiplier),
        )

        assert mock_record

        mock_isolate = create_monopartite_isolate(
            scratch_repo,
            otu_before,
            isolate_name=IsolateName(type=IsolateNameType.ISOLATE, value="mock"),
            record=mock_record,
        )

        assert mock_isolate is None

        otu_after = scratch_repo.get_otu_by_taxid(270478)

        assert mock_record.accession not in otu_after.accessions


class TestAddMultipartiteIsolate:
    """Test that new multipartite isolates follow plan length limits."""

    @pytest.mark.parametrize("sequence_length_multiplier", [1.0, 1.03, 0.98])
    def test_ok(self, scratch_repo, sequence_length_multiplier: float):
        """Test that sequences within recommended length variance
        are added without issue."""
        otu_before = scratch_repo.get_otu_by_taxid(3158377)

        mock_assigned_records = create_mock_isolate_records(
            otu_before,
            sequence_length_multiplier,
        )

        assert mock_assigned_records

        mock_isolate = create_multipartite_isolate(
            scratch_repo,
            otu_before,
            isolate_name=IsolateName(type=IsolateNameType.ISOLATE, value="mock"),
            assigned_records=mock_assigned_records,
        )

        otu_after = scratch_repo.get_otu_by_taxid(3158377)

        assert mock_isolate.id in otu_after.isolate_ids

        assert mock_isolate.accessions.issubset(otu_after.accessions)

    @pytest.mark.parametrize("sequence_length_multiplier", [0.5, 1.035, 20.0])
    def test_fail(self, scratch_repo, sequence_length_multiplier: float):
        """Test that sequences that exceed recommended length variance
        are automatically rejected.
        """
        otu_before = scratch_repo.get_otu_by_taxid(3158377)

        mock_assigned_records = create_mock_isolate_records(
            otu_before,
            sequence_length_multiplier,
        )

        assert mock_assigned_records

        mock_isolate = create_multipartite_isolate(
            scratch_repo,
            otu_before,
            isolate_name=IsolateName(type=IsolateNameType.ISOLATE, value="mock"),
            assigned_records=mock_assigned_records,
        )

        assert mock_isolate is None

        otu_after = scratch_repo.get_otu_by_taxid(3158377)

        assert otu_after.isolate_ids == otu_before.isolate_ids

        assert otu_after.accessions == otu_before.accessions

import pytest
from faker import Faker

from ref_builder.otu.utils import check_sequence_length
from tests.fixtures.providers import SequenceProvider


faker = Faker()
faker.add_provider(SequenceProvider)


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

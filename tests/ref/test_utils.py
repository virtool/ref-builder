from pathlib import Path

import pytest

from ref_builder.utils import format_json, pad_zeroes


class TestPadZeroes:
    def test_ok(self):
        """Test that integers are padded with zeroes to 8 digit strings."""
        assert pad_zeroes(112) == "00000112"
        assert pad_zeroes(1) == "00000001"
        assert pad_zeroes(87654321) == "87654321"

    @pytest.mark.parametrize("number", [100000000, 12345678901234567890])
    def test_too_big(self, number: int):
        """Test that integers larger than 8 digits are returned as strings."""
        with pytest.raises(ValueError) as e:
            pad_zeroes(number)

        assert str(e.value) == "Number is too large to pad"


def test_format_json(tmp_path: Path):
    """Test that a mis-formatted JSON file is formatted in place."""
    with open(tmp_path / "test.json", "w") as f:
        f.write(
            '{"name": "bob",    "active": true, "greeting": "hello",'
            '"street": "downing", \n\n "number": 10}',
        )

    format_json(tmp_path / "test.json")

    assert (tmp_path / "test.json").read_text() == (
        "{\n"
        '  "name": "bob",\n'
        '  "active": true,\n'
        '  "greeting": "hello",\n'
        '  "street": "downing",\n'
        '  "number": 10\n'
        "}"
    )

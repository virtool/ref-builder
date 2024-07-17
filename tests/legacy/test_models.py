import pytest

from ref_builder.legacy.models import (
    LegacyIsolate,
    LegacyOTU,
    LegacySchemaSegment,
    LegacySequence,
)


@pytest.fixture()
def _example_sequence():
    return {
        "_id": "abcd1234",
        "accession": "aj18907.1 ",
        "definition": "This is a sufficiently long definition ",
        "host": "",
        "segment": "RNA 1",
        "sequence": "ASTTCGRKVNACTTATchTACGAGCACYTAGACM   ",
    }


@pytest.fixture()
def _example_isolate(_example_sequence: dict):
    return {
        "id": "bnmh9021",
        "default": True,
        "sequences": [_example_sequence],
        "source_type": "isolate",
        "source_name": "Pine",
    }


@pytest.fixture()
def _example_otu(_example_isolate: dict):
    return {
        "_id": "abcd1234",
        "abbreviation": "GymV1  ",
        "isolates": [_example_isolate],
        "name": "Gymnosperm virus 1 ",
        "taxid": 12345,
        "schema": [
            {
                "molecule": "ssRNA",
                "name": "RNA 1",
                "required": True,
            },
        ],
    }


class TestIsolate:
    def test_ok(self, _example_isolate: dict):
        """Test that a valid isolate passes."""
        isolate = LegacyIsolate(**_example_isolate)

        assert isolate.id == "bnmh9021"
        assert isolate.default is True
        assert len(isolate.sequences) == 1
        assert isolate.source_type == "isolate"
        assert isolate.source_name == "Pine"

    def test_source_type_invalid(self, _example_isolate: dict):
        """Test that an isolate with an invalid source type raises a ValueError."""
        _example_isolate["source_type"] = "invalid"

        with pytest.raises(ValueError) as e:
            LegacyIsolate(**_example_isolate)

        assert (
            "Input should be 'clone', 'culture', 'isolate', 'genbank', 'genotype', "
            "'serotype', 'strain' or 'variant'" in str(e.value)
        )

    def test_source_name_empty(self, _example_isolate: dict):
        """Test that an isolate with an empty source name raises a ValueError."""
        _example_isolate["source_name"] = ""

        with pytest.raises(ValueError) as e:
            LegacyIsolate(**_example_isolate)

        assert "String should have at least 1 character" in str(e.value)

    def test_too_few_sequences(self, _example_isolate: dict):
        """Test that an isolate with no sequences raises a ValueError."""
        _example_isolate["sequences"] = []

        with pytest.raises(ValueError) as e:
            LegacyIsolate(**_example_isolate)

        assert "List should have at least 1 item after validation, not 0" in str(
            e.value,
        )


class TestSequence:
    def test_ok(self, _example_sequence: dict):
        """"""
        sequence = LegacySequence(**_example_sequence)

        assert sequence.id == "abcd1234"
        assert sequence.accession == "AJ18907.1"
        assert sequence.definition == "This is a sufficiently long definition"
        assert sequence.host == ""
        assert sequence.segment == "RNA 1"
        assert sequence.sequence == "ASTTCGRKVNACTTATCHTACGAGCACYTAGACM"

    def test_definition_too_short(self, _example_sequence: dict):
        """Test that a definition that is too short raises a ValueError."""
        _example_sequence["definition"] = "short"

        with pytest.raises(ValueError) as e:
            LegacySequence(**_example_sequence)

        assert "String should have at least 10 characters" in str(e.value)

    def test_accession_without_dot(self, _example_sequence: dict):
        """Test that an accession with no dot raises a ValueError."""
        _example_sequence["accession"] = "aj18907"

        with pytest.raises(ValueError) as e:
            LegacySequence(**_example_sequence)

        assert "String should match pattern '^\\w+\\.\\d+$'" in str(e.value)

    def test_sequence_bad_characters(self, _example_sequence: dict):
        """Test that a sequence containing a 'Z' raises a ValueError."""
        _example_sequence["sequence"] = "ASTTCGSZATAGAGAGM"

        with pytest.raises(ValueError) as e:
            LegacySequence(**_example_sequence)

        assert (
            "String should match pattern '^[ATCGRYKMSWBDHVNatcgrykmswbdhvn]+$'"
            in str(e.value)
        )


class TestOTU:
    def test_ok(self, _example_otu: dict):
        """Test that a valid OTU passes and whitespace is stripped."""
        otu = LegacyOTU(**_example_otu)

        assert otu.id == "abcd1234"
        assert otu.abbreviation == "GymV1"
        assert otu.name == "Gymnosperm virus 1"
        assert otu.otu_schema == [
            LegacySchemaSegment(
                molecule="ssRNA",
                name="RNA 1",
                required=True,
            ),
        ]
        assert otu.taxid == 12345

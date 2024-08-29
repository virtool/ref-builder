from uuid import UUID

import pytest
from syrupy import SnapshotAssertion
from syrupy.filters import props

from ref_builder.otu.update import create_schema_from_records
from ref_builder.schema import SegmentName, OTUSchema, parse_segment_name


@pytest.mark.parametrize(
    "accessions",
    [
        ["NC_024301"],
        [
            "NC_010314",
            "NC_010315",
            "NC_010316",
            "NC_010317",
            "NC_010318",
            "NC_010319",
        ],
    ],
)
def test_create_schema_from_records(
    accessions: list[str],
    scratch_ncbi_client,
    snapshot: SnapshotAssertion,
):
    records = scratch_ncbi_client.fetch_genbank_records(accessions)
    auto_schema = create_schema_from_records(records)

    assert type(auto_schema) is OTUSchema

    assert auto_schema.model_dump() == snapshot(exclude=props("id"))

    for segment in auto_schema.segments:
        assert type(segment.id) is UUID


class TestSegmentNameParser:
    @pytest.mark.parametrize(
        "expected_result, test_strings",
        [
            (SegmentName(prefix="DNA", key="A"), ["DNA A", "DNA_A", "DNA-A"]),
            (SegmentName(prefix="RNA", key="BN"), ["RNA BN", "RNA_BN", "RNA-BN"]),
            (SegmentName(prefix="DNA", key="U3"), ["DNA U3", "DNA U3", "DNA-U3"]),
        ],
    )
    def test_ok(self, expected_result: str, test_strings: list[str]):
        assert (
            parse_segment_name(test_strings[0])
            ==
            parse_segment_name(test_strings[1])
            ==
            parse_segment_name(test_strings[2])
            ==
            expected_result
        )

    @pytest.mark.parametrize(
        "fail_case",
        ["", "*V/", "51f9a0bc-7b3b-434f-bf4c-f7abaa015b8d"]
    )
    def test_fail(self, fail_case: str):
        with pytest.raises(ValueError):
            parse_segment_name(fail_case)

import pytest
from syrupy import SnapshotAssertion

from ref_builder.otu import create_schema_from_records
from ref_builder.schema import OTUSchema


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

    assert auto_schema.model_dump() == snapshot

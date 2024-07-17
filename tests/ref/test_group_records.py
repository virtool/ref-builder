import pytest
from syrupy import SnapshotAssertion

from ref_builder.otu import group_genbank_records_by_isolate


class TestGroupRecords:
    @pytest.mark.parametrize(
        "accessions",
        [
            ["AB017504", "MH200607", "MK431779", "NC_003355"],
            ["NC_036587", "MT240513", "MT240490"],
            ["NC_003493", "Y11023", "GU256539"],
        ],
    )
    def test_group_records_by_isolate_success(
        self,
        accessions,
        scratch_ncbi_client,
        snapshot: SnapshotAssertion,
    ):
        records = scratch_ncbi_client.fetch_genbank_records(accessions)

        assert records

        grouped_records = group_genbank_records_by_isolate(records)

        assert grouped_records

        for source_key in grouped_records:
            assert source_key == snapshot
            assert grouped_records[source_key] == snapshot

    @pytest.mark.parametrize("accessions", [["Y11023"]])
    def test_group_records_by_isolate_failure(self, accessions, scratch_ncbi_client):
        records = scratch_ncbi_client.fetch_genbank_records(accessions)

        assert records

        grouped_records = group_genbank_records_by_isolate(records)

        assert not grouped_records

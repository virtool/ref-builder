import pytest
from syrupy import SnapshotAssertion
from syrupy.filters import props

from ref_builder.ncbi.client import NCBIClient
from ref_builder.ncbi.models import NCBISourceMolType
from ref_builder.otu.utils import create_plan_from_records
from ref_builder.plan import (
    Plan,
    SegmentName,
    extract_segment_name_from_record,
    parse_segment_name,
)
from tests.fixtures.factories import NCBIGenbankFactory, NCBISourceFactory


class TestPlan:
    """Test differentiated isolate plans."""

    def test_create_monopartite_plan_from_records(
        self,
        scratch_ncbi_client: NCBIClient,
        snapshot: SnapshotAssertion,
    ):
        """Test the creation of a monopartite plan from a given record."""
        records = scratch_ncbi_client.fetch_genbank_records(["NC_024301"])
        auto_plan = create_plan_from_records(records, length_tolerance=0.03)

        assert auto_plan.model_dump() == snapshot(exclude=props("id"))

    def test_create_multipartite_plan_from_records(
        self,
        scratch_ncbi_client: NCBIClient,
        snapshot: SnapshotAssertion,
    ):
        """Test the creation of a multipartite plan from a set of records
        representing segments.
        """
        records = scratch_ncbi_client.fetch_genbank_records(
            [
                "NC_010314",
                "NC_010315",
                "NC_010316",
                "NC_010317",
                "NC_010318",
                "NC_010319",
            ]
        )
        auto_plan = create_plan_from_records(records, length_tolerance=0.03)

        assert type(auto_plan) is Plan

        assert auto_plan.model_dump() == snapshot(exclude=props("id"))

    @pytest.mark.parametrize(
        ("accessions", "is_monopartite"),
        [
            (["NC_024301"], True),
            (
                [
                    "NC_010314",
                    "NC_010315",
                    "NC_010316",
                    "NC_010317",
                    "NC_010318",
                    "NC_010319",
                ],
                False,
            ),
        ],
    )
    def test_serialize_plan(
        self,
        accessions: list[str],
        is_monopartite: bool,
        scratch_ncbi_client: NCBIClient,
        snapshot: SnapshotAssertion,
    ):
        records = scratch_ncbi_client.fetch_genbank_records(accessions)

        auto_plan = create_plan_from_records(records, length_tolerance=0.03)

        assert auto_plan.model_dump() == snapshot(exclude=props("id"))

        auto_plan_json = auto_plan.model_dump_json()

        auto_plan_from_json = Plan.model_validate_json(auto_plan_json)

        assert auto_plan_from_json.monopartite == is_monopartite

        assert auto_plan_from_json == auto_plan


class TestParseSegmentName:
    """Test segment name normalization."""

    @pytest.mark.parametrize("delimiter", [" ", "_", "-"])
    @pytest.mark.parametrize("key", ["A", "BN", "U3"])
    @pytest.mark.parametrize("prefix", ["DNA", "RNA"])
    def test_ok(self, delimiter: str, key: str, prefix: str):
        """Test that parse_segment_name() correctly parses prefix-key segment names."""
        assert parse_segment_name(f"{prefix}{delimiter}{key}") == SegmentName(
            prefix=prefix, key=key
        )

    @pytest.mark.parametrize(
        "string",
        ["", "*V/", "51f9a0bc-7b3b-434f-bf4c-f7abaa015b8d"],
    )
    def test_none(self, string: str) -> None:
        """Test that parse_segment_name() returns expected ValueError for bad input."""
        assert parse_segment_name(string) is None


class TestExtractSegmentNameFromRecord:
    @pytest.mark.parametrize("key", ["1", "A", "BN", "U3"])
    @pytest.mark.parametrize(
        "mol_type", [NCBISourceMolType.GENOMIC_RNA, NCBISourceMolType.GENOMIC_DNA]
    )
    def test_ok(
        self,
        key: str,
        mol_type: NCBISourceMolType,
        ncbi_genbank_factory: NCBIGenbankFactory,
        ncbi_source_factory: NCBISourceFactory,
    ) -> None:
        """Test that a full segment name can be parsed
        when the source table segment field contains a key with no prefix.
        """
        record = ncbi_genbank_factory.build(
            source=ncbi_source_factory.build(
                mol_type=mol_type,
                segment=key,
            ),
        )

        match mol_type:
            case NCBISourceMolType.GENOMIC_DNA:
                suffix = "DNA"
            case NCBISourceMolType.GENOMIC_RNA:
                suffix = "RNA"
            case _:
                raise ValueError(f"{mol_type} may not be a valid NCBISourceMolType.")

        assert extract_segment_name_from_record(record) == SegmentName(
            prefix=suffix, key=key
        )

    @pytest.mark.parametrize("key", ["", "V#", "51f9a0bc-7b3b-434f-bf4c-f7abaa015b8d"])
    @pytest.mark.parametrize(
        "mol_type", [NCBISourceMolType.GENOMIC_RNA, NCBISourceMolType.GENOMIC_DNA]
    )
    def test_parse_from_record_fail(
        self,
        key: str,
        mol_type: NCBISourceMolType,
        ncbi_genbank_factory: NCBIGenbankFactory,
        ncbi_source_factory: NCBISourceFactory,
    ):
        """Test that ``None`` is returned if no segment name is found."""
        record = ncbi_genbank_factory.build(
            source=ncbi_source_factory.build(
                mol_type=mol_type,
                segment=key,
            ),
        )

        assert extract_segment_name_from_record(record) is None

import pytest
from pydantic import TypeAdapter
from syrupy import SnapshotAssertion
from syrupy.filters import props

from ref_builder.ncbi.client import NCBIClient
from ref_builder.ncbi.models import NCBISourceMolType
from ref_builder.otu.utils import create_isolate_plan_from_records
from ref_builder.plan import (
    Plan,
    MonopartitePlan,
    MultipartitePlan,
    SegmentName,
    get_multipartite_segment_name,
    parse_segment_name,
)
from tests.fixtures.factories import NCBIGenbankFactory, NCBISourceFactory


plan_typer = TypeAdapter(Plan)


class TestPlan:
    """Test differentiated isolate plans."""

    def test_create_monopartite_plan_from_records(
        self,
        scratch_ncbi_client: NCBIClient,
        snapshot: SnapshotAssertion,
    ):
        """Test the creation of a monopartite plan from a given record."""
        records = scratch_ncbi_client.fetch_genbank_records(["NC_024301"])
        auto_plan = create_isolate_plan_from_records(records, length_tolerance=0.03)

        assert type(auto_plan) is MonopartitePlan

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
        auto_plan = create_isolate_plan_from_records(records, length_tolerance=0.03)

        assert type(auto_plan) is MultipartitePlan

        assert auto_plan.model_dump() == snapshot(exclude=props("id"))

    @pytest.mark.parametrize(
        ("accessions", "plan_type"),
        [
            (["NC_024301"], MonopartitePlan),
            (
                [
                    "NC_010314",
                    "NC_010315",
                    "NC_010316",
                    "NC_010317",
                    "NC_010318",
                    "NC_010319",
                ],
                MultipartitePlan,
            ),
        ],
    )
    def test_serialize_plan(
        self,
        accessions: list[str],
        plan_type: MonopartitePlan | MultipartitePlan,
        scratch_ncbi_client: NCBIClient,
        snapshot: SnapshotAssertion,
    ):
        records = scratch_ncbi_client.fetch_genbank_records(accessions)

        auto_plan = create_isolate_plan_from_records(records, length_tolerance=0.03)

        assert auto_plan.model_dump() == snapshot(exclude=props("id"))

        auto_plan_json = auto_plan.model_dump_json()

        reconstituted_plan = plan_typer.validate_json(auto_plan_json)

        assert type(reconstituted_plan) is plan_type

        assert reconstituted_plan == auto_plan


class TestSegmentNameParser:
    """Test segment name normalization."""

    @pytest.mark.parametrize(
        ("expected_result", "test_strings"),
        [
            (SegmentName(prefix="DNA", key="A"), ["DNA A", "DNA_A", "DNA-A"]),
            (SegmentName(prefix="RNA", key="BN"), ["RNA BN", "RNA_BN", "RNA-BN"]),
            (SegmentName(prefix="DNA", key="U3"), ["DNA U3", "DNA U3", "DNA-U3"]),
        ],
    )
    def test_ok(self, expected_result: str, test_strings: list[str]):
        """Test that parse_segment_name() correctly parses prefix-key segment names."""
        assert (
            parse_segment_name(test_strings[0])
            == parse_segment_name(test_strings[1])
            == parse_segment_name(test_strings[2])
            == expected_result
        )

    @pytest.mark.parametrize(
        "fail_case",
        ["", "*V/", "51f9a0bc-7b3b-434f-bf4c-f7abaa015b8d"],
    )
    def test_fail(self, fail_case: str):
        """Test that parse_segment_name() returns expected ValueError for bad input."""
        with pytest.raises(
            ValueError,
            match=r"\s* is not a valid segment name$",
        ):
            parse_segment_name(fail_case)

    @pytest.mark.parametrize(
        ("expected_result", "test_name", "test_moltype"),
        [
            (SegmentName(prefix="DNA", key="A"), "A", NCBISourceMolType.GENOMIC_DNA),
            (SegmentName(prefix="RNA", key="BN"), "BN", NCBISourceMolType.GENOMIC_RNA),
            (SegmentName(prefix="DNA", key="U3"), "U3", NCBISourceMolType.GENOMIC_DNA),
        ],
    )
    def test_parse_from_record_ok(
        self,
        expected_result: SegmentName,
        test_name: str,
        test_moltype: NCBISourceMolType,
    ):
        """Test that a full segment name can be parsed
        when the source table segment field contains a key with no prefix.
        """
        dummy_record = NCBIGenbankFactory.build(
            source=NCBISourceFactory.build(
                mol_type=test_moltype,
                segment=test_name,
            ),
        )

        assert get_multipartite_segment_name(dummy_record) == expected_result

    @pytest.mark.parametrize(
        ("test_name", "test_moltype"),
        [
            ("", NCBISourceMolType.GENOMIC_DNA),
            ("V#", NCBISourceMolType.GENOMIC_RNA),
            ("51f9a0bc-7b3b-434f-bf4c-f7abaa015b8d", NCBISourceMolType.GENOMIC_DNA),
        ],
    )
    def test_parse_from_record_fail(
        self,
        test_name: str,
        test_moltype: NCBISourceMolType,
    ):
        """Test that get_multipartite_segment_name()
        does not return an invalid segment name.
        """
        dummy_record = NCBIGenbankFactory.build(
            source=NCBISourceFactory.build(
                mol_type=test_moltype,
                segment=test_name,
            ),
        )

        with pytest.raises(ValueError, match=r"\s* is not a valid segment name$"):
            get_multipartite_segment_name(dummy_record)

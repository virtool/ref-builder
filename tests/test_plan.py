from uuid import UUID, uuid4

import pytest
from pydantic import ValidationError
from syrupy import SnapshotAssertion

from ref_builder.ncbi.models import NCBISourceMolType
from ref_builder.plan import (
    Plan,
    Segment,
    SegmentName,
    SegmentRule,
    extract_segment_name_from_record,
    parse_segment_name,
)
from tests.fixtures.factories import NCBIGenbankFactory, NCBISourceFactory
from tests.fixtures.utils import uuid_matcher


class TestPlan:
    example: dict

    @pytest.fixture(autouse=True)
    def _setup(self):
        self.example = {
            "id": uuid4(),
            "segments": [
                {
                    "id": uuid4(),
                    "name": None,
                    "rule": "required",
                    "length": 100,
                    "length_tolerance": 0.01,
                }
            ],
        }

    def test_ok(self, snapshot: SnapshotAssertion):
        """Test that a plan can be created."""
        segment = self.example["segments"][0]

        plan = Plan.model_validate(self.example)

        assert (
            Plan(
                id=self.example["id"],
                segments=[
                    Segment(
                        id=segment["id"],
                        name=segment["name"],
                        rule=SegmentRule(segment["required"]),
                        length=segment["length"],
                        length_tolerance=segment["length_tolerance"],
                    )
                ],
            ).model_dump()
            == plan.model_dump()
            == snapshot(matcher=uuid_matcher)
        )

        assert plan.monopartite

    def test_new(self, snapshot: SnapshotAssertion):
        """Test that the new method creates a new plan without a UUID being provided."""
        plan = Plan.new(
            [
                Segment.new(100, 0.01, SegmentName(prefix="DNA", key="A")),
                Segment.new(200, 0.02, SegmentName(prefix="DNA", key="B")),
            ]
        )

        assert isinstance(plan.id, UUID)

        assert plan.model_dump() == snapshot(
            matcher=uuid_matcher,
        )

    @pytest.mark.parametrize("required", ["required", "recommended", "optional"])
    def test_segment_rule(self, required: str):
        """Test that a string segment rule can be parsed into  SegmentRule."""
        plan = Plan.model_validate(
            {
                **self.example,
                "segments": [
                    {
                        **self.example["segments"][0],
                        "required": required,
                    }
                ],
            }
        )

        assert plan.segments[0].rule == SegmentRule(required)

    @pytest.mark.parametrize(
        "name",
        [
            None,
            {"prefix": "DNA", "key": "A                                             "},
        ],
        ids=["two_unnamed", "one_unnamed"],
    )
    def test_multipartite_unnamed(self, name: dict | None):
        """Test that a multipartite path with any unnamed segments fails validation."""
        self.example["segments"].append({**self.example["segments"][0], "name": name})

        with pytest.raises(
            ValidationError,
            match="All segments must have a name in a multipartite plan.",
        ):
            Plan.model_validate(self.example)

    def test_duplicate_segment_names(self):
        """Test that a plan with duplicate segment names fails validation."""
        self.example["segments"] = [
            {**self.example["segments"][0], "name": {"prefix": "DNA", "key": "A"}}
            for _ in range(2)
        ]

        with pytest.raises(
            ValidationError,
            match="Segment names must be unique within a plan.",
        ):
            Plan.model_validate(self.example)

    def test_get_segment_by_id(self):
        """Test that a segment can be retrieved by its ID."""
        segment = self.example["segments"][0]

        plan = Plan.model_validate(self.example)

        assert (
            plan.get_segment_by_id(segment["id"]).model_dump()
            == Segment.model_validate(segment).model_dump()
        )

    def test_get_segment_by_nonexistent_id(self):
        """Test that None is returned when a segment is not found by its ID."""
        plan = Plan.model_validate(self.example)

        assert plan.get_segment_by_id(uuid4()) is None

    def test_monopartite(self):
        """Test that the monopartite property works as expected."""
        assert Plan.model_validate(self.example).monopartite is True

        self.example["segments"] = [
            {**self.example["segments"][0], "name": {"prefix": "DNA", "key": "A"}}
        ]

        assert Plan.model_validate(self.example).monopartite is True

        self.example["segments"].append(
            {**self.example["segments"][0], "name": {"prefix": "DNA", "key": "B"}}
        )

        assert Plan.model_validate(self.example).monopartite is False

    @pytest.mark.parametrize("rule", ["required", "recommended", "optional"])
    def test_required_segments_monopartite(self, rule: str):
        """Test that the required_segments property works as expected when a monopartite
        plan has either one or no required segments.
        """
        self.example["segments"][0]["required"] = rule

        plan = Plan.model_validate(self.example)

        if rule == "required":
            assert plan.not_required_segments == []
            assert plan.required_segments == plan.segments

        else:
            assert plan.not_required_segments == plan.segments
            assert plan.required_segments == []

    def test_required_segments_multipartite(self):
        """Test that the required_segments property works as expected."""
        self.example["segments"] = [
            {
                **self.example["segments"][0],
                "required": rule,
                "name": {"prefix": "DNA", "key": key},
            }
            for key, rule in zip(
                "ABCD", ["required", "recommended", "required", "optional"], strict=True
            )
        ]

        plan = Plan.model_validate(self.example)

        assert plan.not_required_segments == [plan.segments[1], plan.segments[3]]
        assert plan.required_segments == [plan.segments[0], plan.segments[2]]


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

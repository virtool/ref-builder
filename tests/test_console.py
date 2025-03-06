import pytest
from faker import Faker
from syrupy import SnapshotAssertion

from ref_builder.console import (
    console,
    print_isolate,
    print_isolate_as_json,
    print_otu,
    print_otu_as_json,
    print_otu_list,
)
from ref_builder.models import Molecule, MolType, OTUMinimal, Strandedness, Topology
from ref_builder.plan import Plan, Segment, SegmentName, SegmentRule
from ref_builder.resources import RepoIsolate, RepoOTU, RepoSequence
from ref_builder.utils import Accession, IsolateName, IsolateNameType
from tests.fixtures.factories import OTUMinimalFactory, IsolateFactory
from tests.fixtures.providers import AccessionProvider, SequenceProvider


class TestPrintOTUList:
    """Tests for the ``print_otu_list`` function."""

    def test_ok(
        self,
        capsys: pytest.CaptureFixture,
        otu_minimal_factory: OTUMinimalFactory,
        snapshot: SnapshotAssertion,
    ) -> None:
        """Test that listed OTUs are printed."""
        print_otu_list(otu_minimal_factory.build() for _ in range(5))

        assert capsys.readouterr().out == snapshot

    def test_empty(self) -> None:
        """Test that an empty list of OTUs is printed."""
        with console.capture() as capture:
            print_otu_list(OTUMinimal(**otu) for otu in [])

        assert capture.get() == "No OTUs found\n"


class TestPrintIsolate:
    """Test isolate console output."""

    def test_ok(self, snapshot: SnapshotAssertion):
        """Test that an isolate is printed as expected by ``print_isolate``."""
        fake = Faker(["en_US"])
        fake.add_provider(AccessionProvider)
        fake.add_provider(SequenceProvider)
        fake.seed_instance(8801)

        plan = Plan.new(
            [
                Segment(
                    id=fake.uuid4(),
                    length=1099,
                    length_tolerance=0.03,
                    name=SegmentName("DNA", "R"),
                    rule=SegmentRule.REQUIRED,
                ),
                Segment(
                    id=fake.uuid4(),
                    length=1074,
                    length_tolerance=0.03,
                    name=SegmentName("DNA", "M"),
                    rule=SegmentRule.REQUIRED,
                ),
                Segment(
                    id=fake.uuid4(),
                    length=1087,
                    length_tolerance=0.03,
                    name=SegmentName("DNA", "S"),
                    rule=SegmentRule.REQUIRED,
                ),
            ],
        )

        isolate = RepoIsolate(**IsolateFactory.build_on_plan(plan).model_dump())

        with console.capture() as capture:
            print_isolate(isolate, plan)

        assert capture.get() == snapshot

    def test_json_ok(self, snapshot: SnapshotAssertion):
        """Test that an isolate is printed as expected by ``print_isolate_as_json``."""
        fake = Faker(["en_US"])
        fake.add_provider(AccessionProvider)
        fake.add_provider(SequenceProvider)
        fake.seed_instance(8801)

        plan = Plan.new(
            [
                Segment(
                    id=fake.uuid4(),
                    length=1099,
                    length_tolerance=0.03,
                    name=SegmentName("DNA", "R"),
                    rule=SegmentRule.REQUIRED,
                ),
                Segment(
                    id=fake.uuid4(),
                    length=1074,
                    length_tolerance=0.03,
                    name=SegmentName("DNA", "M"),
                    rule=SegmentRule.REQUIRED,
                ),
                Segment(
                    id=fake.uuid4(),
                    length=1087,
                    length_tolerance=0.03,
                    name=SegmentName("DNA", "S"),
                    rule=SegmentRule.REQUIRED,
                ),
            ],
        )

        isolate = RepoIsolate(**IsolateFactory.build_on_plan(plan).model_dump())

        with console.capture() as capture:
            print_isolate_as_json(isolate)

        assert capture.get() == snapshot


class TestPrintOTU:
    """Test OTU console output."""

    def test_print_otu(self, snapshot: SnapshotAssertion):
        """Test that an OTU is printed as expected by ``print_otu``."""
        fake = Faker(["en_US"])
        fake.add_provider(AccessionProvider)
        fake.add_provider(SequenceProvider)
        fake.seed_instance(8801)

        otu = RepoOTU(
            id=fake.uuid4(),
            acronym="BabAV",
            excluded_accessions=set(),
            legacy_id="",
            molecule=Molecule(
                strandedness=Strandedness.SINGLE,
                topology=Topology.LINEAR,
                type=MolType.RNA,
            ),
            name="Babuvirus abacae",
            plan=Plan.new(
                [
                    Segment(
                        id=fake.uuid4(),
                        length=1099,
                        length_tolerance=0.03,
                        name=SegmentName("DNA", "R"),
                        rule=SegmentRule.REQUIRED,
                    ),
                    Segment(
                        id=fake.uuid4(),
                        length=1074,
                        length_tolerance=0.03,
                        name=SegmentName("DNA", "M"),
                        rule=SegmentRule.REQUIRED,
                    ),
                    Segment(
                        id=fake.uuid4(),
                        length=1087,
                        length_tolerance=0.03,
                        name=SegmentName("DNA", "S"),
                        rule=SegmentRule.REQUIRED,
                    ),
                ],
            ),
            representative_isolate=None,
            taxid=438782,
            isolates=[],
        )

        for _ in range(2):
            sequences = [
                RepoSequence(
                    id=fake.uuid4(),
                    accession=Accession.from_string(fake.accession() + ".1"),
                    definition=fake.sentence(),
                    legacy_id=None,
                    sequence=fake.sequence(),
                    segment=segment.id,
                )
                for segment in otu.plan.segments
            ]

            otu.isolates.append(
                RepoIsolate(
                    id=fake.uuid4(),
                    legacy_id=None,
                    name=IsolateName(type=IsolateNameType.ISOLATE, value=fake.word()),
                    sequences=sequences,
                )
            )

        with console.capture() as capture:
            print_otu(otu)

        assert capture.get() == snapshot

    def test_as_json_ok(self, snapshot: SnapshotAssertion):
        """Test that an OTU is printed as expected by ``print_otu_as_json``."""
        fake = Faker(["en_US"])
        fake.add_provider(AccessionProvider)
        fake.add_provider(SequenceProvider)
        fake.seed_instance(8801)

        otu = RepoOTU(
            id=fake.uuid4(),
            acronym="BabAV",
            excluded_accessions=set(),
            legacy_id="",
            molecule=Molecule(
                strandedness=Strandedness.SINGLE,
                topology=Topology.LINEAR,
                type=MolType.RNA,
            ),
            name="Babuvirus abacae",
            plan=Plan(
                id=fake.uuid4(),
                segments=[
                    Segment(
                        id=fake.uuid4(),
                        length=1099,
                        length_tolerance=0.03,
                        name=SegmentName("DNA", "R"),
                        rule=SegmentRule.REQUIRED,
                    ),
                    Segment(
                        id=fake.uuid4(),
                        length=1074,
                        length_tolerance=0.03,
                        name=SegmentName("DNA", "M"),
                        rule=SegmentRule.REQUIRED,
                    ),
                    Segment(
                        id=fake.uuid4(),
                        length=1087,
                        length_tolerance=0.03,
                        name=SegmentName("DNA", "S"),
                        rule=SegmentRule.REQUIRED,
                    ),
                ],
            ),
            representative_isolate=None,
            taxid=438782,
            isolates=[],
        )

        with console.capture() as capture:
            print_otu_as_json(otu)

        assert capture.get() == snapshot

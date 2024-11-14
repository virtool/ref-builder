import pytest
from faker import Faker
from syrupy import SnapshotAssertion

from ref_builder.console import console, print_otu, print_otu_list
from ref_builder.models import Molecule, MolType, OTUMinimal, Strandedness, Topology
from ref_builder.plan import MultipartitePlan, Segment, SegmentName, SegmentRule
from ref_builder.repo import Repo
from ref_builder.resources import RepoIsolate, RepoOTU, RepoSequence
from ref_builder.utils import Accession, IsolateName, IsolateNameType
from tests.fixtures.factories import OTUMinimalFactory
from tests.fixtures.providers import AccessionProvider, SequenceProvider


class TestPrintOTUList:
    """Tests for the ``print_otu_list`` function."""

    def test_ok(
        self, otu_minimal_factory: OTUMinimalFactory, snapshot: SnapshotAssertion
    ) -> None:
        """Test that listed OTUs are printed."""
        with console.capture() as capture:
            print_otu_list(otu_minimal_factory.build() for _ in range(5))

        assert capture.get() == snapshot

    def test_empty(self) -> None:
        """Test that an empty list of OTUs is printed."""
        with console.capture() as capture:
            print_otu_list(OTUMinimal(**otu) for otu in [])

        assert capture.get() == "No OTUs found\n"


@pytest.mark.parametrize("taxid", [345184, 3158377])
def test_print_otu(scratch_repo: Repo, snapshot: SnapshotAssertion, taxid: int):
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
            strandedness=Strandedness.SINGLE, topology=Topology.LINEAR, type=MolType.RNA
        ),
        name="Babuvirus abacae",
        plan=MultipartitePlan.new(
            [
                Segment(
                    id=fake.uuid4(),
                    length=1099,
                    name=SegmentName("DNA", "R"),
                    required=SegmentRule.REQUIRED,
                ),
                Segment(
                    id=fake.uuid4(),
                    length=1074,
                    name=SegmentName("DNA", "M"),
                    required=SegmentRule.REQUIRED,
                ),
                Segment(
                    id=fake.uuid4(),
                    length=1087,
                    name=SegmentName("DNA", "S"),
                    required=SegmentRule.REQUIRED,
                ),
            ],
        ),
        repr_isolate=None,
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
                segment=str(segment.name),
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

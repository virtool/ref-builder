"""Factories for generating quasi-realistic NCBISource and NCBIGenbank data."""

from faker.providers import lorem
from polyfactory import PostGenerated, Use
from polyfactory.decorators import post_generated
from polyfactory.factories.pydantic_factory import ModelFactory
from pydantic.v1 import UUID4

from ref_builder.models import MolType, OTUMinimal
from ref_builder.ncbi.models import NCBIGenbank, NCBISource, NCBISourceMolType
from ref_builder.otu.models import IsolateBase, OTUBase
from ref_builder.plan import (
    Plan,
    Segment,
    SegmentName,
    SegmentRule,
)
from ref_builder.resources import RepoIsolate, RepoSequence
from ref_builder.utils import Accession, IsolateName, IsolateNameType
from tests.fixtures.providers import (
    AccessionProvider,
    BusinessProvider,
    OrganismProvider,
    SegmentProvider,
    SequenceProvider,
)

DNA_MOLTYPES = {
    NCBISourceMolType.GENOMIC_DNA,
    NCBISourceMolType.OTHER_DNA,
    NCBISourceMolType.UNASSIGNED_DNA,
}
"""NCBISourceMolTypes that map to MolType.DNA"""


ModelFactory.__faker__.add_provider(AccessionProvider)
ModelFactory.__faker__.add_provider(BusinessProvider)
ModelFactory.__faker__.add_provider(OrganismProvider)
ModelFactory.__faker__.add_provider(SegmentProvider)
ModelFactory.__faker__.add_provider(SequenceProvider)
ModelFactory.__faker__.add_provider(lorem)


class NCBISourceFactory(ModelFactory[NCBISource]):
    """NCBISource Factory with quasi-realistic data."""

    host = Use(ModelFactory.__faker__.host)
    """A fake host name for the virus."""

    @classmethod
    def isolate(cls) -> str:
        """Raw isolate name faker."""
        if cls.__faker__.boolean(80):
            delimiter = cls.__faker__.random_element(["", "-", "_"])
            components = cls.__faker__.random_elements(
                [
                    cls.__faker__.country().replace(" ", ""),
                    cls.__faker__.last_name(),
                    cls.__faker__.country_code(),
                    str(cls.__faker__.random_int(0, 9)),
                    str(cls.__faker__.random_int(10, 99)),
                    str(cls.__faker__.random_int(100, 9999)),
                ],
                2,
                unique=True,
            )

            return delimiter.join(components)

        return ""

    organism = Use(ModelFactory.__faker__.organism)
    """A fake name for the virus."""

    taxid = Use(ModelFactory.__faker__.random_int, min=1000, max=999999)

    @classmethod
    def segment(cls) -> str:
        """Raw segment name faker."""
        if cls.__faker__.boolean(80):
            return (
                f"{cls.__faker__.segment_prefix()}"
                f"{cls.__faker__.segment_delimiter()}"
                f"{cls.__faker__.segment_key()}"
            )
        return cls.__faker__.segment_key()

    @classmethod
    def clone(cls) -> str:
        """Raw clone name faker."""
        if cls.__faker__.boolean(10):
            delimiter = cls.__faker__.random_element(["-", "_", " ", "/"])
            return delimiter.join(cls.__faker__.words(2))

        return ""

    @classmethod
    def strain(cls) -> str:
        """Raw strain name faker."""
        if cls.__faker__.boolean(10):
            delimiter = cls.__faker__.random_element(["-", "_", " ", "/"])
            return delimiter.join(cls.__faker__.words(2))

        return ""

    proviral = Use(ModelFactory.__faker__.boolean, chance_of_getting_true=5)
    """Pseudorandom proviral flag for DNA records only."""

    @post_generated
    @classmethod
    def macronuclear(cls, mol_type: NCBISourceMolType) -> bool:
        """Pseudorandom macronuclear flag for DNA records only."""
        if mol_type in DNA_MOLTYPES:
            return cls.__faker__.boolean(5)

        return False

    @post_generated
    @classmethod
    def focus(cls, mol_type: NCBISourceMolType) -> bool:
        """Pseudorandom focus flag for DNA records only.
        Mutually exclusive with transgenic.
        """
        if mol_type in DNA_MOLTYPES:
            return cls.__faker__.boolean(5)

        return False

    @post_generated
    @classmethod
    def transgenic(cls, focus: bool) -> bool:
        """Transgenic flag set to False if focus is True."""
        if focus and cls.__faker__.boolean(5):
            return not focus

        return False


class NCBIGenbankFactory(ModelFactory[NCBIGenbank]):
    """NCBIGenbank Factory with quasi-realistic data."""

    source = Use(NCBISourceFactory.build)

    accession = Use(ModelFactory.__faker__.accession)
    """A fake accession that conforms to NCBI standards."""

    @post_generated
    @classmethod
    def accession_version(cls, accession: str) -> str:
        """Raw accession_version faker."""
        return f"{accession}.{cls.__faker__.random_int(1, 3)}"

    @post_generated
    @classmethod
    def moltype(cls, source: NCBISource) -> MolType:
        """Map moltype field to source.moltype equivalent."""
        if source.mol_type in DNA_MOLTYPES:
            return MolType.DNA

        try:
            return {
                NCBISourceMolType.GENOMIC_RNA: MolType.RNA,
                NCBISourceMolType.MRNA: MolType.MRNA,
                NCBISourceMolType.TRANSCRIBED_RNA: MolType.RNA,
                NCBISourceMolType.VIRAL_CRNA: MolType.CRNA,
                NCBISourceMolType.TRNA: MolType.TRNA,
                NCBISourceMolType.OTHER_RNA: MolType.RNA,
            }[source.mol_type]
        except KeyError as err:
            raise ValueError(
                f"Source moltype {source.mol_type} cannot be matched to MolType",
            ) from err

    @post_generated
    @classmethod
    def organism(cls, source: NCBISource) -> str:
        """Organism faker."""
        return source.organism

    @classmethod
    def sequence(cls) -> str:
        """Sequence faker."""
        return cls.__faker__.sequence()

    @post_generated
    @classmethod
    def taxid(cls, source: NCBISource) -> int:
        """Match taxid field to source.taxid."""
        return source.taxid


def derive_acronym(_: str, values: dict[str, str]) -> str:
    """Derive an acronym from an OTU name."""
    name = values["name"]
    return "".join([part[0].upper() for part in name.split(" ")])


class SegmentFactory(ModelFactory[Segment]):
    """Segment Factory with quasi-realistic data."""

    ModelFactory.__faker__.add_provider(SequenceProvider)

    length = Use(ModelFactory.__faker__.sequence_length)
    """Generate a quasi-realistic length for a sequence."""

    @classmethod
    def name(cls) -> SegmentName | None:
        """Generate a quasi-realistic segment name or null."""
        if cls.__faker__.random_int(0, 10) > 5:
            return SegmentName(
                prefix=cls.__faker__.random_element(["DNA", "RNA"]),
                key=cls.__faker__.segment_key(),
            )

        return None


class SequenceFactory(ModelFactory[RepoSequence]):
    """Sequence factory with quasi-realistic data."""

    id = Use(ModelFactory.__faker__.uuid4)
    """Generate a UUID."""

    definition = Use(ModelFactory.__faker__.sentence)
    """Generate a mock sentence to serve as the definition field."""

    legacy_id = Use(ModelFactory.__faker__.legacy_id)
    """Generate an 8-character unique identifier as used in virtool-cli."""

    segment = Use(ModelFactory.__faker__.uuid4)
    """Generate a quasi-realistic mock segment string."""

    sequence = Use(ModelFactory.__faker__.sequence)
    """Generate a quasi-realistic mock genomic sequence."""

    @classmethod
    def accession(cls) -> Accession:
        """Generate a quasi-realistic accession."""
        ModelFactory.__faker__.add_provider(AccessionProvider)
        return Accession(key=ModelFactory.__faker__.accession(), version=1)


class PlanFactory(ModelFactory[Plan]):
    """A Polyfactory that generates valid instances of :class:`Plan`.

    The factory generates a random number of segments, with a 75% chance of generating a
    monopartite plan.
    """

    @classmethod
    def segments(cls) -> list[Segment]:
        """Return a set of quasi-realistic segments."""
        # The segment represent a monopartite OTU 75% of the time.
        if cls.__faker__.random_int(0, 3):
            return [SegmentFactory.build(name=None, required=SegmentRule.REQUIRED)]

        segment_count = cls.__faker__.random_int(2, 5)
        segment_name_keys = "ABCDEF"

        return [
            SegmentFactory.build(
                name=SegmentName(prefix="DNA", key=segment_name_keys[i]),
                required=SegmentRule.REQUIRED,
            )
            for i in range(segment_count)
        ]


class IsolateFactory(ModelFactory[IsolateBase]):
    """Isolate factory with quasi-realistic data."""

    ModelFactory.__faker__.add_provider(BusinessProvider)

    id = Use(ModelFactory.__faker__.uuid4, cast_to=None)
    """Generate a UUID."""

    legacy_id = Use(ModelFactory.__faker__.legacy_id)
    """Generate an 8-character unique identifier as used in virtool-cli."""

    @classmethod
    def name(cls) -> IsolateName:
        """Generate a quasi-realistic isolate name."""
        return IsolateName(
            type=IsolateNameType.ISOLATE,
            value=cls.__faker__.word(part_of_speech="noun").capitalize(),
        )

    @classmethod
    def sequences(cls) -> list[RepoSequence]:
        """Generate between 1 and 6 sequences with numerically sequential accessions."""
        sequence_count = cls.__faker__.random_int(1, 6)

        return [
            SequenceFactory.build(accession=Accession(key=accession, version=1))
            for accession in cls.__faker__.accessions(sequence_count)
        ]


class OTUFactory(ModelFactory[OTUBase]):
    """OTU Factory with quasi-realistic data."""

    ModelFactory.__faker__.add_provider(BusinessProvider)
    ModelFactory.__faker__.add_provider(OrganismProvider)

    acronym = PostGenerated(derive_acronym)
    """Generate an acronym for the OTU derived from its name."""

    @classmethod
    def excluded_accessions(cls) -> set[str]:
        """Generate a set of excluded accessions."""
        return set()

    @post_generated
    @classmethod
    def isolates(cls, plan: Plan) -> list[IsolateBase]:
        """Derive a list of isolates from a plan."""
        isolates = []

        for _ in range(cls.__faker__.random_int(2, 5)):
            sequences = [
                SequenceFactory.build(segment=segment.id) for segment in plan.segments
            ]

            isolates.append(IsolateFactory.build(sequences=sequences))

        return isolates

    id = Use(ModelFactory.__faker__.uuid4, cast_to=None)
    """Generate a UUID."""

    legacy_id = Use(ModelFactory.__faker__.legacy_id)
    """Generate an 8-character unique identifier as used in virtool-cli."""

    name = Use(ModelFactory.__faker__.organism)
    """Generate a realistic name for a plant virus."""

    plan = Use(PlanFactory.build)
    """Generate a quasi-realistic plan for the OTU."""

    @post_generated
    @classmethod
    def representative_isolate(cls, isolates: list[RepoIsolate]) -> UUID4:
        """Derive a representative isolate from a list of OTUs."""
        return cls.__faker__.random_element(isolates).id

    @post_generated
    @classmethod
    def sequences(cls, isolates: list[RepoIsolate]) -> list[RepoSequence]:
        """Derive a list of sequences from a list of isolates."""
        return [sequence for isolate in isolates for sequence in isolate.sequences]

    taxid = Use(ModelFactory.__faker__.random_int, min=1000, max=999999)
    """A realistic taxonomy ID."""


class OTUMinimalFactory(ModelFactory[OTUMinimal]):
    """OTUMinimal Factory with quasi-realistic data."""

    ModelFactory.__faker__.add_provider(BusinessProvider)
    ModelFactory.__faker__.add_provider(OrganismProvider)

    acronym = PostGenerated(derive_acronym)
    """An acronym for the OTU derived from its name."""

    legacy_id = Use(ModelFactory.__faker__.legacy_id)
    """Generate a realistic 8-character ``legacy_id`` for the OTU."""

    name = Use(ModelFactory.__faker__.organism)
    """Generate a realistic name for the OTU."""

    taxid = Use(ModelFactory.__faker__.random_int, min=1000, max=999999)
    """A realistic taxonomy ID."""

"""Factories for generating quasi-realistic NCBISource and NCBIGenbank data."""

from faker import Faker
from polyfactory import PostGenerated, Use
from polyfactory.decorators import post_generated
from polyfactory.factories.pydantic_factory import ModelFactory
from polyfactory.pytest_plugin import register_fixture
from pydantic.v1 import UUID4

from ref_builder.models import MolType, OTUMinimal
from ref_builder.ncbi.models import NCBIGenbank, NCBISource, NCBISourceMolType
from ref_builder.otu.models import IsolateBase, OTUBase
from ref_builder.plan import (
    MonopartitePlan,
    MultipartitePlan,
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


class NCBISourceFactory(ModelFactory[NCBISource]):
    """NCBISource Factory with quasi-realistic data."""

    __faker__ = Faker()
    __random_seed__ = 10
    __faker__.add_provider(OrganismProvider)
    __faker__.add_provider(SegmentProvider)

    @classmethod
    def taxid(cls) -> int:
        """Taxon ID faker."""
        return cls.__faker__.random_int(1000, 999999)

    @classmethod
    def organism(cls) -> str:
        """Organism name faker."""
        return cls.__faker__.organism()

    @classmethod
    def host(cls) -> str:
        """Pathogen host name faker."""
        return cls.__faker__.host().capitalize()

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

    @classmethod
    def segment(cls) -> str:
        """Raw segment name faker."""
        if cls.__faker__.boolean(80):
            return (
                cls.__faker__.segment_prefix()
                + cls.__faker__.segment_delimiter()
                + cls.__faker__.segment_key()
            )

        return cls.__faker__.segment_key()

    @classmethod
    def clone(cls) -> str:
        """Raw clone name faker."""
        delimiter = cls.__faker__.random_element(["-", "_", " ", "/"])
        if cls.__faker__.boolean(10):
            return delimiter.join(cls.__faker__.words(2))

        return ""

    @classmethod
    def strain(cls) -> str:
        """Raw strain name faker."""
        delimiter = cls.__faker__.random_element(["-", "_", " ", "/"])
        if cls.__faker__.boolean(10):
            return delimiter.join(cls.__faker__.words(2))

        return ""

    @classmethod
    def proviral(cls) -> bool:
        """Pseudorandom proviral flag."""
        return cls.__faker__.boolean(5)

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

    __faker__ = Faker()
    __faker__.add_provider(AccessionProvider)
    __faker__.add_provider(SequenceProvider)

    source = NCBISourceFactory

    @classmethod
    def accession(cls) -> str:
        """Raw accession faker."""
        return cls.__faker__.accession()

    @classmethod
    def sequence(cls) -> str:
        """Sequence faker."""
        return cls.__faker__.sequence()

    @post_generated
    @classmethod
    def accession_version(cls, accession: str) -> str:
        """Raw accession_version faker."""
        return f"{accession}.{cls.__faker__.random_int(1, 3)}"

    @post_generated
    @classmethod
    def organism(cls, source: NCBISource) -> str:
        """Match organism field to source.organism."""
        return source.organism

    @post_generated
    @classmethod
    def taxid(cls, source: NCBISource) -> int:
        """Match taxid field to source.taxid."""
        return source.taxid

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

    @staticmethod
    def build_from_metadata(
        plan: MonopartitePlan | MultipartitePlan,
        moltype: MolType,
        name: str,
        taxid: int,
    ) -> list[NCBIGenbank]:
        """Return a list of mock records matching given OTU metadata."""
        source_moltype = None

        match moltype:
            case MolType.DNA:
                source_moltype = NCBISourceMolType.GENOMIC_DNA
            case MolType.RNA:
                source_moltype = NCBISourceMolType.GENOMIC_RNA
            case MolType.CRNA:
                source_moltype = NCBISourceMolType.VIRAL_CRNA
            case MolType.MRNA:
                source_moltype = NCBISourceMolType.MRNA
            case MolType.TRNA:
                source_moltype = NCBISourceMolType.TRANSCRIBED_RNA

        if source_moltype is None:
            raise ValueError(
                f"MolType {moltype} cannot be matched to NCBISourceMolType"
            )

        if isinstance(plan, MonopartitePlan):
            source = NCBISourceFactory.build(
                moltype=source_moltype,
                organism=name,
                segment="",
                taxid=taxid,
            )

            return [NCBIGenbankFactory.build(source=source)]

        mock_records = []

        for segment in plan.required_segments:
            source = NCBISourceFactory.build(
                moltype=source_moltype,
                organism=name,
                segment=str(segment.name),
                taxid=taxid,
            )

            mock_records.append(NCBIGenbankFactory.build(source=source))

        return mock_records


def derive_acronym(_: str, values: dict[str, str]) -> str:
    """Derive an acronym from an OTU name."""
    name = values["name"]
    return "".join([part[0].upper() for part in name.split(" ")])


class SequenceFactory(ModelFactory[RepoSequence]):
    """Sequence factory with quasi-realistic data."""

    __faker__ = Faker()

    __faker__.add_provider(AccessionProvider)
    __faker__.add_provider(BusinessProvider)
    __faker__.add_provider(SegmentProvider)
    __faker__.add_provider(SequenceProvider)

    id = Use(__faker__.uuid4, cast_to=None)
    """Generate a UUID."""

    definition = Use(__faker__.sentence)
    """Generate a mock sentence to serve as the definition field."""

    legacy_id = Use(__faker__.legacy_id)
    """Generate an 8-character unique identifier as used in virtool-cli."""

    segment = Use(__faker__.segment)
    """Generate a quasi-realistic mock segment string."""

    sequence = Use(__faker__.sequence)
    """Generate a quasi-realistic mock genomic sequence."""

    @classmethod
    def accession(cls) -> Accession:
        """Generate a quasi-realistic accession."""
        return Accession(key=cls.__faker__.accession(), version=1)


class IsolateFactory(ModelFactory[IsolateBase]):
    """Isolate factory with quasi-realistic data."""

    __faker__ = Faker()

    __faker__.add_provider(AccessionProvider)
    __faker__.add_provider(BusinessProvider)

    id = Use(__faker__.uuid4, cast_to=None)
    """Generate a UUID."""

    legacy_id = Use(__faker__.legacy_id)
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
        n_sequences = cls.__faker__.random_int(1, 6)

        accessions = cls.__faker__.accessions(n_sequences)

        sequences = []
        for accession in accessions:
            sequences.append(
                SequenceFactory.build(accession=Accession(key=accession, version=1))
            )

        return sequences


class OTUFactory(ModelFactory[OTUBase]):
    """OTU Factory with quasi-realistic data."""

    __faker__ = Faker()

    __faker__.add_provider(AccessionProvider)
    __faker__.add_provider(BusinessProvider)
    __faker__.add_provider(OrganismProvider)
    __faker__.add_provider(SequenceProvider)

    __random_seed__ = 21

    acronym = PostGenerated(derive_acronym)
    """Generate an acronym for the OTU derived from its name."""

    id = Use(__faker__.uuid4, cast_to=None)
    """Generate a UUID."""

    legacy_id = Use(__faker__.legacy_id)
    """Generate an 8-character unique identifier as used in virtool-cli."""

    name = Use(__faker__.organism)
    """Generate a realistic name for a plant virus."""

    @classmethod
    def excluded_accessions(cls) -> set[str]:
        """Generate a set of excluded accessions."""
        return set()

    @classmethod
    def taxid(cls) -> int:
        """Generate a realistic taxonomy ID."""
        return cls.__faker__.random_int(1000, 999999)

    @classmethod
    def plan(cls) -> MultipartitePlan:
        """Generate a multipartite plan with two segments."""
        return MultipartitePlan.new(
            [
                Segment(
                    id=cls.__faker__.uuid4(),
                    length=1099,
                    length_tolerance=0.03,
                    name=SegmentName("DNA", "A"),
                    required=SegmentRule.REQUIRED,
                ),
                Segment(
                    id=cls.__faker__.uuid4(),
                    length=1074,
                    length_tolerance=0.03,
                    name=SegmentName("DNA", "B"),
                    required=SegmentRule.REQUIRED,
                ),
            ],
        )

    @post_generated
    @classmethod
    def isolates(cls, plan: MultipartitePlan) -> list[IsolateBase]:
        """Derive a list of isolates from a plan."""
        isolates = []

        for _ in range(cls.__faker__.random_int(2, 5)):
            sequences = [
                SequenceFactory.build(segment=str(segment.name))
                for segment in plan.segments
            ]

            isolates.append(IsolateFactory.build(sequences=sequences))

        return isolates

    @post_generated
    @classmethod
    def repr_isolate(cls, isolates: list[RepoIsolate]) -> UUID4:
        """Derive a representative isolate from a list of OTUs."""
        return cls.__faker__.random_element(isolates).id

    @post_generated
    @classmethod
    def sequences(cls, isolates: list[RepoIsolate]) -> list[RepoSequence]:
        """Derive a list of sequences from a list of isolates."""
        return [sequence for isolate in isolates for sequence in isolate.sequences]


@register_fixture
class OTUMinimalFactory(ModelFactory[OTUMinimal]):
    """OTUMinimal Factory with quasi-realistic data."""

    __faker__ = Faker()
    __faker__.add_provider(BusinessProvider)
    __faker__.add_provider(OrganismProvider)

    __random_seed__ = 20

    acronym = PostGenerated(derive_acronym)
    """An acronym for the OTU derived from its name."""

    @classmethod
    def legacy_id(cls) -> str:
        """Generate a realistic 8-character ``legacy_id`` for the OTU."""
        return cls.__faker__.legacy_id()

    @classmethod
    def name(cls) -> str:
        """Generate a realistic name for the OTU."""
        return cls.__faker__.organism()

    @classmethod
    def taxid(cls) -> int:
        """Generate a realistic taxonomy ID."""
        return cls.__faker__.random_int(1000, 999999)

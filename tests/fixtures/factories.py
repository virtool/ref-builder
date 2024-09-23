from faker import Faker
from polyfactory.decorators import post_generated
from polyfactory.factories.pydantic_factory import ModelFactory

from ref_builder.models import MolType
from ref_builder.ncbi.models import NCBIGenbank, NCBISource, NCBISourceMolType
from tests.fixtures.providers import (
    AccessionProvider,
    OrganismProvider,
    SegmentProvider,
    SequenceProvider,
)

DNA_MOLTYPES = {
    NCBISourceMolType.GENOMIC_DNA,
    NCBISourceMolType.OTHER_DNA,
    NCBISourceMolType.UNASSIGNED_DNA
}


class NCBISourceFactory(ModelFactory[NCBISource]):
    """NCBISource Factory with quasi-realistic data."""

    __faker__ = Faker()
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
            components = cls.__faker__.random_elements([
                cls.__faker__.country().replace(" ", ""),
                cls.__faker__.last_name(),
                cls.__faker__.country_code(),
                str(cls.__faker__.random_int(0, 9)),
                str(cls.__faker__.random_int(10, 99)),
                str(cls.__faker__.random_int(100, 9999)),
            ], 2, unique=True)

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

        match source.mol_type:
            case NCBISourceMolType.GENOMIC_RNA:
                return MolType.RNA
            case NCBISourceMolType.MRNA:
                return MolType.MRNA
            case NCBISourceMolType.TRANSCRIBED_RNA:
                return MolType.RNA
            case NCBISourceMolType.VIRAL_CRNA:
                return MolType.CRNA
            case NCBISourceMolType.TRNA:
                return MolType.TRNA
            case NCBISourceMolType.OTHER_RNA:
                return MolType.RNA

        raise ValueError(f"Source moltype {source.mol_type} cannot be matched to MolType")

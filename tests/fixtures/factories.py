from faker import Faker
from polyfactory import Use
from polyfactory.decorators import post_generated
from polyfactory.factories.pydantic_factory import ModelFactory

from ref_builder.models import MolType
from ref_builder.utils import Accession
from ref_builder.ncbi.models import NCBIGenbank, NCBISource, NCBISourceMolType
from tests.fixtures.providers import (
    AccessionProvider, OrganismProvider, SequenceProvider, SourceProvider
)

DNA_MOLTYPES = {
    NCBISourceMolType.GENOMIC_DNA,
    NCBISourceMolType.OTHER_DNA,
    NCBISourceMolType.UNASSIGNED_DNA
}


class NCBISourceFactory(ModelFactory[NCBISource]):
    __faker__ = Faker()
    __faker__.add_provider(OrganismProvider)
    __faker__.add_provider(SourceProvider)

    @classmethod
    def taxid(cls) -> int:
        return cls.__faker__.random_int(1000, 999999)

    @classmethod
    def organism(cls) -> str:
        return cls.__faker__.organism()

    @classmethod
    def host(cls) -> str:
        return cls.__faker__.host().capitalize()

    @classmethod
    def isolate(cls) -> str:
        if cls.__faker__.boolean(80):
            return cls.__faker__.text(20).replace(".", "")
        return ""

    @classmethod
    def segment(cls):
        if cls.__faker__.boolean(80):
            return (
                cls.__faker__.segment_prefix()
                + cls.__faker__.segment_delimiter()
                + cls.__faker__.segment_key()
            )

        return cls.__faker__.segment_key()

    @classmethod
    def clone(cls) -> str:
        if cls.__faker__.boolean(10):
            return cls.__faker__.text(20).replace(".", "")

        return ""

    @classmethod
    def strain(cls):
        if cls.__faker__.boolean(20):
            return cls.__faker__.text(20).replace(".", "")

        return ""

    @classmethod
    def proviral(cls) -> bool:
        return cls.__faker__.boolean(5)

    @post_generated
    @classmethod
    def macronuclear(cls, mol_type: NCBISourceMolType) -> bool:
        if mol_type in DNA_MOLTYPES:
            return cls.__faker__.boolean(5)
        return False

    @post_generated
    @classmethod
    def focus(cls, mol_type: NCBISourceMolType) -> bool:
        """Mutually exclusive with transgenic"""
        if mol_type in DNA_MOLTYPES:
            return cls.__faker__.boolean(5)
        return False

    @post_generated
    @classmethod
    def transgenic(cls, focus: bool) -> bool:
        """Mutually exclusive with focus"""
        if focus:
            return not focus
        return False


class NCBIGenbankFactory(ModelFactory[NCBIGenbank]):
    __faker__ = Faker()
    __faker__.add_provider(AccessionProvider)
    __faker__.add_provider(SequenceProvider)

    source = NCBISourceFactory

    @classmethod
    def accession(cls) -> str:
        return cls.__faker__.accession()

    @classmethod
    def sequence(cls):
        return cls.__faker__.sequence()

    @post_generated
    @classmethod
    def accession_versioned(cls, accession: str):
        return f"{accession}.{cls.__faker__.random_int(1, 3)}"

    @post_generated
    @classmethod
    def organism(cls, source: NCBISource) -> str:
        return source.organism

    @post_generated
    @classmethod
    def taxid(cls, source: NCBISource) -> int:
        return source.taxid

    @post_generated
    @classmethod
    def moltype(self, source: NCBISource) -> MolType:
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


if __name__ == '__main__':
    records = NCBIGenbankFactory.batch(10)

    for record in records:
        print(record)

    # sources = NCBISourceFactory.batch(10)
    #
    # for source in sources:
    #     print(source)
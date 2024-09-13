from faker import Faker
from polyfactory import Use
from polyfactory.decorators import post_generated
from polyfactory.factories.pydantic_factory import ModelFactory

from ref_builder.ncbi.models import NCBIGenbank, NCBISource, NCBISourceMolType
from ref_builder.utils import Accession
from tests.fixtures.providers import OrganismProvider, SourceProvider

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


if __name__ == '__main__':
    sources = NCBISourceFactory.batch(10)

    for source in sources:
        print(source)
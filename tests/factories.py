from typing import Dict

from faker import Faker
from pydantic import BaseModel, Field
from polyfactory import Use, PostGenerated
from polyfactory.factories import DataclassFactory
from polyfactory.factories.pydantic_factory import ModelFactory
from typing import Annotated

from ref_builder.ncbi.models import NCBISource
from ref_builder.schema import Segment, OTUSchema
from ref_builder.snapshotter.models import OTUSnapshotOTU, OTUSnapshotIsolate, OTUSnapshotSequence
from ref_builder.utils import IsolateName, IsolateNameType
from ref_builder.models import Strandedness, MolType, Topology


SEQUENCE_CHARACTERS = 'ATCGRYKMSWBDHVN'
SEQUENCE_MIN = 100
SEQUENCE_MAX = 1500


def generate_sequence(generator) -> str:
    return ''.join(generator.choices(
        population=SEQUENCE_CHARACTERS,
        weights=[50]*4 + [1]*(len(SEQUENCE_CHARACTERS)-4),
        k=generator.randint(SEQUENCE_MIN, SEQUENCE_MAX),
    ))


def generate_legacy_id(generator) -> str | None:
    if generator.boolean(50):
        return generator.password(
            length=8,
            special_chars=False,
            digits=True,
            upper_case=True,
            lower_case=False,
        )

    return None


def construct_versioned_accession(
    name: str, values: Dict[str, str | int], version: int
) -> str:
    if name == "accession_version":
        return f"{values["accession"]}.{version}"

    raise ValueError("Field name is not accession_version")


class SegmentFactory(ModelFactory[Segment]):
    __faker__ = Faker(locale="la")

    @classmethod
    def name(cls) -> str:
        return cls.__faker__.text(10).replace(".", "")

    @classmethod
    def length(cls):
        return cls.__random__.randint(SEQUENCE_MIN, SEQUENCE_MAX)


class OTUSchemaFactory(ModelFactory[OTUSchema]):
    @classmethod
    def segments(cls) -> list[Segment]:
        return SegmentFactory.batch(size=cls.__random__.randint(1, 6))


class OTUFactory(ModelFactory[OTUSnapshotOTU]):
    __faker__ = Faker(locale="la")

    schema = OTUSchemaFactory

    @classmethod
    def acronym(cls) -> str:
        if cls.__faker__.boolean(70):
            return cls.__faker__.lexify(text="????")

        return ""

    @classmethod
    def legacy_id(cls) -> str | None:
        return generate_legacy_id(cls.__faker__)

    @classmethod
    def name(cls) -> str:
        return cls.__faker__.text(20).replace(".", "")

    @classmethod
    def taxid(cls):
        return cls.__random__.randint(10000, 9999999)


class IsolateNameFactory(DataclassFactory[IsolateName]):
    __faker__ = Faker(locale="la")

    @classmethod
    def type(cls):
        return cls.__random__.choice(list(IsolateNameType))

    @classmethod
    def value(cls) -> str:
        return cls.__faker__.text(15).replace(".", "")


class IsolateFactory(ModelFactory[OTUSnapshotIsolate]):
    name = IsolateNameFactory

    @classmethod
    def legacy_id(cls) -> str | None:
        return generate_legacy_id(cls.__faker__)


class SequenceFactory(ModelFactory[OTUSnapshotSequence]):
    __faker__ = Faker(locale="la")

    sequence = Use(generate_sequence, generator=__faker__.random)
    accession_version = PostGenerated(
        construct_versioned_accession, version=__faker__.random_element(elements=(1, 2))
    )

    @classmethod
    def accession(cls) -> str:
        gen_length = 7 if cls.__faker__.boolean(50) else 9

        return cls.__faker__.password(
            length=gen_length,
            special_chars=False,
            digits=True,
            upper_case=True,
            lower_case=False
        )

    @classmethod
    def definition(cls):
        return cls.__faker__.text(50).replace(".", "")

    @classmethod
    def segment(cls) -> str:
        return cls.__faker__.text(10).replace(".", "")

    @classmethod
    def legacy_id(cls) -> str | None:
        return generate_legacy_id(cls.__faker__)


class NCBISourceFactory(ModelFactory[NCBISource]):
    __faker__ = Faker(locale="la")

    @classmethod
    def segment(cls) -> str:
        return cls.__faker__.text(10).replace(".", "")

    @classmethod
    def organism(cls):
        return cls.__faker__.text(20).replace(".", "")

    @classmethod
    def host(cls):
        return cls.__faker__.text(20).replace(".", "")

    @classmethod
    def isolate(cls) -> str:
        if cls.__faker__.boolean(80):
            return cls.__faker__.text(20).replace(".", "")

        return ""

    @classmethod
    def clone(cls) -> str:
        if cls.__faker__.boolean(10):
            return cls.__faker__.text(20).replace(".", "")

        return ""

    @classmethod
    def strain(cls):
        if cls.__faker__.boolean(40):
            return cls.__faker__.text(20).replace(".", "")

        return ""

    @classmethod
    def proviral(cls) -> bool:
        return cls.__faker__.boolean(5)

    @classmethod
    def macronuclear(cls) -> bool:
        return cls.__faker__.boolean(5)

    @classmethod
    def focus(cls) -> bool:
        return cls.__faker__.boolean(5)

    @classmethod
    def transgenic(cls) -> bool:
        return cls.__faker__.pybool(5)


class NCBIRecord(BaseModel):
    accession: str
    accession_version: str
    strandedness: Strandedness
    moltype: MolType
    topology: Topology
    definition: str
    organism: str
    sequence: Annotated[
        str,
        Field(pattern=r"^[ATCGRYKMSWBDHVNatcgrykmswbdhvn]+$"),
    ]
    source: NCBISource
    comment: str = ""


class NCBIRecordFactory(ModelFactory[NCBIRecord]):
    __faker__ = Faker(locale="la")

    sequence = Use(generate_sequence, generator=__faker__.random)
    source = NCBISourceFactory
    accession_version = PostGenerated(
        construct_versioned_accession, version=__faker__.random_element(elements=(1, 2))
    )

    @classmethod
    def accession(cls) -> str:
        gen_length = 7 if cls.__faker__.boolean(50) else 9

        return cls.__faker__.password(
            length=gen_length,
            special_chars=False,
            digits=True,
            upper_case=True,
            lower_case=False
        )

    @classmethod
    def definition(cls):
        return cls.__faker__.text(50).replace(".", "")

    @classmethod
    def organism(cls):
        return cls.__faker__.text(20).replace(".", "")

    @classmethod
    def host(cls):
        return cls.__faker__.text(20).replace(".", "")

    @classmethod
    def comment(cls):
        if cls.__faker__.boolean(10):
            return cls.__faker__.text(100)

        return ""


if __name__ == '__main__':
    otu_instance = OTUFactory.build()
    print(otu_instance)

    source_instance = NCBISourceFactory.build()
    print(source_instance)

    record_instance = NCBIRecordFactory.build()
    print(record_instance)

    isolate_instance = IsolateFactory.build()
    print(isolate_instance)

    sequence_instance = SequenceFactory.build()
    print(sequence_instance)

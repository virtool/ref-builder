from pydantic import BaseModel, Field
from polyfactory import Use
from polyfactory.factories import DataclassFactory
from polyfactory.factories.pydantic_factory import ModelFactory
from typing import Annotated

from ref_builder.ncbi.models import NCBIGenbank, NCBISource
from ref_builder.schema import Segment, OTUSchema
from ref_builder.snapshotter.models import OTUSnapshotOTU, OTUSnapshotIsolate, OTUSnapshotSequence
from ref_builder.utils import IsolateName, IsolateNameType
from ref_builder.models import Strandedness, MolType, Topology


SEQUENCE_CHARACTERS = 'ATCGRYKMSWBDHVN'


def generate_sequence() -> str:
    return ''.join(DataclassFactory.__random__.choices(
        population=SEQUENCE_CHARACTERS,
        weights=[50]*4 + [1]*(len(SEQUENCE_CHARACTERS)-4),
        k=DataclassFactory.__random__.randint(100, 1500),
    ))


def generate_bool(p: float = 0.05) -> bool:
    dice = DataclassFactory.__random__.random()
    return True if dice < p else False


class SegmentFactory(ModelFactory[Segment]):
    @classmethod
    def length(cls):
        return cls.__random__.randint(100, 1500)


class OTUSchemaFactory(ModelFactory[OTUSchema]):
    segments = Use(SegmentFactory.batch, size=ModelFactory.__random__.randint(1, 6))


class OTUFactory(ModelFactory[OTUSnapshotOTU]):
    schema = OTUSchemaFactory

    @classmethod
    def taxid(cls):
        return cls.__random__.randint(10000, 9999999)


class IsolateNameFactory(DataclassFactory[IsolateName]):
    @classmethod
    def type(cls):
        return cls.__random__.choice(list(IsolateNameType))


class IsolateFactory(ModelFactory[OTUSnapshotIsolate]):
    name = IsolateNameFactory


class SequenceFactory(ModelFactory[OTUSnapshotSequence]):
    sequence = Use(generate_sequence)


class NCBISourceFactory(ModelFactory[NCBISource]):
    proviral = Use(generate_bool, 0.05)
    macronuclear = Use(generate_bool, 0.05)
    focus = Use(generate_bool, 0.05)
    transgenic = Use(generate_bool,0.05)


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
    sequence = Use(generate_sequence)
    source = NCBISourceFactory



if __name__ == '__main__':
    # otu_instance = OTUFactory.build()
    #
    # print(otu_instance)

    # source_instance = NCBISourceFactory.build()
    #
    # print(source_instance)

    record_instance = NCBIRecordFactory.build()
    print(record_instance)

    # isolate_instance = IsolateFactory.coverage()
    # print(isolate_instance)

    # sequence_instance = SequenceFactory.build()
    # print(sequence_instance)



from polyfactory import Use
from polyfactory.factories import DataclassFactory
from polyfactory.factories.pydantic_factory import ModelFactory

from ref_builder.schema import Segment, OTUSchema
from ref_builder.snapshotter.models import OTUSnapshotOTU, OTUSnapshotIsolate, OTUSnapshotSequence
from ref_builder.utils import IsolateName, IsolateNameType


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
    # __model__ = OTUSnapshotSequence
    @classmethod
    def sequence(cls):
        return ''.join(
            cls.__random__.choices(
                ['G', 'A', 'C', 'T', 'Y'],
                weights=[5, 5, 5, 5, 0.5],
                k=cls.__random__.randint(100, 1500),
            )
        )


if __name__ == '__main__':
    otu_instance = OTUFactory.build()


    print(otu_instance)

    # isolate_instance = IsolateFactory.coverage()
    #
    # sequence_instance = SequenceFactory.coverage()
    # print(isolate_instance)
    #
    # print(sequence_instance)


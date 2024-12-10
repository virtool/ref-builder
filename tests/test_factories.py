"""Tests for data factories used in testing."""

from uuid import uuid4

from ref_builder.ncbi.models import NCBISource
from ref_builder.resources import RepoIsolate, RepoOTU, RepoSequence
from ref_builder.utils import Accession
from tests.fixtures.factories import (
    IsolateFactory,
    SequenceFactory,
    OTUFactory,
    NCBIGenbankFactory,
    NCBISourceFactory,
)


class TestNCBIFactories:
    def test_source_factory(self):
        """Test that NCBISourceFactory creates valid mock NCBISource objects."""
        assert all(
            [
                isinstance(dummy_source, NCBISource)
                for dummy_source in NCBISourceFactory.coverage()
            ]
        )

    def test_genbank_factory(self):
        """Test that NCBIGenbankFactory creates mock NCBIGenbank records
        and ensure said records are valid RepoSequence input."""

        for dummy_record in NCBIGenbankFactory.coverage():
            assert RepoSequence(
                id=uuid4(),
                accession=Accession.from_string(dummy_record.accession_version),
                definition=dummy_record.definition,
                legacy_id=None,
                segment=dummy_record.source.segment,
                sequence=dummy_record.sequence,
            )


def test_sequence_factory():
    """Test that SequenceFactory creates valid mock sequence data."""
    for dummy_sequence in SequenceFactory.coverage():
        assert RepoSequence.model_validate(dummy_sequence.model_dump())


def test_isolate_factory():
    """Test that IsolateFactory creates valid mock isolate data."""
    for dummy_isolate in IsolateFactory.coverage():
        assert RepoIsolate.model_validate(dummy_isolate.model_dump())


def test_otu_factory():
    """Test that OTUFactory creates valid mock OTU data."""
    for dummy_otu in OTUFactory.coverage():
        assert RepoOTU.model_validate(dummy_otu.model_dump())

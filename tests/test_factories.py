"""Tests for data factories used in testing."""

from uuid import uuid4

import pytest

from ref_builder.ncbi.models import NCBIGenbank, NCBISource
from ref_builder.resources import RepoIsolate, RepoOTU, RepoSequence
from ref_builder.otu.models import IsolateBase, OTUBase
from ref_builder.utils import Accession
from tests.fixtures.factories import (
    IsolateFactory,
    SequenceFactory,
    OTUFactory,
    NCBIGenbankFactory,
    NCBISourceFactory,
)


class TestNCBIFactories:
    @pytest.mark.parametrize("dummy_source", NCBISourceFactory.coverage(10))
    def test_source_factory(self, dummy_source: NCBISource):
        """Test that NCBISourceFactory creates valid input."""
        assert dummy_source

    @pytest.mark.parametrize("dummy_record", NCBIGenbankFactory.coverage(10))
    def test_genbank_factory(self, dummy_record: NCBIGenbank):
        """Test that NCBISourceFactory creates valid input."""
        assert dummy_record

        assert RepoSequence(
            id=uuid4(),
            accession=Accession.from_string(dummy_record.accession_version),
            definition=dummy_record.definition,
            legacy_id=None,
            segment=dummy_record.source.segment,
            sequence=dummy_record.sequence,
        )


@pytest.mark.parametrize("dummy_sequence", SequenceFactory.coverage(10))
def test_sequence_factory(dummy_sequence: RepoSequence):
    assert dummy_sequence

    print(dummy_sequence)

    assert RepoSequence.model_validate(dummy_sequence.model_dump())


@pytest.mark.parametrize("dummy_isolate", IsolateFactory.coverage(5))
def test_isolate_factory(dummy_isolate: IsolateBase):
    assert dummy_isolate

    assert RepoIsolate.model_validate(dummy_isolate.model_dump())


@pytest.mark.parametrize("dummy_otu", OTUFactory.coverage(10))
def test_otu_factory(dummy_otu: OTUBase):
    assert dummy_otu

    assert RepoOTU.model_validate(dummy_otu.model_dump())

"""Tests for data factories used in testing."""

from uuid import uuid4

import pytest

from ref_builder.ncbi.models import NCBIGenbank, NCBISource
from ref_builder.resources import RepoSequence
from ref_builder.utils import Accession
from tests.fixtures.factories import NCBIGenbankFactory, NCBISourceFactory


@pytest.mark.parametrize("dummy_source", NCBISourceFactory.coverage(10))
def test_source_factory(dummy_source: NCBISource):
    """Test that NCBISourceFactory creates valid input."""
    assert dummy_source


@pytest.mark.parametrize("dummy_record", NCBIGenbankFactory.coverage(10))
def test_genbank_factory(dummy_record: NCBIGenbank):
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

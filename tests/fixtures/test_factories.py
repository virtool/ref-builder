from uuid import uuid4

import pytest
import structlog
from factories import NCBIGenbankFactory, NCBISourceFactory

from ref_builder.ncbi.models import NCBISource
from ref_builder.resources import RepoIsolate, RepoSequence
from ref_builder.utils import Accession

test_logger = structlog.get_logger("test")


@pytest.mark.parametrize("dummy_source", NCBISourceFactory.coverage(10))
def test_source_factory(dummy_source: NCBISource):
    """Test that NCBISourceFactory creates valid input."""
    assert dummy_source

    test_logger.info("Dummy source data generated.", source=dummy_source.model_dump())


@pytest.mark.parametrize("dummy_record", NCBIGenbankFactory.coverage(10))
def test_genbank_factory(dummy_record: NCBIGenbankFactory):
    """Test that NCBISourceFactory creates valid input."""
    assert dummy_record

    test_logger.info("Dummy record data generated.", record=dummy_record.model_dump())

    sequence = RepoSequence(
        id=uuid4(),
        accession=Accession.from_string(dummy_record.accession_version),
        definition=dummy_record.definition,
        legacy_id=None,
        segment=dummy_record.source.segment,
        sequence=dummy_record.sequence,
    )

    assert sequence

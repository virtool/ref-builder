"""Tests for data factories used in testing."""

import uuid
from uuid import uuid4

from syrupy import SnapshotAssertion

from ref_builder.ncbi.models import NCBISource
from ref_builder.resources import RepoIsolate, RepoOTU, RepoSequence
from ref_builder.utils import Accession
from tests.fixtures.factories import (
    IsolateFactory,
    NCBIGenbankFactory,
    NCBISourceFactory,
    OTUFactory,
    SequenceFactory,
)


def test_ncbi_source_factory(
    ncbi_source_factory: NCBISourceFactory, snapshot: SnapshotAssertion
):
    """Test that NCBISourceFactory creates valid mock NCBISource objects."""
    assert all(
        isinstance(dummy_source, NCBISource)
        for dummy_source in ncbi_source_factory.coverage()
    )


def test_ncbi_genbank_factory(
    ncbi_genbank_factory: NCBIGenbankFactory, snapshot: SnapshotAssertion
):
    """Test that NCBIGenbankFactory creates valid fake Genbank records."""
    records = list(ncbi_genbank_factory.coverage())

    assert [
        RepoSequence(
            id=uuid4(),
            accession=Accession.from_string(record.accession_version),
            definition=record.definition,
            legacy_id=None,
            segment=uuid.uuid4(),
            sequence=record.sequence,
        )
        for record in records
    ]

    assert [record.model_dump() for record in records] == snapshot


def test_sequence_factory(sequence_factory: SequenceFactory):
    """Test that SequenceFactory creates valid mock sequence data."""
    assert all(
        RepoSequence.model_validate(sequence.model_dump())
        for sequence in sequence_factory.coverage()
    )


def test_isolate_factory():
    """Test that IsolateFactory creates valid mock isolate data."""
    assert all(
        RepoIsolate.model_validate(isolate.model_dump())
        for isolate in IsolateFactory.coverage()
    )


def test_otu_factory(otu_factory: OTUFactory):
    """Test that OTUFactory creates valid mock OTU data."""
    assert all(
        RepoOTU.model_validate(otu.model_dump()) for otu in otu_factory.coverage()
    )

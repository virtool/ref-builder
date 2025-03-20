"""Tests for :class:`.EventIndex`."""

import uuid
from pathlib import Path

import arrow
import pytest

from ref_builder.index import EventIndexItem, Index
from ref_builder.models import OTUMinimal
from ref_builder.resources import RepoOTU

SNAPSHOT_AT_EVENT = (
    31,
    22,
    21,
    24,
    11,
    32,
    32,
    30,
)
"""Hardcoded at_event values for the snapshotter_otus fixture."""


@pytest.fixture()
def index(indexable_otus: list[RepoOTU], tmp_path: Path) -> Index:
    """An index with eight OTUs already cached."""
    _index = Index(tmp_path / "index.db")

    for otu, at_event in zip(indexable_otus, SNAPSHOT_AT_EVENT, strict=True):
        _index.upsert_otu(otu, at_event)

    return _index


class TestDeleteOTU:
    """Test the `delete_otu` method of the Snapshotter class."""

    def test_ok(self, index: Index, indexable_otus: list[RepoOTU]):
        """Test that an OTU can be deleted from the index."""
        otu = indexable_otus[2]

        index.delete_otu(otu.id)

        assert index.get_event_ids_by_otu_id(otu.id) is None
        assert index.get_id_by_name(otu.name) is None
        assert index.get_id_by_taxid(otu.taxid) is None

    def test_not_found(self, index: Index, indexable_otus: list[RepoOTU]):
        """Test that nothing happens when the OTU ID is not found."""
        index.delete_otu(uuid.uuid4())
        assert index.otu_ids == {otu.id for otu in indexable_otus}


def test_iter_otus(index: Index, indexable_otus: list[RepoOTU]):
    """Test that the index iterates over all OTUs ordered by name."""
    assert list(index.iter_minimal_otus()) == sorted(
        [
            OTUMinimal(
                acronym=otu.acronym,
                id=otu.id,
                legacy_id=otu.legacy_id,
                name=otu.name,
                taxid=otu.taxid,
            )
            for otu in indexable_otus
        ],
        key=lambda otu: otu.name,
    )


class TestLoadSnapshot:
    """Test the `load_snapshot` method of the Snapshotter class."""

    def test_ok(self, index: Index, indexable_otus: list[RepoOTU]):
        """Test that we can load a snapshot from the index."""
        for at_event, otu in zip(SNAPSHOT_AT_EVENT, indexable_otus, strict=False):
            snapshot = index.load_snapshot(otu.id)

            assert snapshot.at_event == at_event
            assert snapshot.otu == otu

    def test_not_found(self, index: Index):
        """Test that `None` is returned when the OTU ID is not found."""
        assert index.load_snapshot(uuid.uuid4()) is None


def test_otu_ids(index: Index, indexable_otus: list[RepoOTU]):
    """Test that stored OTU IDs are correct."""
    assert index.otu_ids == {otu.id for otu in indexable_otus}


class TestGetIDByLegacyID:
    """Test the `get_id_by_legacy_id` method of the Snapshotter class."""

    def test_ok(self, index: Index, indexable_otus: list[RepoOTU]):
        """Test that the correct OTU ID is retrieved by legacy ID."""
        for otu in indexable_otus:
            if otu.legacy_id is not None:
                assert index.get_id_by_legacy_id(otu.legacy_id) == otu.id

    def test_not_found(self, index: Index):
        """Test that `None` is returned when the legacy ID is not found."""
        assert index.get_id_by_legacy_id("not found") is None


class TestEvents:
    """Test the event index functionality of the repository Index."""

    def test_ok(self, index: Index, indexable_otus: list[RepoOTU]):
        """Test that we can set and get events IDs for an OTU."""
        otu = indexable_otus[1]

        index.add_event_id(100, otu.id, arrow.utcnow().naive)
        index.add_event_id(101, otu.id, arrow.utcnow().naive)
        index.add_event_id(104, otu.id, arrow.utcnow().naive)

        assert index.get_event_ids_by_otu_id(otu.id) == EventIndexItem(
            event_ids=[100, 101, 104],
            otu_id=otu.id,
        )

    def test_otu_id_not_found(self, index: Index):
        """Test that we get ``None`` when an OTU ID is not found."""
        assert index.get_event_ids_by_otu_id(uuid.uuid4()) is None

    def test_get_first_timestamp_ok(self, index: Index, indexable_otus: list[RepoOTU]):
        """Test ``.get_first_timestamp_by_otu_id()`` retrieves the first timestamp."""
        otu = indexable_otus[1]

        first_timestamp = arrow.utcnow().naive

        index.add_event_id(100, otu.id, first_timestamp)

        assert index.get_first_timestamp_by_otu_id(otu.id) == first_timestamp

        second_timestamp = arrow.utcnow().naive

        index.add_event_id(101, otu.id, second_timestamp)

        assert index.get_first_timestamp_by_otu_id(otu.id) == first_timestamp

    def test_get_latest_timestamp_ok(self, index: Index, indexable_otus: list[RepoOTU]):
        """Test ``.get_latest_timestamp_by_otu_id()`` retrieves the latest timestamp."""
        otu = indexable_otus[1]

        index.add_event_id(100, otu.id, arrow.utcnow().naive)

        first_timestamp = arrow.utcnow().naive

        index.add_event_id(101, otu.id, first_timestamp)

        assert index.get_latest_timestamp_by_otu_id(otu.id) == first_timestamp

        second_timestamp = arrow.utcnow().naive

        index.add_event_id(104, otu.id, second_timestamp)

        assert second_timestamp > first_timestamp

        assert index.get_latest_timestamp_by_otu_id(otu.id) == second_timestamp


class TestGetIDByPartial:
    """Test the `get_id_by_partial` method of the Snapshotter class."""

    def test_ok(self, index: Index, indexable_otus: list[RepoOTU]):
        """Test that the correct OTU ID is retrieved by a truncated partial."""
        for otu in indexable_otus:
            assert otu.id == index.get_id_by_partial(str(otu.id)[:8])

    def test_not_found(self, index: Index):
        """Test that `None` is returned when no matching ID is not found."""
        assert index.get_id_by_partial("00000000") is None


class TestGetIDByName:
    """Test the `get_id_by_name` method of the Snapshotter class."""

    def test_ok(self, index: Index, indexable_otus: list[RepoOTU]):
        """Test that the correct OTU ID is retrieved by name."""
        for otu in indexable_otus:
            assert otu.id == index.get_id_by_name(otu.name)

    def test_not_found(self, index: Index):
        """Test that `None` is returned when the name is not found."""
        assert index.get_id_by_name("not found") is None


class TestGetIDByTaxid:
    """Test the `get_id_by_taxid` method of the Snapshotter class."""

    def test_ok(self, index: Index, indexable_otus: list[RepoOTU]):
        """Test that the correct OTU ID is retrieved by taxid."""
        for otu in indexable_otus:
            assert index.get_id_by_taxid(otu.taxid) == otu.id

    def test_not_found(self, index: Index):
        """Test that `None` is returned when the taxid is not found."""
        assert index.get_id_by_taxid(999999999999999) is None


def test_get_id_by_isolate_id(index: Index, indexable_otus: list[RepoOTU]):
    """Test the `get_id_by_isolate_id` method of the Index class."""
    for otu in indexable_otus:
        first_isolate = next(iter(otu.isolate_ids))

        assert index.get_id_by_isolate_id(first_isolate) == otu.id


def test_get_isolate_id_by_partial(index: Index, indexable_otus: list[RepoOTU]):
    """Test the `get_isolate_id_by_partial` method of the Index class."""
    for otu in indexable_otus:
        first_isolate_id = next(iter(otu.isolate_ids))

        assert (
            index.get_isolate_id_by_partial(str(first_isolate_id)[:8])
            == first_isolate_id
        )

"""Tests for the Snapshotter class."""

from pathlib import Path

import pytest

from ref_builder.resources import RepoOTU
from ref_builder.snapshotter.snapshotter import Snapshotter

SNAPSHOT_AT_EVENTS = (
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
def snapshotter(snapshotter_otus: list[RepoOTU], tmp_path: Path) -> Snapshotter:
    """A snapshotter with eight OTUs already cached."""
    _snapshotter = Snapshotter(tmp_path / "snapshotter")

    for otu, at_event in zip(snapshotter_otus, SNAPSHOT_AT_EVENTS, strict=True):
        _snapshotter.cache_otu(otu, at_event)

    return _snapshotter


def test_load(snapshotter: Snapshotter, snapshotter_otus: list[RepoOTU]):
    """Test that the snapshotter loads snapshots containing the correct OTUs and
    at_event values.
    """
    for otu, at_event in zip(snapshotter_otus, SNAPSHOT_AT_EVENTS, strict=True):
        snapshot = snapshotter.load_by_id(otu.id)

        assert snapshot.at_event == at_event
        assert snapshot.otu == otu


def test_iter_otus(snapshotter: Snapshotter, snapshotter_otus: list[RepoOTU]):
    """Test that the snapshotter iterates over all OTUs."""
    assert sorted(snapshotter.iter_otus(), key=lambda otu: otu.id) == sorted(
        snapshotter_otus,
        key=lambda otu: otu.id,
    )


def test_otu_ids(snapshotter: Snapshotter, snapshotter_otus: list[RepoOTU]):
    """Test that the snapshotter has the correct OTU IDs."""
    assert set(snapshotter.otu_ids) == {otu.id for otu in snapshotter_otus}


class TestGetIDByLegacyID:
    """Test the `get_id_by_legacy_id` method of the Snapshotter class."""

    def test_ok(self, snapshotter: Snapshotter, snapshotter_otus: list[RepoOTU]):
        """Test that the correct OTU ID is retrieved by legacy ID."""
        for otu in snapshotter_otus:
            if otu.legacy_id is not None:
                assert snapshotter.get_id_by_legacy_id(otu.legacy_id) == otu.id

    def test_not_found(self, snapshotter: Snapshotter):
        """Test that `None` is returned when the legacy ID is not found."""
        assert snapshotter.get_id_by_legacy_id("not found") is None

    def test_none(self, snapshotter: Snapshotter):
        """Test that `None` is returned when the legacy ID is `None`."""
        assert snapshotter.get_id_by_legacy_id(None) is None


class TestGetIDByName:
    """Test the `get_id_by_name` method of the Snapshotter class."""

    def test_ok(self, snapshotter: Snapshotter, snapshotter_otus: list[RepoOTU]):
        """Test that the correct OTU ID is retrieved by name."""
        for otu in snapshotter_otus:
            assert snapshotter.get_id_by_name(otu.name) == otu.id

    def test_not_found(self, snapshotter: Snapshotter):
        """Test that `None` is returned when the name is not found."""
        assert snapshotter.get_id_by_name("not found") is None


class TestGetIDByTaxid:
    """Test the `get_id_by_taxid` method of the Snapshotter class."""

    def test_ok(self, snapshotter: Snapshotter, snapshotter_otus: list[RepoOTU]):
        """Test that the correct OTU ID is retrieved by taxid."""
        for otu in snapshotter_otus:
            assert snapshotter.get_id_by_taxid(otu.taxid) == otu.id

    def test_not_found(self, snapshotter: Snapshotter):
        """Test that `None` is returned when the taxid is not found."""
        assert snapshotter.get_id_by_taxid(999999999999999) is None

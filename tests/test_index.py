"""Tests for :class:`.EventIndex`."""

import uuid
from pathlib import Path

import pytest

from ref_builder.index import EventIndex, EventIndexItem


@pytest.mark.parametrize(
    "event_ids",
    ([1, 5, 20, 99, 110], [5, 20, 1, 99, 110]),
    ids=["sorted", "unsorted"],
)
def test_index(event_ids: list[int], tmp_path: Path):
    """Test that we can set, get, and load cached events for an OTU."""
    path = Path(tmp_path / "index")

    index = EventIndex(path)

    otu_id = uuid.uuid4()

    index.set(otu_id, event_ids, 112)

    assert (
        index.get(otu_id)
        == EventIndex(path).get(otu_id)
        == EventIndexItem(
            112,
            [1, 5, 20, 99, 110],
            otu_id,
        )
    )


def test_index_not_found(tmp_path: Path):
    """Test that we get ``None`` when an OTU ID is not found."""
    path = Path(tmp_path / "index")

    assert EventIndex(path).get(uuid.uuid4()) is None

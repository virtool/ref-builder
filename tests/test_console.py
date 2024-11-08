import pytest
from syrupy import SnapshotAssertion

from ref_builder.console import console, print_otu, print_otu_list
from ref_builder.models import OTUMinimal
from ref_builder.repo import Repo
from tests.fixtures.factories import OTUMinimalFactory


class TestPrintOTUList:
    """Tests for the ``print_otu_list`` function."""

    def test_ok(self, snapshot: SnapshotAssertion) -> None:
        """Test that listed OTUs are printed."""
        factory = OTUMinimalFactory()

        with console.capture() as capture:
            print_otu_list(factory.build() for _ in range(5))

        assert capture.get() == snapshot

    def test_empty(self) -> None:
        """Test that an empty list of OTUs is printed."""
        with console.capture() as capture:
            print_otu_list(OTUMinimal(**otu) for otu in [])

        assert capture.get() == "No OTUs found\n"


@pytest.mark.parametrize("taxid", [345184, 3158377])
def test_print_otu(scratch_repo: Repo, snapshot: SnapshotAssertion, taxid: int):
    """Test that an OTU is printed as expected by ``print_otu``."""
    otu = scratch_repo.get_otu_by_taxid(taxid)

    with console.capture() as capture:
        print_otu(otu)

    assert capture.get() == snapshot

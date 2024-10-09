import pytest

from ref_builder.console import print_otu, print_otu_list
from ref_builder.repo import Repo


def test_list_otus(scratch_repo: Repo):
    print_otu_list(scratch_repo.iter_minimal_otus())


@pytest.mark.parametrize("taxid", [345184, 3158377])
def test_print_otu(scratch_repo: Repo, taxid: int):
    otu = scratch_repo.get_otu_by_taxid(taxid)

    print_otu(otu)
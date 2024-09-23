import pytest
from factories import NCBISourceFactory

from ref_builder.ncbi.models import NCBISource


@pytest.mark.parametrize("dummy_source", NCBISourceFactory.coverage(10))
def test_source_factory(dummy_source: NCBISource):
    """Test that NCBISourceFactory creates valid input."""
    print(dummy_source)

    assert dummy_source.organism

from ref_builder.utils import generate_natural_sort_key


def test_generate_natural_sort_key():
    """Test that the function works as expected with ``sorted()``."""
    unsorted = [
        "RNA 1",
        "RNA 10",
        "RNA 2",
        "DNA 20",
        "RNA 4",
        "DNA 1",
        "RNA 8",
        "RNA 9",
        "RNA 1B",
        "DNA1",
    ]

    assert sorted(unsorted, key=generate_natural_sort_key) == [
        "DNA1",
        "DNA 1",
        "DNA 20",
        "RNA 1",
        "RNA 1B",
        "RNA 2",
        "RNA 4",
        "RNA 8",
        "RNA 9",
        "RNA 10",
    ]

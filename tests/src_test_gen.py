from pathlib import Path
from typing import List

import orjson
from pydantic import BaseModel, TypeAdapter

from ref_builder.repo import Repo, DataType
from ref_builder.otu import create_otu_with_schema, add_sequences


class OTUContents(BaseModel):
    taxid: int
    """The Taxonomy ID of the OTU."""

    otu_schema: list[str]
    """A list of accessions that determine the OTU's schema."""

    contents: list[str]
    """A list of accessions to be added to the OTU."""


otu_contents_list_adapter = TypeAdapter(List[OTUContents])
"""Aids the serialization of a scratch repo's table of contents."""


def generate_src_test_files(temp_path: Path, toc_path: Path) -> dict:
    """Initializes the scratch repo in a temporary path and serializes the events into a dict."""
    print("Generating new test files...")

    temp_scratch_repo = Repo.new(
        data_type=DataType.GENOME, name="src_test", path=temp_path, organism="viruses"
    )

    print(f"Temp repo initialized")

    with open(toc_path, "rb") as f:
        toc = otu_contents_list_adapter.validate_json(f.read())

    for otu_contents in toc:
        otu = create_otu_with_schema(
            repo=temp_scratch_repo,
            taxid=otu_contents.taxid,
            accessions=otu_contents.otu_schema
        )

        add_sequences(
            repo=temp_scratch_repo,
            otu=otu,
            accessions=otu_contents.contents
        )

    data = {}

    for event_file_path in (temp_path / "src").glob("*.json"):
        with open(event_file_path, "r") as f:
            data[event_file_path.name] = orjson.loads(f.read())

    print("Data copied")
    print(data.keys())

    return data

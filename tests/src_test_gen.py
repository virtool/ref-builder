import subprocess
from pathlib import Path
from typing import List

import orjson
from pydantic import BaseModel, TypeAdapter


class OTUContents(BaseModel):
    taxid: int
    otu_schema: list[str]
    contents: list[str]


ref_contents_adapter = TypeAdapter(List[OTUContents])


def generate_src_test_files(temp_path: Path, toc_path: Path) -> dict:
    print("Generating new test files...")

    subprocess.run(
        ["ref-builder", "init"]
        + ["--path", str(temp_path)]
        + ["--data-type", "genome",]
        + ["--name", "src_test",]
        + ["--organism", "viruses"],
        check=False,
    )

    print("Repo initialized")

    with open(toc_path, "rb") as f:
        toc = ref_contents_adapter.validate_json(f.read())

    for otu_contents in toc:
        subprocess.run(
            ["ref-builder", "otu"]
            + ["create", str(otu_contents.taxid)]
            + otu_contents.otu_schema +
            ["--path", str(temp_path)],
            check=False,
        )

        subprocess.run(
            ["ref-builder", "sequences", "create"]
            + ["--taxid", str(otu_contents.taxid)]
            + otu_contents.contents
            + ["--path", str(temp_path)],
            check=False,
        )

    data = {}

    for event_file_path in (temp_path / "src").glob("*.json"):
        with open(event_file_path, "r") as f:
            data[event_file_path.name] = orjson.loads(f.read())

    print("Data copied")
    print(data.keys())

    return data

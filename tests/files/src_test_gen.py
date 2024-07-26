import subprocess
import shutil
from pathlib import Path
from typing import List

from pydantic import BaseModel, TypeAdapter


class OTUContents(BaseModel):
    taxid: int
    otu_schema: list[str]
    contents: list[str]


ref_contents_adapter = TypeAdapter(List[OTUContents])


if __name__ == "__main__":
    root_path = Path(__file__).parent

    src_test_path = root_path / "src_test"

    src_test_temp_path = root_path / "src_test_temp"

    if src_test_temp_path.exists():
        shutil.rmtree(src_test_temp_path)
    src_test_temp_path.mkdir()

    subprocess.run(
        ["ref-builder", "init"]
        + ["--path", str(src_test_temp_path)]
        + ["--data-type", "genome",]
        + ["--name", "src_test",]
        + ["--organism", "viruses"],
        check=False,
    )

    contents_toc_path = root_path / "src_test_contents.json"
    with open(contents_toc_path, "rb") as f:
        toc = ref_contents_adapter.validate_json(f.read())

    for otu_contents in toc:
        subprocess.run(
            ["ref-builder", "otu"]
            + ["create", str(otu_contents.taxid)]
            + otu_contents.otu_schema +
            ["--path", str(src_test_temp_path)],
            check=False,
        )

        subprocess.run(
            ["ref-builder", "sequences", "create"]
            + ["--taxid", str(otu_contents.taxid)]
            + otu_contents.contents
            + ["--path", str(src_test_temp_path)],
            check=False,
        )

    shutil.rmtree(src_test_path)
    shutil.copytree(src_test_temp_path / "src", src_test_path)

    shutil.rmtree(src_test_temp_path)


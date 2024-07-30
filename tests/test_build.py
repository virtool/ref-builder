import json
from pathlib import Path

import arrow
import pytest
import orjson
from syrupy import SnapshotAssertion
from syrupy.filters import props

from ref_builder.build import build_json


@pytest.fixture
def comparison_reference(pytestconfig, tmp_path, scratch_repo) -> dict:
    comparison_reference = pytestconfig.cache.get("comparison_reference", None)

    if comparison_reference:
        return comparison_reference

    output_path = tmp_path / "comparison_reference.json"
    build_json(False, output_path, scratch_repo.path, "2.1.0")

    with open(output_path, "rb") as f:
        comparison_reference = orjson.loads(f.read())

    if comparison_reference:
        pytestconfig.cache.set("comparison_reference", comparison_reference)
        return comparison_reference

    raise ValueError("Could not build comparison reference")


def test_ok(
    scratch_repo, comparison_reference, tmp_path: Path, snapshot: SnapshotAssertion
):
    """Test that the command exits with a 0 exit code."""
    output_path = tmp_path / "reference.json"

    build_json(False, output_path, scratch_repo.path, "2.1.0")

    with open(output_path, "rb") as f:
        built_json = orjson.loads(f.read())

    assert comparison_reference

    assert built_json

    assert built_json == snapshot(exclude=props("created_at", "otus"))

    assert built_json["otus"] == comparison_reference["otus"]

    assert (arrow.utcnow() - arrow.get(built_json["created_at"])).seconds == 0


def test_indent(scratch_path: Path, tmp_path: Path):
    """Test that the indent in the reference.json file is properly set."""
    output_path = tmp_path / "reference.json"
    output_indented_path = tmp_path / "reference_indent.json"

    build_json(False, output_path, scratch_path, "2.1.0")
    build_json(True, output_indented_path, scratch_path, "2.1.0")

    assert {**json.load(output_path.open("rb")), "created_at": ""} == {
        **json.load(output_indented_path.open("rb")),
        "created_at": "",
    }

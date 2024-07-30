import json
from pathlib import Path

import arrow
from syrupy import SnapshotAssertion
from syrupy.filters import props

from ref_builder.build import build_json


def test_ok(
    scratch_repo,
    snapshot: SnapshotAssertion,
    tmp_path: Path,
):
    """Test that the command exits with a 0 exit code."""
    output_path = tmp_path / "reference.json"

    build_json(False, output_path, scratch_repo.path, "2.1.0")

    with open(output_path) as f:
        built_json = json.load(f)

    assert built_json == snapshot(exclude=props("created_at", "otus"))
    assert built_json["name"] == "2.1.0"

    for built_otu in built_json["otus"]:
        repo_otu = scratch_repo.get_otu(built_otu["_id"])
        assert repo_otu is not None

        assert built_otu == snapshot(exclude=props("_id", "isolates"))

        for built_isolate in built_otu["isolates"]:
            repo_isolate = repo_otu.get_isolate(built_isolate["id"])
            assert repo_isolate is not None

            assert built_isolate == snapshot(exclude=props("id", "sequences"))

            for built_sequence in built_isolate["sequences"]:
                assert built_sequence["_id"]
                assert built_sequence == snapshot(exclude=props("_id"))

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

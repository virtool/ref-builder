from pathlib import Path

import arrow
import orjson

from ref_builder.models import Molecule, Strandedness
from ref_builder.plan import SegmentRule
from ref_builder.repo import Repo


def _get_molecule_string(molecule: Molecule) -> str:
    """Return a string representation of the molecule."""
    string = ""

    if molecule.strandedness == Strandedness.SINGLE:
        string += "ss"
    else:
        string += "ds"

    if "DNA" in molecule.type:
        string += "DNA"
    else:
        string += "RNA"

    return string


def build_json(indent: bool, output_path: Path, path: Path, version: str) -> None:
    """Build a Virtool reference JSON file from a data directory.

    :param indent: whether to indent the JSON output
    :param output_path: The path to write the output JSON file to
    :param path: The path to a reference repository
    :param version: the version string to include in the reference.json file
    """
    repo = Repo(path)

    otus = []

    for otu in repo.iter_otus():
        isolates = []

        for isolate in otu.isolates:
            sequences = []

            for sequence in isolate.sequences:
                segment_name = otu.plan.get_segment_by_id(sequence.segment).name

                sequences.append(
                    {
                        "_id": sequence.legacy_id or str(sequence.id),
                        "accession": str(sequence.accession),
                        "definition": sequence.definition,
                        "host": "",
                        "segment": "Unnamed"
                        if segment_name is None
                        else str(segment_name),
                        "sequence": sequence.sequence,
                    }
                )

            isolates.append(
                {
                    "id": isolate.legacy_id or isolate.id,
                    "default": isolate.id == otu.representative_isolate,
                    "sequences": sequences,
                    "source_name": isolate.name.value
                    if isolate.name is not None
                    else "unknown",
                    "source_type": str(isolate.name.type)
                    if isolate.name is not None
                    else "unknown",
                },
            )

        isolates.sort(key=lambda x: x["id"])

        molecule = _get_molecule_string(otu.molecule)

        schema = [
            {
                "molecule": molecule,
                "name": str(segment.name) if segment.name else "Unnamed",
                "required": segment.required == SegmentRule.REQUIRED,
            }
            for segment in otu.plan.segments
        ]

        otus.append(
            {
                "_id": otu.legacy_id or otu.id,
                "abbreviation": otu.acronym,
                "isolates": isolates,
                "name": otu.name,
                "schema": schema,
                "taxid": otu.taxid,
            },
        )

    otus.sort(key=lambda x: x["name"])

    with open(output_path, "wb") as f:
        f.write(
            orjson.dumps(
                {
                    "created_at": arrow.utcnow().isoformat(),
                    "data_type": repo.meta.data_type,
                    "name": version,
                    "organism": repo.meta.organism,
                    "otus": otus,
                },
                option=orjson.OPT_INDENT_2 if indent else 0,
            ),
        )

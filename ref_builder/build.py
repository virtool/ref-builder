from pathlib import Path

import arrow
import orjson

from ref_builder.repo import Repo


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
            sequences = [
                {
                    "_id": sequence.legacy_id or sequence.id,
                    "accession": str(sequence.accession),
                    "definition": sequence.definition,
                    "host": "",
                    "segment": sequence.segment,
                    "sequence": sequence.sequence,
                }
                for sequence in isolate.sequences
            ]

            isolates.append(
                {
                    "id": isolate.legacy_id or isolate.id,
                    "default": isolate.id == otu.repr_isolate,
                    "sequences": sequences,
                    "source_name": isolate.name.value
                    if isolate.name is not None
                    else None,
                    "source_type": isolate.name.type
                    if isolate.name is not None
                    else None,
                },
            )

        isolates.sort(key=lambda x: x["id"])

        otus.append(
            {
                "_id": otu.legacy_id or otu.id,
                "abbreviation": otu.acronym,
                "isolates": isolates,
                "name": otu.name,
                "schema": otu.plan.model_dump(),
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

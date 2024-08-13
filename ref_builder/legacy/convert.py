import json
from pathlib import Path

from ref_builder.console import console
from ref_builder.legacy.utils import iter_legacy_otus
from ref_builder.logs import configure_logger
from ref_builder.models import Molecule
from ref_builder.ncbi.client import NCBIClient
from ref_builder.repo import Repo
from ref_builder.schema import OTUSchema
from ref_builder.utils import DataType, IsolateNameType


def convert_legacy_repo(name: str, path: Path, target_path: Path) -> None:
    """Convert the legacy repo at ``path`` to a new event-sourced repo at
    ``target_path``.

    A name is required for the new repo as legacy repos do not have names.

    :param name: the name of the new repo
    :param path: the path to the legacy repo
    :param target_path: the path to the new repo

    """
    configure_logger(False)

    with (path / "src" / "meta.json").open() as f:
        data = json.load(f)

    repo = Repo.new(
        DataType(data["data_type"]),
        name,
        target_path,
        data["organism"],
    )

    ncbi_client = NCBIClient(False)

    for otu in iter_legacy_otus(path / "src"):
        console.print(f"Converting {otu['name']} ({otu['taxid']})")

        first_unversioned_accession = otu["isolates"][0]["sequences"][0][
            "accession"
        ].split(".")[0]

        record = ncbi_client.fetch_genbank_records([first_unversioned_accession])[0]

        molecule = Molecule(
            strandedness=record.strandedness,
            topology=record.topology,
            type=record.moltype,
        )

        try:
            repo_otu = repo.create_otu(
                otu["abbreviation"],
                otu["_id"],
                otu["name"],
                OTUSchema(molecule=molecule, segments=otu["schema"]),
                otu["taxid"],
            )
        except ValueError as err:
            if "OTU already exists" in str(err):
                console.print(
                    f"[red]Error:[/red] Already exists: {otu['name']} ({otu['taxid']})",
                )
                continue

            raise

        for isolate in otu["isolates"]:
            source_type = isolate["source_type"]

            if source_type == "genbank" and isolate["sequences"][0][
                "accession"
            ].startswith("NC_"):
                source_type = "refseq"

            repo_isolate = repo.create_isolate_from_records(
                repo_otu.id,
                isolate["id"],
                isolate["source_name"],
                IsolateNameType(source_type),
            )

            for sequence in isolate["sequences"]:
                repo.create_sequence(
                    repo_otu.id,
                    repo_isolate.id,
                    sequence["accession"],
                    sequence["definition"],
                    sequence["_id"],
                    sequence["segment"],
                    sequence["sequence"],
                )

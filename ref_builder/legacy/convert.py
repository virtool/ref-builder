import json
import uuid
from pathlib import Path

from ref_builder.console import console
from ref_builder.legacy.utils import iter_legacy_otus
from ref_builder.logs import configure_logger
from ref_builder.models import Molecule
from ref_builder.ncbi.client import NCBIClient
from ref_builder.repo import Repo
from ref_builder.plan import MonopartitePlan, MultipartitePlan
from ref_builder.utils import DataType, IsolateName, IsolateNameType


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

        original_segments = [
            {**segment, "id": uuid.uuid4()} for segment in otu["schema"]
        ]
        if len(original_segments) > 2:
            isolate_plan = MultipartitePlan.model_validate(
                {
                    "id": uuid.uuid4(),
                    "segments": original_segments,
                }
            )
        else:
            isolate_plan = MonopartitePlan(
                id=original_segments[0]["id"], length=len(record.sequence)
            )

        try:
            try:
                repo_otu = repo.create_otu(
                    otu["abbreviation"],
                    otu["_id"],
                    molecule=molecule,
                    name=otu["name"],
                    plan=isolate_plan,
                    taxid=otu["taxid"],
                )
            except ValueError as err:
                if "OTU already exists" in str(err):
                    console.print(
                        f"[red]Error:[/red] Already exists: {otu['name']} ({otu['taxid']})",
                    )
                    continue

                raise
        except ValueError as err:
            if "OTU already exists" in str(err):
                console.print(
                    f"[red]Error:[/red] Already exists: {otu['name']} ({otu['taxid']})",
                )
                continue

            raise

        for isolate in otu["isolates"]:
            if isolate["source_type"] in ("genbank", "refseq"):
                name = None
            else:
                name = IsolateName(
                    IsolateNameType(isolate["source_type"]),
                    isolate["source_name"],
                )

            repo_isolate = repo.create_isolate(
                repo_otu.id,
                isolate["id"],
                name,
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

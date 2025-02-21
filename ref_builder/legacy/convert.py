import json
from collections import defaultdict
from pathlib import Path

from ref_builder.console import console
from ref_builder.legacy.utils import iter_legacy_otus
from ref_builder.logs import configure_logger
from ref_builder.models import Molecule
from ref_builder.ncbi.client import NCBIClient
from ref_builder.plan import Plan, Segment, SegmentName, SegmentRule
from ref_builder.repo import Repo
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

        sequences_by_segment = defaultdict(list)

        for isolate in otu["isolates"]:
            for legacy_sequence in isolate["sequences"]:
                sequences_by_segment[legacy_sequence["segment"]].append(legacy_sequence)

        segments = []

        for segment in otu["schema"]:
            sequences = sequences_by_segment[segment["name"]]

            if segment["name"] in ("Genome", "Partial"):
                name = None
            else:
                try:
                    name = SegmentName.from_string(segment["name"])
                except ValueError:
                    name = SegmentName(key=segment["name"], prefix=molecule.type)

            lengths = [len(sequence["sequence"]) for sequence in sequences]

            mean_length = sum(lengths) / len(sequences)

            # Calculate tolerance based on min, max, and mean observed lengths.
            length_tolerance = max(
                0.01,
                (
                    max(mean_length - min(lengths), max(lengths) - mean_length)
                    / mean_length
                ),
            )

            new_segment = Segment.new(
                length=int(mean_length),
                length_tolerance=round(length_tolerance, 2),
                name=name,
                rule=SegmentRule.REQUIRED
                if segment["required"]
                else SegmentRule.OPTIONAL,
            )

            for legacy_sequence in sequences:
                legacy_sequence["segment"] = new_segment.id

            segments.append(new_segment)

        plan = Plan.new(segments)

        repo_otu = repo.create_otu(
            otu["abbreviation"],
            otu["_id"],
            molecule=molecule,
            name=otu["name"],
            plan=plan,
            taxid=otu["taxid"],
        )

        default_isolate_ids = []

        for isolate in otu["isolates"]:
            if isolate["source_type"] in ("genbank", "partial", "refseq", "unknown"):
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

            if isolate["default"]:
                default_isolate_ids.append(repo_isolate.id)

            for legacy_sequence in isolate["sequences"]:
                repo_sequence = repo.create_sequence(
                    repo_otu.id,
                    legacy_sequence["accession"],
                    legacy_sequence["definition"],
                    legacy_sequence["_id"],
                    legacy_sequence["segment"],
                    legacy_sequence["sequence"],
                )

                repo.link_sequence(
                    repo_otu.id,
                    repo_isolate.id,
                    repo_sequence.id,
                )

        if len(default_isolate_ids) > 1:
            raise ValueError("More than one default isolate found.")

        if default_isolate_ids:
            repo.set_representative_isolate(repo_otu.id, default_isolate_ids[0])
        else:
            raise ValueError("No default isolate found.")

        repo_otu = repo.get_otu(repo_otu.id)

        if not repo_otu:
            raise ValueError("Could not retrieve OTU from repo.")

        representative_isolate = repo_otu.representative_isolate

        if not repo_otu.representative_isolate:
            raise ValueError("Could not retrieve representative isolate from OTU.")

        if not repo_otu.get_isolate(representative_isolate):
            raise ValueError("Could not retrieve representative isolate from OTU.")

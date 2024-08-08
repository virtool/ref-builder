"""Repo level checks for legacy repositories."""

from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

from ref_builder.legacy.models import LegacySourceType
from ref_builder.legacy.utils import console, iter_legacy_otus
from ref_builder.ncbi.client import NCBIClient
from ref_builder.otu import get_isolate_name

type LegacyIsolateName = tuple[LegacySourceType, str]
"""A tuple containing the legacy source type and name."""


@dataclass
class IsolateNameCheckResult:
    """The result of checking an existing legacy isolate name against NCBI."""

    current_name: LegacyIsolateName
    """The current name of the isolate."""

    found_names: list[LegacyIsolateName]
    """The isolate names found for the NCBI accessions in the isolate."""

    @property
    def match(self) -> bool:
        """Whether the name found in NCBI matches the name in the isolate."""
        return self.current_name in self.found_names


def check_isolate_names_against_ncbi(path: Path) -> dict[str, IsolateNameCheckResult]:
    """Find isolates whose names don't match what is found in NCBI.

    The following possible results are reported:

    1. One name found. It doesn't match the current one.
    2. Multiple names found. One matches the current one.
    3. Multiple names found. None match the current one.
    4. No names found.

    If a single name is found, and it matches the current one, the isolate is not
    reported.

    :param path: the path to the legacy reference.
    :returns: a dictionary mapping isolate IDs to the results of the check.
    """
    ncbi_client = NCBIClient(path / ".migration_cache", False)

    result = {}

    for otu in iter_legacy_otus(path / "src"):
        result = {}

        for isolate in otu["isolates"]:
            accessions = [sequence["accession"] for sequence in isolate["sequences"]]

            records = ncbi_client.fetch_genbank_records(accessions)

            found_names: list[LegacyIsolateName] = [
                (LegacySourceType(n.type), n.value)
                for n in {get_isolate_name(record) for record in records}
                if n is not None
            ]

            current_name: LegacyIsolateName = (
                LegacySourceType(isolate["source_type"]),
                isolate["source_name"],
            )

            if len(found_names) == 1 and found_names[0] == current_name:
                continue

            result[isolate["id"]] = IsolateNameCheckResult(
                current_name=current_name,
                found_names=found_names,
            )

    return result


def find_duplicate_accessions(path: Path) -> dict[str, list[str]]:
    """Find accessions that are duplicate in a legacy reference.

    Return a dictionary mapping accessions to lists of the IDs of sequences they appear
    in. Unique accessions are not included in the output.

    :param path: the path to the legacy reference.
    :return: a dictionary mapping accessions to lists of the sequence IDs

    """
    accessions_to_sequence_ids = defaultdict(list)

    for otu in iter_legacy_otus(path / "src"):
        for isolate in otu["isolates"]:
            for sequence in isolate["sequences"]:
                accessions_to_sequence_ids[sequence["accession"]].append(
                    sequence["_id"],
                )

    return {
        accession: sorted(sequence_ids)
        for accession, sequence_ids in accessions_to_sequence_ids.items()
        if len(sequence_ids) > 1
    }


def check_unique_ids(path: Path) -> None:
    """Validate that resource IDs are unique.

    * All OTU _id fields are unique.
    * All isolate id fields are unique.
    * All sequence _id fields are unique.
    """
    otu_ids = []
    isolate_ids = []
    sequence_ids = []

    for otu in iter_legacy_otus(path / "src"):
        otu_ids.append(otu["_id"])

        for isolate in otu["isolates"]:
            isolate_ids.append(isolate["id"])
            sequence_ids += [sequence["_id"] for sequence in isolate["sequences"]]

    if len(otu_ids) != len(set(otu_ids)):
        console.log("Found non-unique OTU ids")

    if len(isolate_ids) != len(set(isolate_ids)):
        console.log("Found non-unique isolate ids")

    if len(sequence_ids) != len(set(sequence_ids)):
        console.log("Found non-unique sequence ids")


def check_unique_otu_abbreviations_and_names(path: Path) -> None:
    """Validate that OTU abbreviations and names are unique."""
    abbreviations = set()
    names = set()

    for otu in iter_legacy_otus(path / "src"):
        abbreviation = otu["abbreviation"]

        if abbreviation:
            if abbreviation in abbreviations:
                console.log(f"Found non-unique OTU abbreviation: {abbreviation}")
            else:
                abbreviations.add(abbreviation)

        name = otu["name"].lower()

        if name in names:
            console.log(f"Found non-unique OTU name: {name}")
        else:
            names.add(name)

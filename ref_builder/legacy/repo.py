"""Repo level checks for legacy repositories."""

from pathlib import Path

from ref_builder.console import console
from ref_builder.legacy.utils import iter_legacy_otus


def check_unique_accessions(path: Path) -> None:
    """Validate that all accession fields are unique."""
    accessions = set()

    for otu in iter_legacy_otus(path / "src"):
        for isolate in otu["isolates"]:
            for sequence in isolate["sequences"]:
                if accession := sequence.get("accession"):
                    if accession in accessions:
                        console.log(f"Found non-unique accession: {accession}")
                    else:
                        accessions.add(accession)


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

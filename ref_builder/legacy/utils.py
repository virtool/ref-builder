import re
from collections.abc import Generator
from dataclasses import dataclass
from pathlib import Path
from random import choice
from string import ascii_letters, ascii_lowercase, digits

import orjson
from pydantic_core import ErrorDetails

from ref_builder.legacy.models import LegacyIsolateSource, LegacySourceType
from ref_builder.ncbi.client import NCBIClient
from ref_builder.ncbi.models import NCBIGenbank


class HandleErrorContext:
    """A context that is passed to a ``ValidationError`` error handler.

    This context is for one error from ``ValidationError.errors``, not for the
    ``ValidationError`` itself.

    """

    def __init__(
        self,
        error: ErrorDetails,
        fix: bool,
        ncbi_client: NCBIClient,
        otu: dict,
        repaired_otu: dict,
    ) -> None:
        self.error = error
        """The details of the Pydantic validation error."""

        self.fix = fix
        """Whether to attempt to fix the error or not."""

        self.ncbi_client = ncbi_client
        """An instance of the NCBI client."""

        self.otu = otu
        """A dictionary representation of the OTU."""

        self.repaired_otu = repaired_otu
        """A copy of the OTU that repairs should be applied to."""

    def update_otu(self, update: dict):
        """Update OTU associated with the error.

        Example:
        -------
        .. code-block:: python

            ctx.update_otu({"name": "new_name"})

        :param update: the update dict

        """
        self.repaired_otu.update(update)

    def update_isolate(self, isolate_index: int, update: dict) -> None:
        """Update a legacy isolate on disk given an update dict.

        Example:
        -------
        .. code-block:: python

           ctx.update_isolate(
               isolate_index,
              {"source_name": source.name, "source_type": source.type.lower()},
           )

        :param isolate_index: the index of the isolate in the OTU
        :param update: the update dict

        """
        self.repaired_otu["isolates"][isolate_index].update(update)

    def update_sequence(
        self,
        isolate_index: int,
        sequence_index: int,
        update: dict,
    ) -> None:
        """Update a legacy sequence on disk given an update dict.

        Example:
        -------
        .. code-block:: python

            ctx.update_sequence(
                isolate_index,
                sequence_index,
                {"sequence": "ATGCATCGTATGCTYACTAGCTAC"},
            )

        """
        self.repaired_otu["isolates"][isolate_index]["sequences"][
            sequence_index
        ].update(
            update,
        )


@dataclass
class ErrorHandledResult:
    """The return value for a legacy validation ``ValidationError`` handler."""

    message: str
    """The message to display to the user.

    Includes ``rich``-compatible markup.
    """

    fixed: bool = False
    """Whether the error was fixed.

    This value is used to show the user whether the error was fixed or not.
    """


def build_legacy_otu(path: Path) -> dict:
    """Build a dict-based representation of an OTU from an OTU directory.

    :param path: the path to the OTU directory
    :return: a dict-based representation of the OTU

    """
    with (path / "otu.json").open("r") as f:
        otu_data = orjson.loads(f.read())

    otu_data["isolates"] = []

    for isolate_dir_path in sorted(path.iterdir()):
        if isolate_dir_path.name == "otu.json":
            continue

        with open(isolate_dir_path / "isolate.json") as f:
            isolate_data = orjson.loads(f.read())

        isolate_data["sequences"] = []

        for sequence_path in sorted(isolate_dir_path.iterdir()):
            if sequence_path.name == "isolate.json":
                continue

            with open(sequence_path) as f:
                sequence_data = orjson.loads(f.read())

            isolate_data["sequences"].append(sequence_data)

        otu_data["isolates"].append(isolate_data)

    return otu_data


def replace_otu(path: Path, otu: dict) -> None:
    """Replace the otu in the path with the data in ``otu``.

    :param path: the path to the otu
    :param otu: the new otu data

    """
    if not path.parent.is_dir():
        raise FileNotFoundError(f"No such OTU parent directory: {path.parent}")

    if not path.is_dir():
        raise FileNotFoundError(f"No such OTU directory: {path}")

    with open(path / "otu.json", "wb") as f:
        f.write(
            orjson.dumps(
                {
                    "_id": otu["_id"],
                    "name": otu["name"],
                    "abbreviation": otu["abbreviation"],
                    "schema": otu["schema"],
                    "taxid": otu["taxid"],
                },
                option=orjson.OPT_INDENT_2,
            ),
        )

    for isolate in sorted(otu["isolates"], key=lambda x: x["id"]):
        isolate_id = isolate["id"]

        with open(path / isolate_id / "isolate.json", "wb") as f:
            f.write(
                orjson.dumps(
                    {
                        "id": isolate["id"],
                        "default": isolate["default"],
                        "source_name": isolate["source_name"],
                        "source_type": isolate["source_type"],
                    },
                    option=orjson.OPT_INDENT_2,
                ),
            )

        for sequence in sorted(isolate["sequences"], key=lambda x: x["_id"]):
            with open(path / isolate_id / f"{sequence['_id']}.json", "wb") as f:
                f.write(
                    orjson.dumps(
                        {
                            "_id": sequence["_id"],
                            "accession": sequence["accession"],
                            "definition": sequence["definition"],
                            "host": sequence["host"],
                            "segment": sequence["segment"],
                            "sequence": sequence["sequence"],
                        },
                        option=orjson.OPT_INDENT_2,
                    ),
                )


def iter_legacy_otus(src_path: Path) -> Generator[dict, None, None]:
    """Iterate over all OTUs in a legacy repository.

    :param src_path: The path to the legacy repository.
    :return: a generator that yields OTU data dictionaries.
    """
    for alpha_dir_path in sorted(src_path.iterdir()):
        if alpha_dir_path.name == "meta.json":
            continue

        for otu_dir_path in sorted(alpha_dir_path.iterdir()):
            yield build_legacy_otu(otu_dir_path)


def generate_random_alphanumeric(
    length: int = 8,
    mixed_case: bool = False,
    excluded: list | set | None = None,
) -> str:
    """Generates a random string composed of letters and numbers.

    :param length: the length of the string.
    :param mixed_case: included alpha characters will be mixed case instead of lowercase
    :param excluded: strings that may not be returned.
    :return: a random alphanumeric string.
    """
    excluded_set = set() if excluded is None else set(excluded)

    characters = digits + (ascii_letters if mixed_case else ascii_lowercase)

    candidate = "".join([choice(characters) for _ in range(length)])

    if candidate not in excluded_set:
        return candidate

    return generate_random_alphanumeric(length=length, excluded=excluded_set)


def generate_otu_dirname(name: str, otu_id: str = "") -> str:
    """Takes in a human-readable string, replaces whitespace and symbols
    and adds the Virtool hash id as a suffix

    :param name: Human-readable, searchable name of the OTU
    :param otu_id: ID hash of OTU
    :return: A directory name in the form of 'converted_otu_name--taxid'
    """
    no_plus = name.replace("+", "plus ")
    no_symbols = re.split(r"[():/-]+", no_plus)
    joined = " ".join(no_symbols)
    no_whitespace = re.sub(r"[\s]+", "_", joined)

    dirname = no_whitespace.lower()
    dirname += "--" + otu_id

    return dirname


def generate_unique_ids(
    n: int = 1,
    length: int = 8,
    mixed_case: bool = False,
    excluded: list | set | None = None,
) -> set:
    """Generates an n-length list of unique alphanumeric IDs.

    :param n: The number of strings to be generated
    :param mixed_case: included alpha characters will be mixed case instead of lowercase
    :param length: The length of each string
    :param excluded: List of alphanumeric strings that should be excluded from generation
    """
    if excluded is None:
        excluded_set = set()
    else:
        excluded_set = set(excluded)

    new_uniques = set()
    while len(new_uniques) < n:
        new_uniques.add(generate_random_alphanumeric(length, mixed_case, excluded_set))

    return new_uniques


def is_v1(src_path: Path) -> bool:
    """Returns True if the given reference directory is formatted under version 1 guidelines.
    Determines if virtool ref migrate should be run before using this version of virtool-cli

    :param src_path: Path to a src database directory
    :return: Boolean value depending on whether alphabetized bins are found.
    """
    alpha_bins = list(src_path.glob("[a-z]"))

    return bool(alpha_bins)


def extract_isolate_source(
    genbank_records: list[NCBIGenbank],
) -> LegacyIsolateSource:
    """Extract a legacy isolate source from a set of Genbank records associated with the
    isolate.
    """
    for record in genbank_records:
        if record.source.isolate:
            return LegacyIsolateSource(
                name=record.source.isolate,
                type=LegacySourceType.ISOLATE,
            )

        if record.source.strain:
            return LegacyIsolateSource(
                name=record.source.strain,
                type=LegacySourceType.STRAIN,
            )

        if record.source.clone:
            return LegacyIsolateSource(
                name=record.source.clone,
                type=LegacySourceType.CLONE,
            )

    accessions = sorted(
        (record.accession for record in genbank_records if record.accession),
        key=lambda x: int(x.replace("NC_", "").replace(".1", "")),
    )

    return LegacyIsolateSource(
        name=accessions[0].upper(),
        type=LegacySourceType.GENBANK,
    )

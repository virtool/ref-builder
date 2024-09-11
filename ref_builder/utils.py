from dataclasses import dataclass
from enum import StrEnum
from pathlib import Path
from typing import NamedTuple

import orjson

ZERO_PADDING_MAX = 99999999
"""The maximum number that can be padded with zeroes in event IDs and filenames."""


class Accession(NamedTuple):
    """A Genbank accession number."""

    key: str
    """The accession key.

    In the accession "MN908947.3", the key is "MN908947".
    """

    version: int
    """The version number.

    In the accession "MN908947.3", the version is 3.
    """

    @classmethod
    def from_string(cls, string: str) -> "Accession":
        """Create an Accession from a raw accession string."""
        accession_parts = string.split(".")

        if len(accession_parts) == 2:
            try:
                return Accession(
                    key=accession_parts[0],
                    version=int(accession_parts[1]),
                )
            except ValueError:
                raise ValueError(
                    f"Raw accession {string} does not include a valid version.",
                )

        msg = f"Raw accession {string} is not versioned."
        raise ValueError(msg)

    def __str__(self) -> str:
        """Return the accession as a string."""
        return f"{self.key}.{self.version}"


class DataType(StrEnum):
    """Possible data types for a reference repository."""

    BARCODE = "barcode"
    GENOME = "genome"


class IsolateNameType(StrEnum):
    """Possible types for isolate names.

    **Ordered by priority**. Do not reorder attributes.

    Isolate name types were previously called "source types". They are referred to this
    way in Virtool.
    """

    ISOLATE = "isolate"
    STRAIN = "strain"
    CLONE = "clone"
    VARIANT = "variant"
    GENOTYPE = "genotype"
    SEROTYPE = "serotype"


@dataclass(frozen=True)
class IsolateName:
    """A name for an isolate.

    The isolate name consists of a type and a value.

    For example, in the isolate name "Isolate PPSMV2-Badnapur", the type is "isolate"
    and the value is "PPSMV2-Badnapur.

    In the Genbank record for a sequence that belongs to this isolate, the source table
    contains:

    .. code-block:: text

        /isolate="PPSMV2-Badnapur"

    """

    type: IsolateNameType
    """The type of sub-species categorization."""

    value: str
    """The name of this subcategory."""

    def __str__(self) -> str:
        """Return the isolate name as a formatted string."""
        return f"{self.type.capitalize()} {self.value}"


def format_json(path: Path) -> None:
    """Format the JSON file at ``path`` in place.

    :param path: the path to the JSON file to format
    """
    with path.open("rb") as f:
        data = orjson.loads(f.read())

    with path.open("wb") as f:
        f.write(orjson.dumps(data, option=orjson.OPT_INDENT_2))


def pad_zeroes(number: int) -> str:
    """Pad a number with zeroes to make it 8 characters long.

    :param number: the number to pad
    :return: the padded number
    """
    if number > ZERO_PADDING_MAX:
        raise ValueError("Number is too large to pad")

    return str(number).zfill(8)

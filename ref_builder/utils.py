import re
from dataclasses import dataclass
from enum import StrEnum
from pathlib import Path

import orjson

ZERO_PADDING_MAX = 99999999
"""The maximum number that can be padded with zeroes in event IDs and filenames."""

GENBANK_ACCESSION_PATTERN = re.compile(pattern=r"^[A-Z]{1,2}[0-9]{5,6}$")
REFSEQ_ACCESSION_PATTERN = re.compile(pattern=r"^NC_[0-9]{6}$")


@dataclass(frozen=True)
class Accession:
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
        try:
            key, string_version = string.split(".")
        except ValueError as e:
            if "not enough values to unpack" in str(e):
                raise ValueError(
                    "Accession does not contain two parts delimited by a period."
                )

            raise

        try:
            version = int(string_version)
        except ValueError as e:
            if "invalid literal for int() with base 10" in str(e):
                raise ValueError("Accession version is not an integer.")

            raise

        return Accession(key=key, version=version)

    def __eq__(self, other: "Accession") -> bool:
        if isinstance(other, Accession):
            return self.key == other.key and self.version == other.version

    def __lt__(self, other: "Accession") -> bool:
        if isinstance(other, Accession):
            return self.key < other.key or self.version < other.version

    def __gt__(self, other: "Accession") -> bool:
        if isinstance(other, Accession):
            return self.key > other.key or self.version > other.version

    def __str__(self) -> str:
        """Return the accession as a string."""
        return f"{self.key}.{self.version}"


class ExcludedAccessionAction(StrEnum):
    """Possible actions that can be taken on the excluded/allowed status of an accession."""

    ALLOW = "allow"
    EXCLUDE = "exclude"


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


def is_accession_key_valid(accession_key: str) -> bool:
    """Return True if the given accession is a valid Genbank or RefSeq accession."""
    return (
        GENBANK_ACCESSION_PATTERN.match(accession_key) is not None
        or REFSEQ_ACCESSION_PATTERN.match(accession_key) is not None
    )


def get_accession_key(raw: str) -> str:
    """Parse a string to check if it follows Genbank or RefSeq accession parameters and
    return the key part only.
    """
    if is_accession_key_valid(raw):
        return raw

    try:
        versioned_accession = Accession.from_string(raw)
    except ValueError:
        raise ValueError("Invalid accession key")

    if GENBANK_ACCESSION_PATTERN.match(
        versioned_accession.key
    ) or REFSEQ_ACCESSION_PATTERN.match(versioned_accession.key):
        return versioned_accession.key

    raise ValueError("Invalid accession key")


def generate_natural_sort_key(string: str) -> list[int | str]:
    """Generate a natural order sorting key for a string.

    This list: ["1", "10", "2"] will be sorted as ["1", "2", "10"], as opposed to
    lexographical order, which would result in ["1", "10", "2"].

    :param string: the string to convert to a sorting key
    :return: the sorting key
    """

    def _convert(string_: str) -> int | str:
        return int(string_) if string_.isdigit() else string_.lower()

    return [_convert(c) for c in re.split("([0-9]+)", string)]


def is_refseq(accession_key: str) -> bool:
    """Return True if accession is RefSeq."""
    return re.match(REFSEQ_ACCESSION_PATTERN, accession_key) is not None


def pad_zeroes(number: int) -> str:
    """Pad a number with zeroes to make it 8 characters long.

    :param number: the number to pad
    :return: the padded number
    """
    if number > ZERO_PADDING_MAX:
        raise ValueError("Number is too large to pad")

    return str(number).zfill(8)

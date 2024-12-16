"""Command-line interface for reference builder."""

import glob
from pathlib import Path

import click
import structlog

from ref_builder.cli.otu import otu
from ref_builder.build import build_json
from ref_builder.legacy.convert import convert_legacy_repo
from ref_builder.legacy.utils import iter_legacy_otus
from ref_builder.legacy.validate import validate_legacy_repo
from ref_builder.logs import configure_logger
from ref_builder.ncbi.client import NCBIClient
from ref_builder.options import (
    legacy_repo_path_option,
    path_option,
)
from ref_builder.repo import Repo
from ref_builder.utils import DataType, format_json

PRECACHE_BUFFER_SIZE = 450
"""The number of Genbank records to fetch in a single request when pre-caching."""

logger = structlog.get_logger()


@click.group()
@click.option("--debug", is_flag=True, help="Show debug logs")
@click.option("-v", "--verbose", "verbosity", count=True)
def entry(debug: bool, verbosity: int) -> None:
    """Build and maintain reference sets of pathogen genome sequences."""
    if debug:
        verbosity = 2

    configure_logger(verbosity)


@entry.command()
@click.option(
    "--data-type",
    help="The type of data the reference will contain (eg. genome)",
    required=True,
    type=click.Choice(DataType),
)
@click.option(
    "--name",
    help="A name for the reference",
    required=True,
    type=str,
)
@click.option(
    "--organism",
    default="",
    help="The organism the reference is for (eg. virus)",
    type=str,
)
@click.option(
    "--path",
    default=".",
    help="The path to initialize the repository at",
    type=click.Path(path_type=Path),
)
def init(data_type: DataType, name: str, organism: str, path: Path) -> None:
    """Create a new reference repository.

    If a path is not provided, the repository will be created in the current directory.

    Repositories cannot be created in non-empty directories.
    """
    Repo.new(data_type, name, path, organism)


entry.add_command(otu)


@entry.group()
def legacy() -> None:
    """Validate and convert legacy references."""


@legacy.command()
@click.option(
    "--name",
    help="the name to set for the new repository",
    required=True,
    type=str,
)
@click.option(
    "--target-path",
    help="the path for the new converted repository",
    required=True,
    type=click.Path(exists=False, file_okay=False, path_type=Path),
)
@legacy_repo_path_option
def convert(name: str, path: Path, target_path: Path) -> None:
    """Convert a legacy reference repository to a new event-sourced repository."""
    convert_legacy_repo(name, path, target_path)


@legacy.command()
@legacy_repo_path_option
def precache(path: Path) -> None:
    """Pre-cache all accessions in a legacy reference repository."""
    ncbi = NCBIClient(False)

    buffer = []

    for otu_ in iter_legacy_otus(path / "src"):
        for isolate in otu_["isolates"]:
            for sequence in isolate["sequences"]:
                buffer.append(sequence["accession"])

                if len(buffer) > PRECACHE_BUFFER_SIZE:
                    ncbi.fetch_genbank_records(buffer)
                    buffer = []

    if buffer:
        ncbi.fetch_genbank_records(buffer)


@legacy.command(name="format")
@legacy_repo_path_option
def reformat(path: Path) -> None:
    """Format a legacy reference repository.

    Re-formats every JSON file in a legacy reference repository to have a consistent
    format.
    """
    src_path = path / "src"

    for json_path in glob.glob(str(src_path / "**/*.json"), recursive=True):
        format_json(Path(json_path))


@legacy.command()
@click.option(
    "--fix",
    is_flag=True,
    help="Attempt to fix errors in-place",
)
@click.option(
    "--limit",
    default=0,
    help="Exit if this many OTUs have errors",
)
@click.option(
    "--no-ok",
    is_flag=True,
    help="Don't print anything if an otu is valid",
)
@legacy_repo_path_option
def validate(fix: bool, limit: int, no_ok: bool, path: Path) -> None:
    """Validate a legacy reference repository."""
    validate_legacy_repo(fix, limit, no_ok, path)


@entry.command()
@click.option(
    "-i",
    "--indent",
    is_flag=True,
    help="Indent the output JSON file",
)
@click.option(
    "-V",
    "--version",
    default="",
    type=str,
    help="A version string to include in the output file",
)
@click.option(
    "-o",
    "--target-path",
    required=True,
    type=click.Path(exists=False, file_okay=True, path_type=Path),
    help="The path to write the reference.json file to",
)
@path_option
def build(indent: bool, path: Path, target_path: Path, version: str) -> None:
    """Build a Virtool reference.json file from the reference repository."""
    build_json(indent, target_path, path, version)

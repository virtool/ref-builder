"""Command-line interface for reference builder."""

import glob
import sys
from pathlib import Path
from uuid import UUID

import click

from ref_builder.build import build_json
from ref_builder.console import print_otu, print_otu_list
from ref_builder.legacy.utils import iter_legacy_otus
from ref_builder.legacy.validate import validate_legacy_repo
from ref_builder.logging import configure_logger
from ref_builder.ncbi.client import NCBIClient
from ref_builder.options import debug_option, ignore_cache_option, path_option
from ref_builder.otu import (
    add_sequences,
    create_otu,
    create_otu_with_schema,
    update_otu,
)
from ref_builder.repo import Repo
from ref_builder.utils import DataType, format_json


@click.group()
def entry() -> None:
    """Build and maintain reference sets of pathogen genome sequences."""


@entry.command()
@click.option(
    "--data-type",
    help="the type of data the reference contains (eg. genome)",
    required=True,
    type=click.Choice(DataType),
)
@click.option(
    "--name",
    help="the type of data the reference contains (eg. genome)",
    required=True,
    type=str,
)
@click.option(
    "--organism",
    default="",
    help="the organism the reference is for (eg. virus)",
    type=str,
)
@click.option(
    "--path",
    default=".",
    help="the path to initialize the repository at",
    type=click.Path(path_type=Path),
)
def init(data_type: DataType, name: str, organism: str, path: Path) -> None:
    """Create a new event-sourced repo."""
    Repo.new(data_type, name, path, organism)


@entry.group()
def otu() -> None:
    """Manage OTUs."""


@otu.command(name="create")
@click.argument("TAXID", type=int)
@click.argument(
    "accessions_",
    metavar="ACCESSIONS",
    nargs=-1,
    type=str,
)
@click.option("--autofill/--no-fill", default=False)
@ignore_cache_option
@debug_option
@path_option
def otu_create(
    debug: bool,
    ignore_cache: bool,
    path: Path,
    taxid: int,
    accessions_: list[str],
    autofill: bool,
) -> None:
    """Create a new OTU for the given Taxonomy ID and accessions."""
    configure_logger(debug)

    repo = Repo(path)

    if accessions_:
        try:
            new_otu = create_otu_with_schema(
                repo,
                taxid,
                accessions_,
                ignore_cache=False,
            )
        except ValueError as e:
            click.echo(e, err=True)
            sys.exit(1)

    else:
        click.echo("No accessions given. Creating without schema...", err=True)
        try:
            new_otu = create_otu(repo, taxid, ignore_cache=False)
        except ValueError as e:
            click.echo(e, err=True)
            sys.exit(1)

    if autofill:
        update_otu(repo, new_otu, ignore_cache=ignore_cache)


@otu.command(name="update")
@click.argument("TAXID", type=int)
@debug_option
@ignore_cache_option
@path_option
def otu_update(debug: bool, ignore_cache: bool, path: Path, taxid: int) -> None:
    """Update an existing OTU with the latest data from NCBI."""
    configure_logger(debug)

    repo = Repo(path)

    otu_ = repo.get_otu_by_taxid(taxid)

    if otu_ is None:
        click.echo(f"OTU not found for Taxonomy ID {taxid}.", err=True)
        click.echo(f'Run "virtool otu create {taxid} --autofill" instead.')
        sys.exit(1)

    update_otu(repo, otu_, ignore_cache=ignore_cache)


@otu.command(name="get")
@click.argument("IDENTIFIER", type=str)
@path_option
def otu_get(identifier: str, path: Path) -> int:
    """Get an OTU by its Taxonomy ID."""
    try:
        identifier = int(identifier)
        otu_ = Repo(path).get_otu_by_taxid(identifier)
    except ValueError:
        otu_ = Repo(path).get_otu(UUID(identifier))

    if otu_ is None:
        click.echo("No OTU found.", err=True)
        return 1

    print_otu(otu_)

    return 0


@otu.command(name="list")
@path_option
def otu_list(path: Path) -> None:
    """List all OTUs in the repository."""
    print_otu_list(Repo(path).iter_otus())


@entry.group()
def sequences() -> None:
    """Manage sequences."""


@sequences.command(name="create")
@click.argument(
    "accessions_",
    metavar="ACCESSIONS",
    nargs=-1,
    type=str,
)
@click.option("--taxid", type=int, required=True)
@ignore_cache_option
@path_option
@debug_option
def create_sequence(
    accessions_: list[str],
    debug: bool,
    ignore_cache: bool,
    path: Path,
    taxid: int,
) -> None:
    """Fetch and write the data for the given NCBI accessions to an OTU.

    virtool add accessions --taxid 2697049 MN996528.1 --path [repo_path]

    """
    configure_logger(debug)

    repo = Repo(path)

    existing_otu = repo.get_otu_by_taxid(taxid)
    if existing_otu is None:
        click.echo(f"OTU not found for Taxonomy ID {taxid}.", err=True)
        click.echo(f'Run "virtool otu create {taxid} --path {path} --autofill" instead')
        sys.exit(1)

    add_sequences(repo, existing_otu, accessions=accessions_, ignore_cache=ignore_cache)


@entry.group()
def legacy() -> None:
    """Validate and convert legacy references."""


@legacy.command()
@click.option(
    "--path",
    help="the name to the legacy reference repository",
    required=True,
    type=click.Path(exists=True, file_okay=False, path_type=Path),
)
def precache(path: Path) -> None:
    """Pre-cache all accessions in a legacy reference repository."""
    ncbi = NCBIClient(path / ".migration_cache", False)

    buffer = []

    for otu_ in iter_legacy_otus(path / "src"):
        for isolate in otu_["isolates"]:
            for sequence in isolate["sequences"]:
                buffer.append(sequence["accession"])

                if len(buffer) > 450:
                    ncbi.fetch_genbank_records(buffer)
                    buffer = []

    if buffer:
        ncbi.fetch_genbank_records(buffer)


@legacy.command(name="format")
@click.option(
    "--path",
    help="the name to the legacy reference repository",
    required=True,
    type=click.Path(exists=True, file_okay=False, path_type=Path),
)
def reformat(path: Path) -> None:
    """Format a legacy reference repository.

    Re-formats every JSON file in a legacy reference repository to have a consistent
    format.
    """
    src_path = path / "src"

    for json_path in glob.glob(str(src_path / "**/*.json"), recursive=True):
        format_json(Path(json_path))


@legacy.command()
@debug_option
@click.option(
    "--fix",
    is_flag=True,
    help="attempt to fix errors in-place",
)
@click.option(
    "--limit",
    default=0,
    help="exit if this many otus have errors",
)
@click.option(
    "--no-ok",
    is_flag=True,
    help="don't print anything if an otu is valid",
)
@click.option(
    "--path",
    required=True,
    type=click.Path(exists=True, file_okay=False, path_type=Path),
    help="the path to a legacy reference directory",
)
def validate(debug: bool, fix: bool, limit: int, no_ok: bool, path: Path) -> None:
    """Validate a legacy reference repository."""
    configure_logger(debug)

    validate_legacy_repo(fix, limit, no_ok, path)


@entry.command()
@click.option(
    "-i",
    "--indent",
    is_flag=True,
    help="auto-indent the output JSON file",
)
@click.option(
    "-V",
    "--version",
    default="",
    type=str,
    help="a version string to include in the reference.json file",
)
@click.option(
    "-o",
    "--output-path",
    required=True,
    type=click.Path(exists=False, file_okay=True, path_type=Path),
    help="the path to write the reference.json file to",
)
@click.option(
    "-p",
    "--path",
    required=True,
    type=click.Path(exists=True, file_okay=False, path_type=Path),
    help="the path to the reference repository",
)
def build(output_path: Path, path: Path, indent: bool, version: str) -> None:
    """Build a Virtool reference.json file from a reference repository."""
    build_json(indent, output_path, path, version)

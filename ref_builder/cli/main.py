"""Command-line interface for reference builder."""

import glob
import sys
from pathlib import Path
from uuid import UUID

import click
from click import Context

from ref_builder.build import build_json
from ref_builder.cli.validate import validate_no_duplicate_accessions
from ref_builder.console import print_otu, print_otu_list
from ref_builder.legacy.convert import convert_legacy_repo
from ref_builder.legacy.utils import iter_legacy_otus
from ref_builder.legacy.validate import validate_legacy_repo
from ref_builder.logs import configure_logger
from ref_builder.ncbi.client import NCBIClient
from ref_builder.options import (
    ignore_cache_option,
    legacy_repo_path_option,
    path_option,
)
from ref_builder.otu.create import create_otu
from ref_builder.otu.update import (
    add_and_name_isolate,
    add_genbank_isolate,
    add_unnamed_isolate,
    auto_update_otu,
    exclude_accessions_from_otu,
    promote_otu_accessions,
    set_representative_isolate,
)
from ref_builder.otu.utils import RefSeqConflictError
from ref_builder.repo import Repo
from ref_builder.utils import DataType, IsolateName, IsolateNameType, format_json

PRECACHE_BUFFER_SIZE = 450
"""The number of Genbank records to fetch in a single request when pre-caching."""


@click.group()
@click.option("--debug", is_flag=True, help="Show debug logs")
def entry(debug: bool) -> None:
    """Build and maintain reference sets of pathogen genome sequences."""
    configure_logger(debug)


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


@entry.group()
def otu() -> None:
    """Manage OTUs."""


@otu.command(name="create")
@click.argument("TAXID", type=int)
@click.argument(
    "accessions_",
    callback=validate_no_duplicate_accessions,
    metavar="ACCESSIONS",
    nargs=-1,
    type=str,
    required=True,
)
@click.option("--acronym", type=str, default="")
@ignore_cache_option
@path_option
def otu_create(
    accessions_: list[str],
    acronym: str,
    ignore_cache: bool,
    path: Path,
    taxid: int,
) -> None:
    """Create a new OTU.

    OTUs are created from a list of accessions and a taxonomy ID. The associated Genbank
    records are fetched and used to create the first isolate and a plan.
    """
    repo = Repo(path)

    if len(accessions_) != len(set(accessions_)):
        click.echo("Duplicate accessions were provided.", err=True)
        sys.exit(1)

    try:
        create_otu(
            repo,
            taxid,
            accessions_,
            acronym=acronym,
            ignore_cache=ignore_cache,
        )
    except ValueError as e:
        click.echo(e, err=True)
        sys.exit(1)


@otu.command(name="get")
@click.argument("IDENTIFIER", type=str)
@path_option
def otu_get(identifier: str, path: Path) -> None:
    """Get an OTU by its unique ID or taxonomy ID."""
    try:
        identifier = int(identifier)
        otu_ = Repo(path).get_otu_by_taxid(identifier)
    except ValueError:
        otu_ = Repo(path).get_otu(UUID(identifier))

    if otu_ is None:
        click.echo("OTU not found.", err=True)
        sys.exit(1)

    print_otu(otu_)


@otu.command(name="list")
@path_option
def otu_list(path: Path) -> None:
    """List all OTUs in the repository."""
    print_otu_list(Repo(path).iter_minimal_otus())


@otu.group(invoke_without_command=True)
@click.argument("TAXID", type=int)
@path_option
@click.pass_context
def update(ctx: Context, path: Path, taxid: int) -> None:
    """Update the specified OTU with new data."""
    repo = Repo(path)

    ctx.ensure_object(dict)
    ctx.obj = {
        "REPO": repo,
        "TAXID": taxid,
    }

    if not repo.get_otu_id_by_taxid(taxid):
        click.echo(f"OTU {taxid} not found.", err=True)
        sys.exit(1)


@update.command(name="automatic")  # type: ignore
@ignore_cache_option
@click.pass_context
def otu_auto_update(ctx: Context, ignore_cache: bool) -> None:
    """Automatically update an OTU with the latest data from NCBI."""
    repo = ctx.obj["REPO"]
    taxid = ctx.obj["TAXID"]

    otu_ = repo.get_otu_by_taxid(taxid)

    if otu_ is None:
        click.echo(f"OTU not found for Taxonomy ID {taxid}.", err=True)
        click.echo(f'Run "virtool otu create {taxid} --autofill" instead.')
        sys.exit(1)

    auto_update_otu(repo, otu_, ignore_cache=ignore_cache)


@update.command(name="promote")  # type: ignore
@ignore_cache_option
@click.pass_context
def otu_promote_accessions(
    ctx: Context,
    ignore_cache: bool,
) -> None:
    """Promote all RefSeq accessions within this OTU."""
    repo = ctx.obj["REPO"]

    otu_ = repo.get_otu_by_taxid(ctx.obj["TAXID"])

    promote_otu_accessions(repo, otu_, ignore_cache)


@update.command(name="isolate")  # type: ignore
@click.argument(
    "accessions_",
    callback=validate_no_duplicate_accessions,
    metavar="ACCESSIONS",
    nargs=-1,
    type=str,
    required=True,
)
@click.option(
    "--unnamed",
    is_flag=True,
    default=False,
    help="ignore isolate names in Genbank sources",
)
@click.option(
    "--name",
    type=(IsolateNameType, str),
    help='an overriding name for the isolate, e.g. "isolate ARWV1"',
)
@ignore_cache_option
@click.pass_context
def isolate_create(
    ctx: Context,
    accessions_: list[str],
    ignore_cache: bool,
    name: tuple[IsolateNameType, str] | None,
    unnamed: bool,
) -> None:
    """Create a new isolate using the given accessions."""
    repo = ctx.obj["REPO"]
    taxid = ctx.obj["TAXID"]

    otu_ = repo.get_otu_by_taxid(taxid)
    if otu_ is None:
        click.echo(f"OTU {taxid} not found.", err=True)
        sys.exit(1)

    if unnamed:
        add_unnamed_isolate(
            repo,
            otu_,
            accessions_,
            ignore_cache=ignore_cache,
        )

    if name is not None:
        isolate_name_type, isolate_name_value = name
        isolate_name_ = IsolateName(type=isolate_name_type, value=isolate_name_value)

        try:
            add_and_name_isolate(
                repo,
                otu_,
                accessions_,
                ignore_cache=ignore_cache,
                isolate_name=isolate_name_,
            )
        except RefSeqConflictError as err:
            click.echo(
                f"{err.isolate_name} already exists, but RefSeq items may be "
                "promotable.",
            )
            sys.exit(1)

    try:
        add_genbank_isolate(
            repo,
            otu_,
            accessions_,
            ignore_cache=ignore_cache,
        )
    except RefSeqConflictError as err:
        click.echo(
            f"{err.isolate_name} already exists, but RefSeq items may be promotable,",
        )
        sys.exit(1)

    except ValueError as e:
        click.echo(e, err=True)
        sys.exit(1)


@update.command(name="exclude")  # type: ignore
@click.argument(
    "accessions_",
    metavar="ACCESSIONS",
    nargs=-1,
    type=str,
    required=True,
)
@click.pass_context
def accession_exclude(
    ctx: Context,
    accessions_: list[str],
) -> None:
    """Exclude the given accessions from this OTU.

    Any duplicate accessions will be ignored. Only one exclusion per unique accession
    will be made.
    """
    repo = ctx.obj["REPO"]
    taxid = ctx.obj["TAXID"]

    otu_ = repo.get_otu_by_taxid(taxid)

    if otu_ is None:
        click.echo(f"OTU {taxid} not found.", err=True)
        sys.exit(1)

    exclude_accessions_from_otu(repo, otu_, accessions_)


@update.command(name="default")  # type: ignore
@click.argument("ISOLATE_KEY", type=str)
@click.pass_context
def otu_set_representative_isolate(
    ctx: Context,
    isolate_key: str,
) -> None:
    """Update the OTU with a new representative isolate."""
    taxid = ctx.obj["TAXID"]

    repo = ctx.obj["REPO"]

    otu_ = repo.get_otu_by_taxid(taxid)

    if otu_ is None:
        click.echo(f"OTU not found for Taxonomy ID {taxid}.", err=True)
        sys.exit(1)

    try:
        isolate_id = UUID(isolate_key)
    except ValueError:
        parts = isolate_key.split(" ")
        try:
            isolate_name = IsolateName(IsolateNameType(parts[0].lower()), parts[1])
        except ValueError:
            click.echo(f'Error: "{isolate_key}" is not a valid isolate name.', err=True)
            sys.exit(1)

        isolate_id = otu_.get_isolate_id_by_name(isolate_name)

    if isolate_id is None:
        click.echo("Isolate could not be found in this OTU.", err=True)
        sys.exit(1)

    set_representative_isolate(repo, otu_, isolate_id)


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

import sys
from pathlib import Path
from uuid import UUID

import click
from click import Context
import structlog

from ref_builder.cli.validate import validate_no_duplicate_accessions
from ref_builder.console import print_otu, print_otu_list
from ref_builder.options import ignore_cache_option, path_option
from ref_builder.otu.create import create_otu
from ref_builder.otu.isolate import (
    add_unnamed_isolate,
    add_and_name_isolate,
    add_genbank_isolate,
)
from ref_builder.otu.modify import (
    exclude_accessions_from_otu,
    set_representative_isolate,
    add_segments_to_plan,
    resize_monopartite_plan,
)
from ref_builder.otu.update import auto_update_otu, promote_otu_accessions
from ref_builder.otu.utils import RefSeqConflictError
from ref_builder.plan import SegmentRule, SegmentName
from ref_builder.repo import Repo
from ref_builder.utils import IsolateNameType, IsolateName

logger = structlog.get_logger()


@click.group(name="otu")
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

    if ctx.invoked_subcommand is None:
        click.echo(ctx.get_help())


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


@update.group()
@click.pass_context
def plan(ctx: Context) -> None:
    """Add to and replace isolate plans for this OTU."""


@plan.command(name="extend")
@click.argument(
    "accessions_",
    metavar="ACCESSIONS",
    nargs=-1,
    type=str,
    required=True,
)
@click.option(
    "--optional",
    is_flag=True,
    help="Set additional segments as fully optional",
)
@click.pass_context
@ignore_cache_option
def plan_extend_segment_list(
    ctx: Context,
    accessions_: list[str],
    optional: bool,
    ignore_cache: bool,
) -> None:
    """Add recommended or optional segments to the OTU plan."""
    repo = ctx.obj["REPO"]
    taxid = ctx.obj["TAXID"]

    otu_ = repo.get_otu_by_taxid(taxid)
    if otu_ is None:
        click.echo(f"OTU {taxid} not found.", err=True)
        sys.exit(1)

    try:
        add_segments_to_plan(
            repo,
            otu_,
            rule=SegmentRule.OPTIONAL if optional else SegmentRule.RECOMMENDED,
            accessions=accessions_,
            ignore_cache=ignore_cache,
        )
    except ValueError as e:
        click.echo(e, err=True)


@plan.command(name="multipartite")
@click.argument(
    "accessions_",
    metavar="ACCESSIONS",
    nargs=-1,
    type=str,
    required=True,
)
@click.option(
    "--name",
    type=(str, str),
    help='a name for the originalsegment, e.g. "RNA A"',
)
@click.option(
    "--optional",
    is_flag=True,
    help="Set additional segments as fully optional",
)
@click.pass_context
@ignore_cache_option
def plan_resize_monopartite(
    ctx: Context,
    name: tuple[str, str],
    accessions_: list[str],
    optional: bool,
    ignore_cache: bool,
) -> None:
    """Change a monopartite OTU plan to a multipartite plan."""
    repo = ctx.obj["REPO"]
    taxid = ctx.obj["TAXID"]

    prefix, key = name

    otu_ = repo.get_otu_by_taxid(taxid)
    if otu_ is None:
        click.echo(f"OTU {taxid} not found.", err=True)
        sys.exit(1)

    try:
        resize_monopartite_plan(
            repo,
            otu_,
            name=SegmentName(prefix=prefix, key=key),
            rule=SegmentRule.OPTIONAL if optional else SegmentRule.RECOMMENDED,
            accessions=accessions_,
            ignore_cache=ignore_cache,
        )
    except ValueError as e:
        click.echo(e, err=True)

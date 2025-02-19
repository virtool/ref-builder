import sys
from pathlib import Path
from uuid import UUID

import click

from ref_builder.cli.utils import pass_repo
from ref_builder.cli.validate import validate_no_duplicate_accessions
from ref_builder.options import ignore_cache_option, path_option
from ref_builder.otu.isolate import (
    add_and_name_isolate,
    add_genbank_isolate,
    add_unnamed_isolate,
)
from ref_builder.otu.modify import delete_isolate_from_otu
from ref_builder.otu.utils import RefSeqConflictError
from ref_builder.repo import Repo, locked_repo
from ref_builder.utils import IsolateName, IsolateNameType


@click.group(name="isolate")
@click.pass_context
@path_option
def isolate(ctx: click.Context, path: Path) -> None:
    """Manage isolates."""
    ctx.obj = ctx.with_resource(locked_repo(path))


@isolate.command(name="create")  # type: ignore
@click.argument(
    "accessions_",
    callback=validate_no_duplicate_accessions,
    metavar="ACCESSIONS",
    nargs=-1,
    type=str,
    required=True,
)
@click.option("--taxid", type=int, required=True)
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
@pass_repo
def isolate_create(
    repo: Repo,
    accessions_: list[str],
    ignore_cache: bool,
    name: tuple[IsolateNameType, str] | None,
    taxid: int,
    unnamed: bool,
) -> None:
    """Create a new isolate using the given accessions."""
    otu_ = repo.get_otu_by_taxid(taxid)

    if otu_ is None:
        click.echo(f"OTU {taxid} not found.", err=True)
        sys.exit(1)

    if unnamed:
        with repo.use_transaction():
            add_unnamed_isolate(
                repo,
                otu_,
                accessions_,
                ignore_cache=ignore_cache,
            )
        sys.exit(1)

    if name is not None:
        isolate_name_type, isolate_name_value = name
        isolate_name_ = IsolateName(type=isolate_name_type, value=isolate_name_value)

        try:
            with repo.use_transaction():
                add_and_name_isolate(
                    repo,
                    otu_,
                    accessions_,
                    ignore_cache=ignore_cache,
                    isolate_name=isolate_name_,
                )
            sys.exit(1)

        except RefSeqConflictError as e:
            click.echo(
                f"{e.isolate_name} already exists, but RefSeq items may be "
                "promotable.",
            )
            sys.exit(1)

    try:
        with repo.use_transaction():
            add_genbank_isolate(
                repo,
                otu_,
                accessions_,
                ignore_cache=ignore_cache,
            )
    except RefSeqConflictError as e:
        click.echo(
            f"{e.isolate_name} already exists, but RefSeq items may be promotable,",
        )
        sys.exit(1)


@isolate.command(name="delete")  # type: ignore
@click.option("--taxid", type=int, required=True)
@click.argument("ISOLATE_KEY", type=str)
@pass_repo
def isolate_delete(repo: Repo, isolate_key: str, taxid: int) -> None:
    """Delete isolate."""
    otu_ = repo.get_otu_by_taxid(taxid)

    if otu_ is None:
        click.echo(f"OTU {taxid} not found.", err=True)
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

    delete_isolate_from_otu(repo, otu_, isolate_id)

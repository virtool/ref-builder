import sys
from pathlib import Path

import click

from ref_builder.cli.validate import validate_no_duplicate_accessions
from ref_builder.options import ignore_cache_option
from ref_builder.otu.isolate import add_unnamed_isolate, add_and_name_isolate, add_genbank_isolate
from ref_builder.otu.utils import RefSeqConflictError
from ref_builder.repo import Repo
from ref_builder.utils import IsolateNameType, IsolateName


@click.group(name="isolate")
def isolate() -> None:
    """Manage isolates."""


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
def isolate_create(
    path: Path,
    accessions_: list[str],
    ignore_cache: bool,
    taxid: int,
    name: tuple[IsolateNameType, str] | None,
    unnamed: bool,
) -> None:
    """Create a new isolate using the given accessions."""
    repo = Repo(path)

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

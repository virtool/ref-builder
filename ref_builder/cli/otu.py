import sys
from pathlib import Path
from uuid import UUID

import click
import structlog

from ref_builder.cli.validate import validate_no_duplicate_accessions
from ref_builder.console import print_otu, print_otu_list
from ref_builder.options import ignore_cache_option, path_option
from ref_builder.otu.create import create_otu
from ref_builder.otu.modify import (
    add_segments_to_plan,
    allow_accessions_into_otu,
    exclude_accessions_from_otu,
    rename_plan_segment,
    set_representative_isolate,
)
from ref_builder.otu.update import (
    auto_update_otu,
    batch_update_repo,
    promote_otu_accessions,
)
from ref_builder.plan import SegmentName, SegmentRule
from ref_builder.repo import Repo
from ref_builder.utils import IsolateName, IsolateNameType

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


@otu.command(name="batch-update")
@ignore_cache_option
@path_option
def otu_batch_auto_update(path: Path, ignore_cache: bool) -> None:
    """Update all OTUs with the latest data from NCBI."""
    repo = Repo(path)

    batch_update_repo(repo, ignore_cache)


@otu.command(name="update")
@ignore_cache_option
@path_option
@click.argument("TAXID", type=int)
def otu_auto_update(path: Path, taxid: int, ignore_cache: bool) -> None:
    """Update an OTU with the latest data from NCBI."""
    repo = Repo(path)
    otu_ = repo.get_otu_by_taxid(taxid)

    if otu_ is None:
        click.echo(f"OTU not found for Taxonomy ID {taxid}.", err=True)
        click.echo(f'Run "virtool otu create {taxid} --autofill" instead.')
        sys.exit(1)

    auto_update_otu(repo, otu_, ignore_cache=ignore_cache)


@otu.command(name="promote")
@path_option
@ignore_cache_option
@click.argument("TAXID", type=int)
def otu_promote_accessions(
    path: Path,
    taxid: int,
    ignore_cache: bool,
) -> None:
    """Promote all RefSeq accessions within an OTU."""
    repo = Repo(path)
    otu_ = repo.get_otu_by_taxid(taxid)

    if otu_ is None:
        click.echo(f"OTU not found for Taxonomy ID {taxid}.", err=True)
        click.echo(f'Run "virtool otu create {taxid} --autofill" instead.')
        sys.exit(1)

    promote_otu_accessions(repo, otu_, ignore_cache)


@otu.command(name="exclude-accessions")  # type: ignore
@path_option
@ignore_cache_option
@click.argument("TAXID", type=int)
@click.argument(
    "accessions_",
    metavar="ACCESSIONS",
    nargs=-1,
    type=str,
    required=True,
)
def otu_exclude_accessions(
    path: Path,
    taxid: int,
    accessions_: list[str],
) -> None:
    """Exclude the given accessions from this OTU.

    Any duplicate accessions will be ignored. Only one exclusion per unique accession
    will be made.
    """
    repo = Repo(path)
    otu_ = repo.get_otu_by_taxid(taxid)

    if otu_ is None:
        click.echo(f"OTU {taxid} not found.", err=True)
        sys.exit(1)

    exclude_accessions_from_otu(repo, otu_, accessions_)


@otu.command(name="allow-accessions")  # type: ignore
@path_option
@ignore_cache_option
@click.argument("TAXID", type=int)
@click.argument(
    "accessions_",
    metavar="ACCESSIONS",
    nargs=-1,
    type=str,
    required=True,
)
def otu_allow_accessions(
    path: Path,
    taxid: int,
    accessions_: list[str],
) -> None:
    """Allow the given excluded accessions back into the OTU.

    Any duplicate accessions will be ignored. Only one exclusion per unique accession
    will be made.
    """
    repo = Repo(path)
    otu_ = repo.get_otu_by_taxid(taxid)

    if otu_ is None:
        click.echo(f"OTU {taxid} not found.", err=True)
        sys.exit(1)

    allow_accessions_into_otu(repo, otu_, set(accessions_))


@otu.command(name="set-default-isolate")  # type: ignore
@path_option
@click.argument("TAXID", type=int)
@click.argument("ISOLATE_KEY", type=str)
def otu_set_representative_isolate(
    path: Path,
    taxid: int,
    isolate_key: str,
) -> None:
    """Update the OTU with a new representative isolate."""
    repo = Repo(path)
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

    set_representative_isolate(repo, otu_, isolate_id)


@otu.command(name="extend-plan")
@click.argument("TAXID", type=int)
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
@path_option
@ignore_cache_option
def plan_extend_segment_list(
    path: Path,
    taxid: int,
    accessions_: list[str],
    optional: bool,
    ignore_cache: bool,
) -> None:
    """Add recommended or optional segments to the OTU plan."""
    repo = Repo(path)
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


@otu.command(name="rename-plan-segment")
@click.argument("TAXID", type=int)
@click.argument("SEGMENT_ID", type=str)
@click.argument(
    "segment_name_",
    metavar="SEGMENT_NAME",
    nargs=2,
    type=str,
    required=True,
)
@path_option
def plan_rename_segment(
    path: Path, taxid: int, segment_id: str, segment_name_: tuple[str, str]
) -> None:
    """Give an unnamed segment a name."""
    repo = Repo(path)
    otu_ = repo.get_otu_by_taxid(taxid)

    if otu_ is None:
        click.echo(f"OTU {taxid} not found.", err=True)
        sys.exit(1)

    try:
        rename_plan_segment(
            repo,
            otu_,
            segment_id=UUID(segment_id),
            segment_name=SegmentName(segment_name_[0], segment_name_[1]),
        )
    except ValueError as e:
        click.echo(e, err=True)

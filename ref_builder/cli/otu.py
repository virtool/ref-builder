import datetime
import sys
from pathlib import Path
from uuid import UUID

import click
import structlog

from ref_builder.cli.options import ignore_cache_option, path_option
from ref_builder.cli.utils import (
    get_otu_from_identifier,
    get_otu_isolate_ids_from_identifier,
    pass_repo,
)
from ref_builder.cli.validate import validate_no_duplicate_accessions
from ref_builder.console import (
    print_otu,
    print_otu_as_json,
    print_otu_event_log,
    print_otu_list,
)
from ref_builder.otu.create import create_otu_with_taxid, create_otu_without_taxid
from ref_builder.otu.modify import (
    add_segments_to_plan,
    allow_accessions_into_otu,
    exclude_accessions_from_otu,
    rename_plan_segment,
    set_representative_isolate,
)
from ref_builder.otu.promote import (
    promote_otu_accessions,
    upgrade_outdated_sequences_in_otu,
)
from ref_builder.otu.update import (
    auto_update_otu,
    batch_update_repo,
)
from ref_builder.plan import SegmentName, SegmentRule
from ref_builder.repo import Repo, locked_repo

logger = structlog.get_logger()


@click.group(name="otu")
@path_option
@click.pass_context
def otu(ctx: click.Context, path: Path) -> None:
    """Manage OTUs."""
    ctx.obj = ctx.with_resource(locked_repo(path))


@otu.command(name="create")
@click.option("--taxid", type=int)
@click.argument(
    "accessions_",
    callback=validate_no_duplicate_accessions,
    metavar="ACCESSIONS",
    nargs=-1,
    type=str,
    required=True,
)
@click.option("--acronym", type=str, default="")
@pass_repo
@ignore_cache_option
def otu_create(
    repo: Repo,
    accessions_: list[str],
    acronym: str,
    ignore_cache: bool,
    taxid: int,
) -> None:
    """Create a new OTU.

    OTUs are created from a list of accessions and a taxonomy ID. The associated Genbank
    records are fetched and used to create the first isolate and a plan.
    """
    if len(accessions_) != len(set(accessions_)):
        click.echo("Duplicate accessions were provided.", err=True)
        sys.exit(1)

    if taxid:
        try:
            otu_ = create_otu_with_taxid(
                repo,
                taxid,
                accessions_,
                acronym=acronym,
                ignore_cache=ignore_cache,
            )
        except ValueError as e:
            click.echo(e, err=True)
            sys.exit(1)

    else:
        try:
            otu_ = create_otu_without_taxid(
                repo,
                accessions_,
                acronym=acronym,
                ignore_cache=ignore_cache,
            )
        except ValueError as e:
            click.echo(e, err=True)
            sys.exit(1)

    if otu_ is None:
        click.echo("OTU was not created correctly.", err=True)
        sys.exit(1)

    print_otu(otu_)


@otu.command(name="get")
@click.argument("IDENTIFIER", type=str)
@click.option(
    "--as-json",
    "--json",
    metavar="JSON",
    is_flag=True,
    help="Output in JSON form",
)
@pass_repo
def otu_get(repo: Repo, identifier: str, as_json: bool) -> None:
    """Get an OTU by its unique ID or taxonomy ID.

    IDENTIFIER is a taxonomy ID or unique OTU ID (>8 characters)
    """
    otu_ = get_otu_from_identifier(repo, identifier)

    if as_json:
        print_otu_as_json(otu_)
    else:
        print_otu(otu_)


@otu.command(name="list")
@pass_repo
def otu_list(repo: Repo) -> None:
    """List all OTUs in the repository."""
    print_otu_list(repo.iter_minimal_otus())


@otu.command(name="list-events")
@click.argument("IDENTIFIER", type=str)
@pass_repo
def otu_event_logs(repo: Repo, identifier: str) -> None:
    """Print a log of OTU events to console.

    IDENTIFIER is a taxonomy ID or unique OTU ID (>8 characters)
    """
    otu_ = get_otu_from_identifier(repo, identifier)

    print_otu_event_log(list(repo.iter_otu_events(otu_.id)))


@otu.command(name="batch-update")
@click.option("--precache", is_flag=True, help="Precache all NCBI Nucleotide records.")
@click.option(
    "--precache-batch-size",
    type=int,
    default=250,
    help="Change precache batch size (default 250).",
)
@click.option(
    "--start-date",
    type=click.DateTime(["%Y-%m-%d", "%Y/%m/%d"]),
    help="Exclude records edited before this date",
)
@click.option(
    "--fetch-index-path",
    type=click.Path(
        exists=True,
        file_okay=True,
        path_type=Path,
    ),
    help="Input a file path to a fetch index file.",
)
@ignore_cache_option
@pass_repo
def otu_batch_auto_update(
    repo: Repo,
    precache_batch_size: int,
    precache: bool,
    ignore_cache: bool,
    start_date: datetime.date | None,
    fetch_index_path: Path | None,
) -> None:
    """Update all OTUs with the latest data from NCBI."""
    batch_update_repo(
        repo,
        start_date=start_date,
        fetch_index_path=fetch_index_path,
        chunk_size=precache_batch_size,
        precache_records=precache,
        ignore_cache=ignore_cache,
    )


@otu.command(name="update")
@click.argument("IDENTIFIER", type=str)
@click.option(
    "--start-date",
    type=click.DateTime(["%Y-%m-%d", "%Y/%m/%d"]),
    help="Exclude records edited before this date",
)
@click.option(
    "--fetch-index-path",
    type=click.Path(
        exists=True,
        file_okay=True,
        path_type=Path,
    ),
    help="Input a file path to a fetch index file.",
)
@ignore_cache_option
@pass_repo
def otu_auto_update(
    repo: Repo,
    identifier: str,
    start_date: datetime.date | None,
    fetch_index_path: Path | None,
    ignore_cache: bool,
) -> None:
    """Update an OTU with the latest data from NCBI."""
    otu_ = get_otu_from_identifier(repo, identifier)

    auto_update_otu(
        repo,
        otu_,
        fetch_index_path=fetch_index_path,
        start_date=start_date,
        ignore_cache=ignore_cache,
    )


@otu.command(name="promote")
@click.argument("IDENTIFIER", type=str)
@ignore_cache_option
@pass_repo
def otu_promote_accessions(repo: Repo, identifier: str, ignore_cache: bool) -> None:
    """Promote all RefSeq accessions within an OTU."""
    otu_ = get_otu_from_identifier(repo, identifier)

    promote_otu_accessions(repo, otu_, ignore_cache)


@otu.command(name="upgrade")
@click.argument("IDENTIFIER", type=str)
@click.option(
    "--start-date",
    type=click.DateTime(["%Y-%m-%d", "%Y/%m/%d"]),
    help="Exclude records modified before this date",
)
@ignore_cache_option
@pass_repo
def otu_upgrade_accessions(
    repo: Repo, identifier: str, start_date: datetime.date | None, ignore_cache: bool
) -> None:
    """Upgrade all outdated sequences within an OTU."""
    otu_ = get_otu_from_identifier(repo, identifier)

    upgrade_outdated_sequences_in_otu(
        repo, otu_, modification_date_start=start_date, ignore_cache=ignore_cache
    )


@otu.command(name="exclude-accessions")  # type: ignore
@click.argument("IDENTIFIER", type=str)
@click.argument(
    "accessions_",
    metavar="ACCESSIONS",
    nargs=-1,
    type=str,
    required=True,
)
@pass_repo
def otu_exclude_accessions(
    repo: Repo,
    identifier: str,
    accessions_: list[str],
) -> None:
    """Exclude the given accessions from this OTU.

    Any duplicate accessions will be ignored. Only one exclusion per unique accession
    will be made.

    IDENTIFIER is a taxonomy ID or unique OTU ID (>8 characters)
    """
    otu_ = get_otu_from_identifier(repo, identifier)

    exclude_accessions_from_otu(repo, otu_, accessions_)


@otu.command(name="allow-accessions")  # type: ignore
@click.argument("IDENTIFIER", type=str)
@click.argument(
    "accessions_",
    metavar="ACCESSIONS",
    nargs=-1,
    type=str,
    required=True,
)
@pass_repo
def otu_allow_accessions(
    repo: Repo,
    identifier: str,
    accessions_: list[str],
) -> None:
    """Allow the given excluded accessions back into the OTU.

    Any duplicate accessions will be ignored. Only one exclusion per unique accession
    will be made.

    IDENTIFIER is a taxonomy ID or unique OTU ID (>8 characters)
    """
    otu_ = get_otu_from_identifier(repo, identifier)

    allow_accessions_into_otu(repo, otu_, set(accessions_))


@otu.command(name="set-default-isolate")  # type: ignore
@click.argument("IDENTIFIER", type=str)
@click.option(
    "--isolate-id",
    "isolate_id_",
    type=str,
    required=True,
    help="The id of the isolate to set as default.",
)
@pass_repo
def otu_set_representative_isolate(
    repo: Repo,
    isolate_id_: str,
    identifier: str,
) -> None:
    """Update the OTU with a new representative isolate.

    IDENTIFIER is a taxonomy ID or unique OTU ID (>8 characters)
    """
    otu_ = get_otu_from_identifier(repo, identifier)

    otu_id, isolate_id = get_otu_isolate_ids_from_identifier(repo, isolate_id_)

    set_representative_isolate(repo, otu_, isolate_id)


@otu.command(name="extend-plan")
@click.argument("IDENTIFIER", type=str)
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
@ignore_cache_option
@pass_repo
def plan_extend_segment_list(
    repo: Repo,
    accessions_: list[str],
    ignore_cache: bool,
    identifier: str,
    optional: bool,
) -> None:
    """Add recommended or optional segments to the OTU plan.

    IDENTIFIER is a taxonomy ID or unique OTU ID (>8 characters)
    """
    otu_ = get_otu_from_identifier(repo, identifier)

    if not add_segments_to_plan(
        repo,
        otu_,
        rule=SegmentRule.OPTIONAL if optional else SegmentRule.RECOMMENDED,
        accessions=accessions_,
        ignore_cache=ignore_cache,
    ):
        sys.exit(1)


@otu.command(name="rename-plan-segment")
@click.argument("IDENTIFIER", type=str)
@click.option(
    "--segment-id",
    "segment_id_",
    type=str,
    required=True,
    help="The ID of the segment.",
)
@click.option(
    "--segment-name",
    "--name",
    "segment_name_",
    type=(str, str),
    required=True,
    help="A segment name, e.g. 'RNA B'",
)
@pass_repo
def plan_rename_segment(
    repo: Repo,
    segment_id_: str,
    segment_name_: tuple[str, str],
    identifier: str,
) -> None:
    """Give an unnamed segment a name.

    IDENTIFIER is a taxonomy ID or unique OTU ID (>8 characters)
    """
    otu_ = get_otu_from_identifier(repo, identifier)

    try:
        rename_plan_segment(
            repo,
            otu_,
            segment_id=UUID(segment_id_),
            segment_name=SegmentName(segment_name_[0], segment_name_[1]),
        )
    except ValueError as e:
        click.echo(e, err=True)

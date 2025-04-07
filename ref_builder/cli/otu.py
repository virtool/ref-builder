import datetime
import sys
from pathlib import Path
from uuid import UUID

import click
import structlog

from ref_builder.cli.utils import pass_repo, get_otu_from_identifier
from ref_builder.cli.validate import validate_no_duplicate_accessions
from ref_builder.console import (
    print_otu,
    print_otu_as_json,
    print_otu_list,
    print_otu_event_log,
)
from ref_builder.errors import PartialIDConflictError, InvalidInputError
from ref_builder.options import ignore_cache_option, path_option
from ref_builder.otu.create import create_otu_with_taxid, create_otu_without_taxid
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
)
from ref_builder.otu.promote import promote_otu_accessions
from ref_builder.plan import SegmentName, SegmentRule
from ref_builder.repo import Repo, locked_repo
from ref_builder.utils import IsolateName, IsolateNameType

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
    try:
        otu_id = UUID(identifier)
    except ValueError:
        otu_id = _get_otu_id_from_other_identifier(repo, identifier)

    if otu_id is None:
        click.echo("OTU not found.", err=True)
        sys.exit(1)

    otu_ = repo.get_otu(otu_id)

    if otu_ is None:
        click.echo("OTU not found.", err=True)
        sys.exit(1)

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
    try:
        otu_id = UUID(identifier)
    except ValueError:
        otu_id = _get_otu_id_from_other_identifier(repo, identifier)

    otu_ = repo.get_otu(otu_id)

    if otu_ is None:
        click.echo("OTU not found.", err=True)
        sys.exit(1)

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
    try:
        otu_id = UUID(identifier)
    except ValueError:
        otu_id = _get_otu_id_from_other_identifier(repo, identifier)

    otu_ = repo.get_otu(otu_id)

    if otu_ is None:
        click.echo("OTU not found.", err=True)
        sys.exit(1)

    allow_accessions_into_otu(repo, otu_, set(accessions_))


@otu.command(name="set-default-isolate")  # type: ignore
@click.argument("TAXID", type=int)
@click.argument("ISOLATE_KEY", type=str)
@pass_repo
def otu_set_representative_isolate(
    repo: Repo,
    isolate_key: str,
    taxid: int,
) -> None:
    """Update the OTU with a new representative isolate."""
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
@ignore_cache_option
@pass_repo
def plan_extend_segment_list(
    repo: Repo,
    accessions_: list[str],
    ignore_cache: bool,
    taxid: int,
    optional: bool,
) -> None:
    """Add recommended or optional segments to the OTU plan."""
    otu_ = repo.get_otu_by_taxid(taxid)

    if otu_ is None:
        click.echo(f"OTU {taxid} not found.", err=True)
        sys.exit(1)

    with repo.use_transaction():
        add_segments_to_plan(
            repo,
            otu_,
            rule=SegmentRule.OPTIONAL if optional else SegmentRule.RECOMMENDED,
            accessions=accessions_,
            ignore_cache=ignore_cache,
        )


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
@pass_repo
def plan_rename_segment(
    repo: Repo,
    segment_id: str,
    segment_name_: tuple[str, str],
    taxid: int,
) -> None:
    """Give an unnamed segment a name."""
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


def _get_otu_id_from_other_identifier(repo: Repo, identifier: str) -> UUID | None:
    """Return an OTU id from the repo if identifier matches a single OTU.
    Return None if no matching OTU is found, raise a ValueError if >1 OTU is found.
    """
    try:
        taxid = int(identifier)
    except ValueError:
        pass
    else:
        return repo.get_otu_id_by_taxid(taxid)

    try:
        return repo.get_otu_id_by_partial(identifier)

    except PartialIDConflictError:
        click.echo(
            "Partial ID too short to narrow down results.",
            err=True,
        )

    except InvalidInputError as e:
        click.echo(
            e,
            err=True,
        )

    return None

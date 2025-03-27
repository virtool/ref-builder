import sys
from pathlib import Path

import click

from ref_builder.cli.utils import pass_repo
from ref_builder.console import print_event, print_event_as_json
from ref_builder.options import path_option
from ref_builder.repo import Repo, locked_repo


@click.group(name="events")
@path_option
@click.pass_context
def events(ctx: click.Context, path: Path) -> None:
    """Read events."""
    ctx.obj = ctx.with_resource(locked_repo(path))


@events.command(name="get")
@click.argument("IDENTIFIER", type=int)
@click.option(
    "--as-json",
    "--json",
    metavar="JSON",
    is_flag=True,
    help="Output in JSON form",
)
@pass_repo
def event_get(repo: Repo, identifier: int, as_json: bool) -> None:
    """Get and view an event."""
    if (event := repo.get_event(identifier)) is None:
        click.echo("Event not found.", err=True)
        sys.exit(1)

    if as_json:
        print_event_as_json(event)
    else:
        print_event(event)

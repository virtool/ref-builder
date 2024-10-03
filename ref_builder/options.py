"""Common options for the Click-based CLI."""

from pathlib import Path

import click

ignore_cache_option = click.option(
    "--ignore-cache",
    is_flag=True,
    help="Ignore cached records",
)
"""A click option for disabling cached results in fetch operations."""

path_option = click.option(
    "--path",
    default=".",
    type=click.Path(exists=True, file_okay=False, path_type=Path),
    help="the path to the reference repository",
)
"""A click option for the path to the reference repository."""

import sys
from uuid import UUID

import click

from ref_builder.errors import PartialIDConflictError, InvalidInputError
from ref_builder.repo import Repo
from ref_builder.resources import RepoOTU

pass_repo = click.make_pass_decorator(Repo)


def get_otu_from_identifier(repo: Repo, identifier: str) -> RepoOTU:
    """Return an OTU from the given repo if identifier matches a single OTU."""
    try:
        otu_id = UUID(identifier)
    except ValueError:
        otu_id = _get_otu_id_from_other_identifier(repo, identifier)

    if otu_id is None:
        click.echo("OTU not found.", err=True)
        sys.exit(1)

    if (otu_ := repo.get_otu(otu_id)) is None:
        click.echo("OTU not found.", err=True)
        sys.exit(1)

    return otu_


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

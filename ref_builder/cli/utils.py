import click

from ref_builder.repo import Repo

pass_repo = click.make_pass_decorator(Repo)

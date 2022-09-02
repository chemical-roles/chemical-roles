# -*- coding: utf-8 -*-

"""CLI for Chemical Roles."""

import click

from .chebi import curate_chebi
from .mesh import curate_mesh


@click.group()
def curate():
    """Run the curation CLI."""


curate.add_command(curate_chebi)
curate.add_command(curate_mesh)

if __name__ == "__main__":
    curate()

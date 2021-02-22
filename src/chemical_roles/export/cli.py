# -*- coding: utf-8 -*-

"""CLI for Chemical Roles exporters."""

import click
import pybel


@click.group()
def export():
    """Export the database."""


@export.command()
@click.argument('path')
def bel(path):
    """Write BEL export."""
    from .bel import get_bel
    graph = get_bel()
    pybel.dump(graph, path)


if __name__ == '__main__':
    export()

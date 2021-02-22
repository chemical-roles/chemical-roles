# -*- coding: utf-8 -*-

"""CLI for Chemical Roles exporters."""

import click


@click.group()
def export():
    """Export the database."""


@export.command()
def summary():
    """Rewrite readme and generate new export."""
    from .build import rewrite_repo_readme, write_export
    import seaborn as sns
    sns.set(font_scale=1.3, style='whitegrid')
    rewrite_repo_readme()
    write_export()


@export.command()
@click.argument('path')
def bel(path):
    """Write BEL export."""
    import pybel
    from .bel import get_bel
    graph = get_bel()
    pybel.dump(graph, path)


if __name__ == '__main__':
    export()

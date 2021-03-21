# -*- coding: utf-8 -*-

"""CLI for Chemical Roles exporters."""

import os

import click

from ..constants import DATA


@click.group()
def export():
    """Export the database."""


@export.command(name='all')
@click.pass_context
def export_all(ctx):
    """Export all."""
    ctx.invoke(summary)
    ctx.invoke(obo)
    ctx.invoke(bel)


directory_option = click.option('--directory', default=DATA)


@export.command()
def summary():
    """Rewrite readme and generate new export."""
    from .build import rewrite_repo_readme, write_export
    import seaborn as sns
    sns.set(font_scale=1.3, style='whitegrid')
    rewrite_repo_readme()
    write_export()


@export.command()
@directory_option
def bel(directory):
    """Write BEL export."""
    import pybel
    from .bel import get_bel
    graph = get_bel()
    pybel.dump(graph, os.path.join(directory, 'crog.bel.nodelink.json.gz'))


@export.command()
@directory_option
def obo(directory):
    """Write OBO export."""
    from .obo import get_obo
    o = get_obo()
    o.write_obo(os.path.join(directory, 'crog.obo'))
    o.write_obonet_gz(os.path.join(directory, 'crog.obonet.json.gz'))


if __name__ == '__main__':
    export()

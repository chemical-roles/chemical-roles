# -*- coding: utf-8 -*-

"""CLI for Chemical Roles exporters."""

import os
from pathlib import Path

import click

from ..constants import DATA


@click.group()
def export() -> None:
    """Export the database."""


@export.command(name="all")
@click.pass_context
def export_all(ctx: click.Context) -> None:
    """Export all."""
    ctx.invoke(summary)
    ctx.invoke(obo)
    ctx.invoke(bel)
    ctx.invoke(indra)


directory_option = click.option("--directory", default=DATA, type=Path)


@export.command()
def summary() -> None:
    """Rewrite readme and generate new export."""
    import seaborn as sns

    from .build import rewrite_repo_readme, write_export

    sns.set(font_scale=1.3, style="whitegrid")
    rewrite_repo_readme()
    write_export()


@export.command()
@directory_option
def bel(directory: Path) -> None:
    """Write BEL export."""
    import pybel

    from .bel import get_bel

    graph = get_bel()
    pybel.dump(graph, os.path.join(directory, "crog.bel.nodelink.json.gz"))


@export.command()
@directory_option
def indra(directory: Path) -> None:
    """Write INDRA export."""
    import pybel

    from .bel import get_bel

    graph = get_bel(use_inferred=False, add_evidence=False)
    pybel.to_indra_statements_json_file(
        graph, os.path.join(directory, "crog.indra.json"), sort_keys=True
    )


@export.command()
@directory_option
def obo(directory: Path) -> None:
    """Write OBO export."""
    from .obo import get_obo

    o = get_obo(use_inferred=False, use_sub_roles=False)
    o.write_obo(os.path.join(directory, "crog-slim.obo"))
    o.write_owl(os.path.join(directory, "crog-slim.owl"))

    o = get_obo(use_inferred=True, use_sub_roles=True)
    o.write_obo(os.path.join(directory, "crog.obo"))
    o.write_owl(os.path.join(directory, "crog.owl"))


if __name__ == "__main__":
    export()

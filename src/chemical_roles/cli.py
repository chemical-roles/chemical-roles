# -*- coding: utf-8 -*-

"""CLI for Chemical Roles."""

import click

from .curate.cli import curate
from .export.cli import export
from .lint import lint


@click.group()
def main():
    """Run the Chemical-Roles CLI."""


main.add_command(curate)
main.add_command(export)
main.add_command(lint)

if __name__ == "__main__":
    main()

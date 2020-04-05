# -*- coding: utf-8 -*-

"""Scripts to maintain the integrity of the curated data."""

import logging
import sys

import click

from utils import XREFS_PATH, get_single_mappings, get_xrefs_df, sort_xrefs_df


@click.group()
def main():
    """Run linters."""


@main.command()
def tabs():
    """Find entries missing tabs."""
    errors = set()
    with open(XREFS_PATH) as file:
        for i, line in enumerate(file, start=1):
            line = line.strip()
            if 7 != line.count('\t'):
                print(f'line {i}: Not enough entries - {line}')
                errors.add(i)
    if errors:
        click.secho('Document not clean.', fg='red', bold=True)
        sys.exit(1)


@main.command()
def mappings():
    """Find single mapped entries."""
    df = get_xrefs_df()

    idx = (df['target_db'] == 'pr') & (df['type'] == 'protein')
    errors = get_single_mappings(df, idx)
    if errors:
        click.secho('Some entries only mapped to Protein Ontology', fg='red', bold=True)
        sys.exit(1)

    idx = (df['target_db'] == 'go') & (df['type'] == 'protein complex')
    errors = get_single_mappings(df, idx)
    if errors:
        click.secho('Some complexes only mapped to Gene Ontology', fg='red', bold=True)
        sys.exit(1)


@main.command()
def sort():
    """Sort the entries."""
    sort_xrefs_df()


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    main()

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
    protein_pr_errors = get_single_mappings(df, idx)
    if protein_pr_errors:
        click.secho('Some entries only mapped to Protein Ontology', fg='red', bold=True)
        _p(protein_pr_errors)

    idx = (df['target_db'] == 'go') & (df['type'] == 'protein complex')
    go_complex_errors = get_single_mappings(df, idx)
    if go_complex_errors:
        click.secho('Some complexes only mapped to Gene Ontology', fg='red', bold=True)
        _p(go_complex_errors)

    idx = (df['target_db'] == 'mesh')
    mesh_errors = get_single_mappings(df, idx)
    if mesh_errors:
        click.secho('Some roles only mapped to MeSH', fg='red', bold=True)
        _p(mesh_errors)

    idx = (df['type'] == 'molecular function')
    mf_errors = get_single_mappings(df, idx)
    if mf_errors:
        click.secho('Some roles only mapped to molecular function', fg='red', bold=True)
        _p(mf_errors)

    if any([
        protein_pr_errors,
        go_complex_errors,
        mesh_errors,
        mf_errors,
    ]):
        click.secho('The job is not yet done...', fg='red')


def _p(errors):
    m = max(len(k) for _, _, k in errors)
    for (k_db, k_id, k), vs in errors.items():
        click.echo(f'{k_db}:{k_id} ! {k:{m}} for role{"" if len(vs) == 1 else "s"} {", ".join(vs)}')


@main.command()
def sort():
    """Sort the entries."""
    sort_xrefs_df()


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    main()

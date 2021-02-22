# -*- coding: utf-8 -*-

"""Scripts to maintain the integrity of the curated data."""

import sys

import bioregistry
import click
from more_click import verbose_option

from .resources import XREFS_PATH, get_xrefs_df
from .utils import get_single_mappings, sort_xrefs_df


@click.group()
def lint():
    """Run linters."""


@lint.command()
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


@lint.command()
@verbose_option
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


@lint.command()
def sort():
    """Sort the entries."""
    sort_xrefs_df()


@lint.command()
def validate():
    """Validate identifiers."""
    df = get_xrefs_df()
    for i, (prefix, identifier) in df[['source_db', 'source_id']].iterrows():
        norm_prefix = bioregistry.normalize_prefix(prefix)
        if prefix != norm_prefix:
            raise ValueError(f'invalid source prefix: {prefix} should be {norm_prefix}')
        if not bioregistry.validate(prefix, identifier):
            raise ValueError(
                f'[line {i}] Invalid source curie: {prefix}:{identifier} for pattern {bioregistry.get_pattern(prefix)}',
            )
    for i, (prefix, identifier) in df[['target_db', 'target_id']].iterrows():
        norm_prefix = bioregistry.normalize_prefix(prefix)
        if prefix != norm_prefix:
            raise ValueError(f'invalid target prefix: {prefix} should be {norm_prefix}')
        if not bioregistry.validate(prefix, identifier):
            raise ValueError(
                f'[line {i}] Invalid target curie: {prefix}:{identifier} for pattern {bioregistry.get_pattern(prefix)}',
            )


if __name__ == '__main__':
    lint()

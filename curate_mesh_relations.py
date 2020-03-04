# -*- coding: utf-8 -*-

"""A script to help curate MeSH relations and infer new ones."""

import os

import click
import pandas as pd
from pyobo.sources import mesh

HERE = os.path.abspath(os.path.dirname(__file__))
MESH_XREFS_PATH = os.path.join(HERE, 'resources', 'xrefs.tsv')

BLACKLIST = {
    'D004791',  # Enzyme
}

SUFFIXES = [
    'inhibitor',
    'activator',
    'agonist',
    'antagonist',
    'modulator',
    'suppressor',
    'deactivator',
    'drug',
    'agent',
]
SUFFIXES.extend([
    f'{suffix}s'
    for suffix in SUFFIXES
])


@click.command()
def main():
    """Run the MeSH curation pipeline."""
    xrefs_df = pd.read_csv(MESH_XREFS_PATH, sep='\t')
    mesh_xrefs_df = xrefs_df[xrefs_df['source_db'] == 'mesh']
    curated = set(mesh_xrefs_df['source_id'])

    mesh_obo = mesh.get_obo()
    terms = {
        term.curie: (term, suffix)
        for term in mesh_obo
        if term.identifier not in curated and term.identifier not in BLACKLIST
        for suffix in SUFFIXES
        if term.name.lower().endswith(suffix)
    }

    for i, (curie, (term, suffix)) in enumerate(sorted(terms.items()), start=1):
        print('mesh', term.identifier, term.name, suffix.rstrip('s'), '?', '?', '?', '?', '?', sep='\t')


if __name__ == '__main__':
    main()

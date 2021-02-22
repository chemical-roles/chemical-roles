# -*- coding: utf-8 -*-

"""A script to help curate MeSH relations and infer new ones."""

from typing import Optional, TextIO

import click
import pyobo
from more_click import verbose_option
from tqdm import tqdm

from ..resources import UNCURATED_MESH_PATH, get_xrefs_df
from ..utils import SUFFIXES, yield_gilda

MESH_BLACKLIST = {
    'D004791',  # Enzyme
    'D000074389',  # Therapeutic Index, Drug
    'D000075203',  # Contraindications, Drug
    'D000076742',  # Synthetic Drugs
    'D000078742',  # Substandard Drugs
    'D000078903',  # Catalog, Drug
    'D004305',  # Dose-Response Relationship, Drug
    'D004366',  # Nonprescription Drugs
    'D007202',  # Indicators and Reagents
    'D007880',  # Legislation, Drug
    'D011355',  # Prodrugs
    'D013287',  # Street Drugs
    'D015198',  # Designer Drugs
    'D016147',  # Genes, Tumor Suppressor
    'D016153',  # Genes, Suppressor
    'D019155',  # Veterinary Drugs
}


@click.command(name='mesh')
@verbose_option
@click.option('--show-ungrounded', is_flag=True)
@click.option('--output', type=click.File('w'), default=UNCURATED_MESH_PATH)
def curate_mesh(show_ungrounded: bool, output: Optional[TextIO]):
    """Run the MeSH curation pipeline."""
    xrefs_df = get_xrefs_df()
    mesh_xrefs_df = xrefs_df[xrefs_df['source_db'] == 'mesh']
    curated_mesh_ids = set(mesh_xrefs_df['source_id'])

    terms = {
        identifier: (name, name[:-len(suffix)], suffix.strip('s'))
        for identifier, name in pyobo.get_id_name_mapping('mesh').items()
        if identifier not in curated_mesh_ids and identifier not in MESH_BLACKLIST
        for suffix in SUFFIXES
        if name.lower().endswith(suffix)
    }

    it = sorted(terms.items(), key=lambda t: t[1][0])
    it = tqdm(it, desc='making MeSH curation sheet')
    for _, (identifier, (name, search_text, suffix)) in enumerate(it, start=1):
        for row in yield_gilda('mesh', identifier, name, suffix, search_text, show_ungrounded or output is not None):
            print(*row, sep='\t', file=output)


if __name__ == '__main__':
    curate_mesh()

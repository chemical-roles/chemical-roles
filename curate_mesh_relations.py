# -*- coding: utf-8 -*-

"""A script to help curate MeSH relations and infer new ones."""

import click
import pyobo

from utils import SUFFIXES, get_xrefs_df

BLACKLIST = {
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


@click.command()
def main():
    """Run the MeSH curation pipeline."""
    xrefs_df = get_xrefs_df()
    mesh_xrefs_df = xrefs_df[xrefs_df['source_db'] == 'mesh']
    curated_mesh_ids = set(mesh_xrefs_df['source_id'])

    terms = {
        identifier: (name, suffix.strip('s'))
        for identifier, name in pyobo.get_id_name_mapping('mesh').items()
        if identifier not in curated_mesh_ids and identifier not in BLACKLIST
        for suffix in SUFFIXES
        if name.lower().endswith(suffix)
    }

    for i, (identifier, (name, suffix)) in enumerate(sorted(terms.items(), key=lambda t: t[1][0]), start=1):
        print('mesh', identifier, name, suffix, '?', '?', '?', '?', sep='\t')


if __name__ == '__main__':
    main()

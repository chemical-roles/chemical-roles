# -*- coding: utf-8 -*-

"""Chemical Roles resources."""

import os

import pandas as pd

HERE = os.path.abspath(os.path.dirname(__file__))

RECLASSIFICATION_PATH = os.path.join(HERE, 'reclassification.tsv')
UNCURATED_CHEBI_PATH = os.path.join(HERE, 'uncurated_chebi.tsv')
UNCURATED_MESH_PATH = os.path.join(HERE, 'uncurated_mesh.tsv')
XREFS_PATH = os.path.join(HERE, 'xrefs.tsv')
IRRELEVANT_ROLES_PATH = os.path.join(HERE, 'irrelevant_roles.tsv')
BLACKLIST_ROLES_PATH = os.path.join(HERE, 'blacklist.tsv')


def get_xrefs_df() -> pd.DataFrame:
    """Get xrefs.tsv."""
    return pd.read_csv(XREFS_PATH, sep='\t', comment='#', dtype=str)


def get_irrelevant_roles_df() -> pd.DataFrame:
    """Get irrelevant roles."""
    return pd.read_csv(IRRELEVANT_ROLES_PATH, sep='\t', dtype=str, comment='#', skip_blank_lines=True)


def get_blacklist_roles_df() -> pd.DataFrame:
    """Get roles blacklisted (should not be curated)."""
    return pd.read_csv(BLACKLIST_ROLES_PATH, sep='\t', dtype=str)

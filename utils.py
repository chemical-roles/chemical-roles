# -*- coding: utf-8 -*-

"""Chemical relation curation utilities."""

import logging
import os
from collections import defaultdict
from typing import Mapping, Set

import pandas as pd
import requests

logger = logging.getLogger(__name__)

HERE = os.path.abspath(os.path.dirname(__file__))
RESOURCES_DIRECTORY = os.path.join(HERE, 'resources')
EXPORT_DIRECTORY = os.path.join(HERE, 'export')

#: Path to xrefs.tsv
XREFS_PATH = os.path.join(RESOURCES_DIRECTORY, 'xrefs.tsv')
XREFS_COLUMNS = [
    'source_db', 'source_id', 'source_name',
    'modulation',
    'target_type', 'target_db', 'target_id', 'target_name',
]

IRRELEVANT_ROLES_PATH = os.path.join(RESOURCES_DIRECTORY, 'irrelevant_roles.tsv')
IRRELEVANT_ROLES_COLUMNS = ['database', 'identifier', 'name']

BLACKLIST_ROLES_PATH = os.path.join(RESOURCES_DIRECTORY, 'blacklist.tsv')
BLACKLIST_ROLES_COLUMNS = ['database', 'identifier', 'name']

SUFFIXES = [
    'inhibitor',
    'deactivator',
    'activator',
    'antagonist',
    'agonist',
    'modulator',
    'suppressor',
    'drug',
    'agent',
]
SUFFIXES.extend([
    f'{suffix}s'
    for suffix in SUFFIXES
])


def get_xrefs_df() -> pd.DataFrame:
    """Get xrefs.tsv."""
    return pd.read_csv(XREFS_PATH, sep='\t', comment='#', dtype=str)


def get_irrelevant_roles_df() -> pd.DataFrame:
    """Get irrelevant roles."""
    return pd.read_csv(IRRELEVANT_ROLES_PATH, sep='\t', dtype=str, comment='#', skip_blank_lines=True)


def get_blacklist_roles_df() -> pd.DataFrame:
    """Get roles blacklisted (should not be curated)."""
    return pd.read_csv(BLACKLIST_ROLES_PATH, sep='\t', dtype=str)


def sort_xrefs_df() -> None:
    """Sort xrefs.tsv."""
    df = get_xrefs_df()
    df = df.sort_values(['source_db', 'source_name', 'modulation'])
    df = df.drop_duplicates()
    df.to_csv(XREFS_PATH, index=False, sep='\t')


#: The base URL of the GILDA service
GILDA_URL = 'http://grounding.indra.bio'


def post_gilda(text: str, url: str = GILDA_URL) -> requests.Response:
    """Send text to GILDA."""
    return requests.post(f'{url}/ground', json={'text': text})


def get_single_mappings(df: pd.DataFrame, idx) -> Mapping[str, Set[str]]:
    """Get ChEBI identifiers that are only mapped to one thing based on slicing the dataframe on the given index."""
    errors = defaultdict(set)
    chebi_ids = sorted(df.loc[idx, 'source_id'].unique())
    for chebi_id, sdf in df[df['source_id'].isin(chebi_ids)].groupby(['source_id']):
        if 1 == len(sdf.index):
            target_db, target_id, target_name = list(sdf[['target_db', 'target_id', 'target_name']].values)[0]
            errors[target_db, target_id, target_name].add(chebi_id)
    return dict(errors)

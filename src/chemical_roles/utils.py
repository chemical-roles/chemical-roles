# -*- coding: utf-8 -*-

"""Chemical relation curation utilities."""

import logging
from collections import defaultdict
from typing import Iterable, Mapping, Set, Tuple

import pandas as pd
import requests

from .resources import XREFS_PATH, get_xrefs_df

logger = logging.getLogger(__name__)

XREFS_COLUMNS = [
    'source_db', 'source_id', 'source_name',
    'modulation',
    'target_type', 'target_db', 'target_id', 'target_name',
]

IRRELEVANT_ROLES_COLUMNS = ['database', 'identifier', 'name']

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


GildaTuple = Tuple[str, str, str, str, str, str, str, str]


def yield_gilda(
    source_db: str,
    identifier: str,
    name: str,
    suffix: str,
    search_text: str,
    show_missing: bool,
) -> Iterable[GildaTuple]:
    """Yield results from gilda."""
    results = post_gilda(search_text).json()
    if results:
        for result in results:
            term = result["term"]
            term_db = term['db'].lower()
            term_id = term['id']
            if term_id.lower().startswith(f'{term_db}:'):
                term_id = term_id[len(f'{term_db}:'):]

            if term_db == source_db and term_id == identifier:
                continue
            yield (
                source_db, identifier, name,
                suffix or '?',
                '?', term_db, term_id, term['entry_name'],
            )
    elif show_missing:
        yield source_db, identifier, name, suffix or '?', '?', '?', '?', '?'


def get_single_mappings(df: pd.DataFrame, idx) -> Mapping[str, Set[str]]:
    """Get ChEBI identifiers that are only mapped to one thing based on slicing the dataframe on the given index."""
    errors = defaultdict(set)
    chebi_ids = sorted(df.loc[idx, 'source_id'].unique())
    for chebi_id, sdf in df[df['source_id'].isin(chebi_ids)].groupby(['source_id']):
        if 1 == len(sdf.index):
            target_db, target_id, target_name = list(sdf[['target_db', 'target_id', 'target_name']].values)[0]
            errors[target_db, target_id, target_name].add(chebi_id)
    return dict(errors)

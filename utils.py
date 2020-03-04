# -*- coding: utf-8 -*-

"""Chemical relation curation utilities."""

import os

import pandas as pd
import requests

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


def get_xrefs_df() -> pd.DataFrame:
    """Get xrefs.tsv."""
    return pd.read_csv(XREFS_PATH, sep='\t', comment='#', dtype={'db_id': str})


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

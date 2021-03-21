# -*- coding: utf-8 -*-

"""Export utilities."""

import itertools as itt
import logging
from collections import defaultdict
from functools import lru_cache
from typing import Iterable, List, Mapping, Tuple

import networkx as nx
import pandas as pd
import pyobo
from protmapper import uniprot_client
from protmapper.api import hgnc_id_to_up, hgnc_name_to_id
from pyobo.sources import expasy
from pyobo.sources.chebi import get_chebi_role_to_children
from pyobo.struct import has_member
from tqdm import tqdm

from chemical_roles.resources import get_xrefs_df
from chemical_roles.utils import XREFS_COLUMNS

logger = logging.getLogger(__name__)


@lru_cache(maxsize=4)
def get_relations_df(use_sub_roles: bool = False, use_inferred: bool = True) -> pd.DataFrame:
    """Assemble the relations dataframe."""
    xrefs_df = get_xrefs_df()
    if not use_inferred:
        return xrefs_df

    famplex_id_to_members = _get_famplex()

    logger.info('getting enzyme classes')
    expasy_graph, ec_code_to_children = get_expasy_closure()
    logger.info('getting ec2go')
    ec2go = expasy.get_ec2go()
    logger.info('ec2go has %d elements', len(ec2go))

    rows = list(xrefs_df.values)
    x = defaultdict(list)
    it = tqdm(
        rows,
        total=len(xrefs_df.index),
        desc='inferring over target hierarchies',
    )
    non_chebi_counter = 0
    for source_db, source_id, _source_name, modulation, target_type, target_db, target_id, target_name in it:
        if source_db != 'chebi':
            non_chebi_counter += 1
            continue

        if source_id.startswith(f'{source_db.upper()}:'):
            source_id = source_id[len(source_db) + 1:]
        if target_id.startswith(f'{target_db.upper()}:'):
            target_id = target_id[len(target_db) + 1:]

        if target_db == 'hgnc':
            # Append original
            x[source_db, source_id].append((modulation, 'protein', 'hgnc', target_id, target_name))
            # Append inferred
            for uniprot_id, uniprot_name in get_uniprot_id_names(target_id):
                x[source_db, source_id].append((modulation, 'protein', 'uniprot', uniprot_id, uniprot_name))

        elif target_db == 'fplx':
            # Append original
            x[source_db, source_id].append((modulation, target_type, target_db, target_id, target_name))
            # Append inferred
            for hgnc_id, hgnc_symbol in famplex_id_to_members.get(target_id, []):
                x[source_db, source_id].append((modulation, 'protein', 'hgnc', hgnc_id, hgnc_symbol))
                for uniprot_id, uniprot_name in get_uniprot_id_names(hgnc_id):
                    x[source_db, source_id].append((modulation, 'protein', 'uniprot', uniprot_id, uniprot_name))

        elif target_db == 'eccode':
            children_ec_codes = ec_code_to_children.get(target_id)
            if children_ec_codes is None:
                # this is the case for about 15 entries
                logger.info(f'could not find children of {target_db}:{target_id}')
                continue

            for sub_target_db, sub_target_id, sub_target_name in children_ec_codes:
                target_type = DB_TO_TYPE[sub_target_db]
                x[source_db, source_id].append((
                    modulation, target_type, sub_target_db, sub_target_id, sub_target_name,
                ))

            for go_id, go_name in ec2go.get(target_id, []):
                x[source_db, source_id].append((
                    modulation, 'molecular function', 'go', go_id, go_name,
                ))

        else:
            x[source_db, source_id].append((modulation, target_type, target_db, target_id, target_name))

    logger.info('x mapping: %d/%d', len(x), sum(map(len, x.values())))
    logger.info('skipped %d non-chebi source terms', non_chebi_counter)

    logger.info('inferring over role hiearchies')
    db_to_role_to_chemical_curies = {
        'chebi': get_chebi_role_to_children(),
    }
    for (role_db, role_id), entries in tqdm(sorted(x.items()), desc='inferring over role hierarchies'):
        sub_role_curies = {(role_db, role_id)}

        if role_db == 'chebi' and use_sub_roles:
            sub_role_curies |= {
                pyobo.normalize_curie(c)
                for c in pyobo.get_subhierarchy(role_db, role_id)
            }

        chemical_curies = set(itt.chain.from_iterable(
            db_to_role_to_chemical_curies[sub_role_db].get(sub_role_id, [])
            for sub_role_db, sub_role_id in sub_role_curies
        ))
        if not chemical_curies:
            tqdm.write(f'no inference for {role_db}:{role_id} ! {pyobo.get_name(role_db, role_id)}')
            continue

        for modulation, target_type, target_db, target_id, target_name in entries:
            for chemical_db, chemical_id in chemical_curies:
                rows.append((
                    chemical_db, chemical_id, pyobo.get_name(chemical_db, chemical_id),
                    modulation, target_type, target_db, target_id, target_name,
                ))

    logger.info('inferred df has %d rows', len(rows))
    rv = pd.DataFrame(rows, columns=XREFS_COLUMNS)
    rv.sort_values(XREFS_COLUMNS, inplace=True)
    rv.drop_duplicates(inplace=True)
    return rv


FAMPLEX_RELATIONS_URL = 'https://raw.githubusercontent.com/sorgerlab/famplex/master/relations.csv'
DB_TO_TYPE = {
    'eccode': 'protein family',
    'uniprot': 'protein',
    'prosite': 'protein',
}


def _get_famplex():
    logger.info('loading famplex mapping')
    famplex_id_to_members = defaultdict(list)
    famplex_relations_df = pd.read_csv(FAMPLEX_RELATIONS_URL)
    for source_id, source_name, rel, target_db, target_name in famplex_relations_df.values:
        if source_id.lower() == 'hgnc' and rel == 'isa' and target_db.lower() == 'fplx':
            try:
                hgnc_id = hgnc_name_to_id[source_name]
            except KeyError:
                logger.warning(f'Could not find {source_name} for fplx:{target_name}')
                continue
            famplex_id_to_members[target_name].append((hgnc_id, source_name))

    logger.info('famplex mapping has %d elements', len(famplex_id_to_members))
    return famplex_id_to_members


def get_expasy_closure() -> Tuple[nx.DiGraph, Mapping[str, List[str]]]:
    """Get the ExPASy closure map."""
    _graph = nx.DiGraph()
    expasy_obo = expasy.get_obo()
    for term in expasy_obo:
        for parent_term in term.parents:
            _graph.add_edge(
                (term.prefix, term.identifier, term.identifier),
                (parent_term.prefix, parent_term.identifier, parent_term.identifier),
            )

        for member in term.get_relationships(has_member):
            _graph.add_edge(
                (member.prefix, member.identifier, member.name),
                (term.prefix, term.identifier, term.identifier),
            )

    rv = {
        identifier: list(nx.ancestors(_graph, (db, identifier, name)))
        for db, identifier, name in _graph
        if db == expasy_obo.ontology
    }
    return _graph, rv


def get_uniprot_id_names(hgnc_id: str) -> Iterable[Tuple[str, str]]:
    """Get all of the UniProt identifiers for a given gene."""
    try:
        r = hgnc_id_to_up[str(hgnc_id)]
    except KeyError:
        tqdm.write(f'could not find HGNC:{hgnc_id}')
        return

    for _uniprot_id in r.split(', '):
        yield _uniprot_id, uniprot_client.get_mnemonic(_uniprot_id)

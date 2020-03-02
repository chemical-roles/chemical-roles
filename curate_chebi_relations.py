# -*- coding: utf-8 -*-

"""A script to help curate ChEBI relations and infer new ones."""

import json
import logging
import os
import sys
from collections import defaultdict
from typing import Collection, Iterable, Tuple

import click
import networkx as nx
import pandas as pd
import requests
from protmapper import uniprot_client
from protmapper.api import hgnc_id_to_up
from tabulate import tabulate

from pyobo.sources.expasy import get_obo as get_expasy_obo
from pyobo.utils import get_obo_graph

logger = logging.getLogger(__name__)

HERE = os.path.abspath(os.path.dirname(__file__))
RESOURCES_DIRECTORY = os.path.join(HERE, 'resources')
EXPORT_DIRECTORY = os.path.join(HERE, 'export')

GILDA_URL = 'http://grounding.indra.bio'

FAMPLEX_RELATIONS_URL = 'https://raw.githubusercontent.com/sorgerlab/famplex/master/relations.csv'
FAMPLEX_EQUIVALENCES_URL = 'https://raw.githubusercontent.com/sorgerlab/famplex/master/equivalences.csv'
FAMPLEX_HGNC_SYMBOL_MAP_URL = 'https://raw.githubusercontent.com/sorgerlab/famplex/master/export/hgnc_symbol_map.csv'

RELATIONS_OUTPUT_PATH = os.path.join(EXPORT_DIRECTORY, 'relations.tsv')
RELATIONS_SLIM_OUTPUT_PATH = os.path.join(EXPORT_DIRECTORY, 'relations_slim.tsv')
SUMMARY_OUTPUT_PATH = os.path.join(EXPORT_DIRECTORY, 'relations_summary.tsv')

XREFS_COLUMNS = [
    'source_db', 'source_id', 'source_name',
    'modulation',
    'target_type', 'target_db', 'target_id', 'target_name',
]

XREFS_PATH = os.path.join(RESOURCES_DIRECTORY, 'xrefs.tsv')
RECLASSIFICATION_PATH = os.path.join(RESOURCES_DIRECTORY, 'reclassification.tsv')

BIOCHEMICAL_ROLE_CHEBI_ID = 'CHEBI:52206'
PATHWAY_INHIBITOR_CHEBI_ID = 'CHEBI:76932'
ENZYME_INHIBITOR_CHEBI_ID = 'CHEBI:23924'
AGONIST_CHEBI_ID = 'CHEBI:48705'
INVERSE_AGONIST_CHEBI_ID = 'CHEBI:90847'
INHIBITOR_CHEBI_ID = 'CHEBI:35222'
ANTAGONIST_CHEBI_ID = 'CHEBI:48706'
BLACKLIST = [
    'CHEBI:48001',  # protein synthesis inhibitor
]


def _get_curated_xrefs_df() -> pd.DataFrame:
    return pd.read_csv(XREFS_PATH, sep='\t', comment='#', dtype={'db_id': str})


def _get_inhibitors_reclassification() -> pd.DataFrame:
    return pd.read_csv(RECLASSIFICATION_PATH, sep='\t', comment='#')


def post_gilda(text: str, url: str = GILDA_URL) -> requests.Response:
    """Send text to GILDA."""
    return requests.post(f'{url}/ground', json={'text': text})


def get_expasy_closure():
    """Get the ExPASy closure map."""
    expasy_obo = get_expasy_obo()
    _graph = nx.DiGraph()
    for term in expasy_obo.terms:
        for parent_term in term.parents:
            _graph.add_edge(
                (term.reference.prefix, term.reference.identifier, term.name),
                (parent_term.prefix, parent_term.identifier, parent_term.name),
            )
        for xref in term.xrefs:
            if xref.prefix in {'uniprot', 'prosite'}:
                _graph.add_edge(
                    (xref.prefix, xref.identifier, xref.name),
                    (term.reference.prefix, term.reference.identifier, term.name),
                )

    rv = {
        identifier: list(nx.ancestors(_graph, (db, identifier, name)))
        for db, identifier, name in _graph
        if db == 'ec-code'
    }
    return _graph, rv


def get_enzyme_inhibitor_df(chebi_graph: nx.MultiDiGraph) -> pd.DataFrame:
    expasy_graph, ec_code_to_children = get_expasy_closure()
    rv = []

    source_db = 'chebi'
    for chebi_id, data in chebi_graph.nodes(data=True):
        chebi_name = data['name']

        # Do this as a loop since there is at least one entry that corresponds to several EC codes
        ec_codes = []

        # Requested fix in https://github.com/ebi-chebi/ChEBI/issues/3651
        if chebi_name == 'EC 1.22* (oxidoreductase acting on halogen in donors) inhibitor':
            modulation = 'inhibitor'
            ec_codes.append('1.22.-.-')

        # Not sure why this has two
        elif chebi_name == 'EC 1.1.1.34/EC 1.1.1.88 (hydroxymethylglutaryl-CoA reductase) inhibitor':
            modulation = 'inhibitor'
            ec_codes.append('1.1.1.34')
            ec_codes.append('1.1.1.88')

        # Requested rename in https://github.com/ebi-chebi/ChEBI/issues/3653
        elif chebi_name == 'EC 1.11.1.11 (L-ascorbate peroxidase) inhibitors':
            modulation = 'inhibitor'
            ec_codes.append('1.11.1.11')

        # Requested typo fix in https://github.com/ebi-chebi/ChEBI/issues/3652
        elif chebi_name == 'EC 3.5.5.1 (nitrilase) inhhibitor':
            modulation = 'inhibitor'
            ec_codes.append('3.5.5.1')

        # All other normal cases
        elif chebi_name.startswith('EC '):
            if chebi_name.endswith('inhibitor'):
                modulation = 'inhibitor'
            elif chebi_name.endswith('activator'):
                modulation = 'activator'
            else:
                logger.warning(f'Unhandled suffix: chebi:{chebi_id} ! {chebi_name}')
                continue

            ec_code = chebi_name[len('EC '):].split()[0].replace('*', '-').rstrip('.')

            # Add any remaining dashes
            for _ in range(3 - ec_code.count('.')):
                ec_code += '.-'

            ec_codes.append(ec_code)

        else:
            continue

        for ec_code in ec_codes:
            rv.append((
                source_db, chebi_id, chebi_name, modulation, 'protein family', 'ec-code', ec_code, ec_code,
            ))

            children_ec_codes = ec_code_to_children.get(ec_code)
            if children_ec_codes is None:
                print(f'could not find {ec_code} (for chebi:{chebi_id} ! {chebi_name})')
                continue

            for target_db, target_id, target_name in children_ec_codes:
                target_type = DB_TO_TYPE[target_db]
                rv.append((
                    source_db, chebi_id, chebi_name, modulation, target_type, target_db, target_id, target_name,
                ))

    return pd.DataFrame(rv, columns=XREFS_COLUMNS)


DB_TO_TYPE = {'ec-code': 'protein family', 'uniprot': 'protein', 'prosite': 'protein'}


def suggest_pathway_inhibitor_curation(graph: nx.MultiDiGraph) -> None:
    inhibitors = _get_curated_xrefs_df()
    curated_chebi_ids = set(inhibitors.chebi_id)

    reclassify_df = _get_inhibitors_reclassification()
    reclassify_chebi_ids = set(reclassify_df.chebi_id)

    print(f'Children of {PATHWAY_INHIBITOR_CHEBI_ID} ({graph.nodes[PATHWAY_INHIBITOR_CHEBI_ID]["name"]})')
    for node in nx.ancestors(graph, PATHWAY_INHIBITOR_CHEBI_ID):
        if any(node in group for group in (curated_chebi_ids, reclassify_chebi_ids, BLACKLIST)):
            continue  # we already curated this!
        name = graph.nodes[node]['name']
        if name.endswith('inhibitor') and not name.startswith('EC '):
            print(node, name)
            results = post_gilda(name[:-len(' inhibitor')]).json()
            if results:
                print(json.dumps(results, indent=2))


def suggest_inhibitor_curation(graph: nx.MultiDiGraph) -> None:
    it = (
        set(nx.ancestors(graph, INHIBITOR_CHEBI_ID))
        - set(nx.ancestors(graph, PATHWAY_INHIBITOR_CHEBI_ID))
        - set(nx.ancestors(graph, ENZYME_INHIBITOR_CHEBI_ID))
    )
    for t in _suggest_xrefs_curation(graph, INHIBITOR_CHEBI_ID, it, 'inhibitor'):
        print(*t, sep='\t')


def suggest_agonist_curation(graph: nx.MultiDiGraph) -> None:
    _single_suggest(graph, AGONIST_CHEBI_ID, 'agonist')


def suggest_antagonist_curation(graph: nx.MultiDiGraph) -> None:
    _single_suggest(graph, ANTAGONIST_CHEBI_ID, 'antagonist')


def suggest_inverse_agonist_curation(graph: nx.MultiDiGraph) -> None:
    _single_suggest(graph, INVERSE_AGONIST_CHEBI_ID, 'inverse agonist')


def suggest_activator_curation(graph: nx.MultiDiGraph) -> None:
    _single_suggest(graph, BIOCHEMICAL_ROLE_CHEBI_ID, 'activator')


def _single_suggest(graph, chebi_id, modulation, file=None) -> None:
    it = set(nx.ancestors(graph, chebi_id))
    for t in _suggest_xrefs_curation(graph, chebi_id, it, modulation):
        print(*t, sep='\t', file=file)


def _suggest_xrefs_curation(
    graph: nx.MultiDiGraph,
    top_chebi_id: str,
    it: Collection[str],
    suffix: str,
    show_missing: bool = False
) -> Iterable[Tuple[str, str, str, str, str, str]]:
    """

    :param graph:
    :param top_chebi_id:
    :param it:
    :param suffix: The modulation type (e.g., inhibitor, agonist)
    """
    xrefs_df = _get_curated_xrefs_df()
    curated_chebi_ids = set(xrefs_df.chebi_id)

    if it:
        print(f'\nChildren of {top_chebi_id} ({graph.nodes[top_chebi_id]["name"]})\n')
    for node in it:
        if node in curated_chebi_ids:
            continue
        name = graph.nodes[node]['name']
        if not name.endswith(suffix):
            continue
        results = post_gilda(name[:-len(suffix)].rstrip()).json()
        if results:
            for result in results:
                term = result["term"]
                yield node, name, suffix, term['db'].lower(), term['id'], term['entry_name']
        elif show_missing:
            yield node, name, suffix, '?', '?', '?', '?'


def get_relations_df(graph: nx.MultiDiGraph) -> pd.DataFrame:
    xrefs_df = _get_curated_xrefs_df()

    famplex_map_df = pd.read_csv(FAMPLEX_HGNC_SYMBOL_MAP_URL)
    hgnc_symbol_to_id = dict(famplex_map_df.values)

    famplex_id_to_members = defaultdict(list)
    famplex_relations_df = pd.read_csv(FAMPLEX_RELATIONS_URL)
    for source_id, source_name, rel, target_db, target_name in famplex_relations_df.values:
        if source_id == 'HGNC' and rel == 'isa' and target_db == 'FPLX':
            try:
                hgnc_symbol = hgnc_symbol_to_id[source_name]
            except KeyError:
                logger.warning(f'Could not find {source_name} for fplx:{target_name}')
                continue

            famplex_id_to_members[target_name].append((hgnc_symbol, source_name))

    def _get_uniprot_id_names(hgnc_id: str) -> Iterable[Tuple[str, str]]:
        try:
            r = hgnc_id_to_up[str(hgnc_id)]
        except KeyError:
            _k, _v = list(hgnc_id_to_up.items())[0]
            print(f'could not find {hgnc_id} ({type(hgnc_id)} in dict. Example: {_k} ({type(_k)}), {_v} ({type(_v)})')
            raise

        for _uniprot_id in r.split(', '):
            yield _uniprot_id, uniprot_client.get_mnemonic(_uniprot_id)

    xrefs = defaultdict(list)
    for source_db, source_id, _, modulation, target_type, target_db, target_id, target_name in xrefs_df.values:
        if target_db == 'hgnc':
            for uniprot_id, uniprot_name in _get_uniprot_id_names(target_id):
                xrefs[source_id].append((modulation, 'protein', 'uniprot', uniprot_id, uniprot_name))

        elif target_db == 'fplx':
            xrefs[source_id].append((modulation, target_type, target_db, target_id, target_name))

            for hgnc_id, hgnc_symbol in famplex_id_to_members.get(target_id, []):
                xrefs[source_id].append((modulation, 'protein', 'hgnc', hgnc_id, hgnc_symbol))
                for uniprot_id, uniprot_name in _get_uniprot_id_names(hgnc_id):
                    xrefs[source_id].append((modulation, 'protein', 'uniprot', uniprot_id, uniprot_name))
        else:
            xrefs[source_id].append((modulation, target_type, target_db, target_id, target_name))

    rv = [
        ('chebi', child_chebi_id, child_name, *xref)
        for child_chebi_id, child_name, role_chebi_id in _iterate_roles(graph, xrefs_df['source_id'])
        for xref in xrefs.get(role_chebi_id, [])
    ]
    return pd.DataFrame(rv, columns=XREFS_COLUMNS)


def _iterate_roles(graph: nx.MultiDiGraph, chebi_ids):
    for chebi_id in chebi_ids:
        for child_chebi_id, _ in graph.in_edges(chebi_id):
            child_data = graph.nodes[child_chebi_id]
            child_name = child_data['name']
            relationships = defaultdict(list)
            for r in child_data.get('relationship', []):
                role, role_chebi_id = r.split()
                relationships[role].append(role_chebi_id)
            relationships = dict(relationships)
            for role_chebi_id in relationships.get('has_role', []):
                yield child_chebi_id, child_name, role_chebi_id


@click.command()
@click.option('-s', '--suggest', is_flag=True)
@click.option('-b', '--debug', is_flag=True)
def main(suggest: bool, debug: bool) -> None:
    level = logging.DEBUG if debug else logging.INFO
    logger.setLevel(level)
    logging.basicConfig(level=level)

    # Sort the curated xrefs file
    (
        _get_curated_xrefs_df()
            .sort_values(['source_db', 'source_name', 'modulation'])
            .drop_duplicates()
            .to_csv(XREFS_PATH, index=False, sep='\t')
    )

    # Get the graph
    graph = get_obo_graph('chebi')

    if suggest:
        suggest_activator_curation(graph)
        suggest_pathway_inhibitor_curation(graph)
        suggest_inhibitor_curation(graph)
        suggest_agonist_curation(graph)
        suggest_antagonist_curation(graph)
        suggest_inverse_agonist_curation(graph)
        return sys.exit(0)

    relations_df = get_relations_df(graph)
    enzyme_inhibitor_df = get_enzyme_inhibitor_df(graph)
    df = pd.concat([
        relations_df,
        enzyme_inhibitor_df,
    ]).drop_duplicates()

    columns = [
        'modulation', 'target_type', 'source_db', 'source_id', 'source_name', 'target_db', 'target_id', 'target_name',
    ]
    df[columns].sort_values(columns).to_csv(RELATIONS_OUTPUT_PATH, sep='\t', index=False)

    slim_columns = ['source_db', 'source_id', 'modulation', 'target_type', 'target_db', 'target_id']
    df[slim_columns].sort_values(slim_columns).to_csv(RELATIONS_SLIM_OUTPUT_PATH, sep='\t', index=False)

    summary_df = df.groupby(['source_db', 'modulation', 'target_type', 'target_db']).size().reset_index()
    summary_df.to_csv(SUMMARY_OUTPUT_PATH, sep='\t', index=False)

    print(tabulate(
        summary_df.values,
        ['source_db', 'relation', 'target_type', 'target_db', 'count'],
        tablefmt='github',
    ))


if __name__ == '__main__':
    main()

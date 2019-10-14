"""A script to help curate ChEBI relations and infer new ones."""

import json
import logging
import os
import pickle
from collections import defaultdict
from typing import Collection, Iterable, Tuple
from urllib.request import urlretrieve

import click
import obonet
import pandas as pd
import requests
from bio2bel_expasy.parser import get_expasy_closed_tree
from networkx import MultiDiGraph, ancestors

logger = logging.getLogger(__name__)

HERE = os.path.abspath(os.path.dirname(__file__))
RESOURCES_DIRECTORY = os.path.join(HERE, 'resources')
EXPORT_DIRECTORY = os.path.join(HERE, 'export')

GILDA_URL = 'http://grounding.indra.bio'

CHEBI_RELEASE = '179'
CHEBI_OBO_URL = f'ftp://ftp.ebi.ac.uk/pub/databases/chebi/archive/rel{CHEBI_RELEASE}/ontology/chebi.obo.gz'
CHEBI_OBO_PATH = os.path.join(RESOURCES_DIRECTORY, 'chebi.obo.gz')
CHEBI_OBO_PICKLE_PATH = os.path.join(RESOURCES_DIRECTORY, 'chebi.obo.pickle')

FAMPLEX_RELATIONS_URL = 'https://raw.githubusercontent.com/sorgerlab/famplex/master/relations.csv'
FAMPLEX_EQUIVALENCES_URL = 'https://raw.githubusercontent.com/sorgerlab/famplex/master/equivalences.csv'
FAMPLEX_HGNC_SYMBOL_MAP_URL = 'https://raw.githubusercontent.com/sorgerlab/famplex/master/export/hgnc_symbol_map.csv'

RELATIONS_OUTPUT_PATH = os.path.join(EXPORT_DIRECTORY, 'relations.tsv')
RELATIONS_SLIM_OUTPUT_PATH = os.path.join(EXPORT_DIRECTORY, 'relations_slim.tsv')
SUMMARY_OUTPUT_PATH = os.path.join(EXPORT_DIRECTORY, 'relations_summary.tsv')

XREFS_COLUMNS = ['chebi_id', 'chebi_name', 'modulation', 'entity_type', 'db', 'db_id', 'db_name']

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
    return pd.read_csv(XREFS_PATH, sep='\t', comment='#')


def _get_inhibitors_reclassification() -> pd.DataFrame:
    return pd.read_csv(RECLASSIFICATION_PATH, sep='\t', comment='#')


def post_gilda(text: str, url: str = GILDA_URL) -> requests.Response:
    """Send text to GILDA."""
    return requests.post(f'{url}/ground', json={'text': text})


def get_graph(
    path: str = CHEBI_OBO_PATH,
    pickle_path: str = CHEBI_OBO_PICKLE_PATH,
) -> MultiDiGraph:
    if os.path.exists(pickle_path):
        logger.info('Opening cached pre-parsed OBO')
        with open(pickle_path, 'rb') as file:
            return pickle.load(file)
    if not os.path.exists(path):
        logger.info(f'Downloading from {CHEBI_OBO_URL}')
        urlretrieve(CHEBI_OBO_URL, path)

    logger.info('Parsing OBO')
    graph = obonet.read_obo(path)

    with open(pickle_path, 'wb') as file:
        logger.info('Caching pre-parsed OBO')
        pickle.dump(graph, file)

    return graph


def get_enzyme_inhibitor_df(graph: MultiDiGraph) -> pd.DataFrame:
    expasy = get_expasy_closed_tree()

    rv = []
    for chebi_id, data in graph.nodes(data=True):
        chebi_name = data['name']

        # Do this as a loop since there is at least one entry that corresponds to several EC codes
        ec_codes = []

        # Requested fix in https://github.com/ebi-chebi/ChEBI/issues/3651
        if chebi_name == 'EC 1.22* (oxidoreductase acting on halogen in donors) inhibitor':
            ec_codes.append('1.22.-.-')

        # Not sure why this has two
        elif chebi_name == 'EC 1.1.1.34/EC 1.1.1.88 (hydroxymethylglutaryl-CoA reductase) inhibitor':
            ec_codes.append('1.1.1.34')
            ec_codes.append('1.1.1.88')

        # Requested rename in https://github.com/ebi-chebi/ChEBI/issues/3653
        elif chebi_name == 'EC 1.11.1.11 (L-ascorbate peroxidase) inhibitors':
            ec_codes.append('1.11.1.11')

        # Requested typo fix in https://github.com/ebi-chebi/ChEBI/issues/3652
        elif chebi_name == 'EC 3.5.5.1 (nitrilase) inhhibitor':
            ec_codes.append('3.5.5.1')

        # All other normal cases
        elif chebi_name.startswith('EC '):
            if chebi_name.endswith('inhibitor'):
                modulation = 'inhibitor'
            elif chebi_name.endswith('activator'):
                modulation = 'activator'
            else:
                logger.warning(f'Unhandled suffix: {chebi_id} ! {chebi_name}')
                continue

            ec_code = chebi_name[len('EC '):].split()[0].replace('*', '-').rstrip('.')

            # Add any remaining dashes
            for _ in range(3 - ec_code.count('.')):
                ec_code += '.-'

            ec_codes.append(ec_code)

        else:
            continue

        for ec_code in ec_codes:
            rv.append((chebi_id, chebi_name, modulation, 'enzyme', 'ec-code', ec_code, ec_code))

            expasy_children = expasy.get(ec_code)
            if expasy_children is None:
                print(f'could not find {ec_code} (for chebi:{chebi_id} ! {chebi_name})')
                continue

            for c_db, c_identifier, c_name in expasy_children:
                if c_db == 'ec-code':
                    entity_type = 'enzyme'
                else:
                    entity_type = 'protein'

                rv.append((chebi_id, chebi_name, modulation, entity_type, c_db, c_identifier, c_name or c_identifier))

    return pd.DataFrame(rv, columns=XREFS_COLUMNS)


def suggest_pathway_inhibitor_curation(graph: MultiDiGraph) -> None:
    inhibitors = _get_curated_xrefs_df()
    curated_chebi_ids = set(inhibitors.chebi_id)

    reclassify_df = _get_inhibitors_reclassification()
    reclassify_chebi_ids = set(reclassify_df.chebi_id)

    print(f'Children of {PATHWAY_INHIBITOR_CHEBI_ID} ({graph.nodes[PATHWAY_INHIBITOR_CHEBI_ID]["name"]})')
    for node in ancestors(graph, PATHWAY_INHIBITOR_CHEBI_ID):
        if any(node in group for group in (curated_chebi_ids, reclassify_chebi_ids, BLACKLIST)):
            continue  # we already curated this!
        name = graph.nodes[node]['name']
        if name.endswith('inhibitor') and not name.startswith('EC '):
            print(node, name)
            results = post_gilda(name[:-len(' inhibitor')]).json()
            if results:
                print(json.dumps(results, indent=2))


def suggest_inhibitor_curation(graph: MultiDiGraph) -> None:
    it = (
        set(ancestors(graph, INHIBITOR_CHEBI_ID))
        - set(ancestors(graph, PATHWAY_INHIBITOR_CHEBI_ID))
        - set(ancestors(graph, ENZYME_INHIBITOR_CHEBI_ID))
    )
    for t in _suggest_xrefs_curation(graph, INHIBITOR_CHEBI_ID, it, 'inhibitor'):
        print(*t, sep='\t')


def suggest_agonist_curation(graph: MultiDiGraph) -> None:
    _single_suggest(graph, AGONIST_CHEBI_ID, 'agonist')


def suggest_antagonist_curation(graph: MultiDiGraph) -> None:
    _single_suggest(graph, ANTAGONIST_CHEBI_ID, 'antagonist')


def suggest_inverse_agonist_curation(graph: MultiDiGraph) -> None:
    _single_suggest(graph, INVERSE_AGONIST_CHEBI_ID, 'inverse agonist')


def suggest_activator_curation(graph: MultiDiGraph) -> None:
    _single_suggest(graph, BIOCHEMICAL_ROLE_CHEBI_ID, 'activator')


def _single_suggest(graph, chebi_id, modulation, file=None) -> None:
    it = set(ancestors(graph, chebi_id))
    for t in _suggest_xrefs_curation(graph, chebi_id, it, modulation):
        print(*t, sep='\t', file=file)


def _suggest_xrefs_curation(
    graph: MultiDiGraph,
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
    xrefs = _get_curated_xrefs_df()
    curated_chebi_ids = set(xrefs.chebi_id)

    print(f'Children of {top_chebi_id} ({graph.nodes[top_chebi_id]["name"]})')
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


def get_relations_df(graph: MultiDiGraph) -> pd.DataFrame:
    xrefs_df = _get_curated_xrefs_df()

    famplex_map_df = pd.read_csv(FAMPLEX_HGNC_SYMBOL_MAP_URL)
    hgnc_symbol_to_id = dict(famplex_map_df.values)

    famplex_id_to_members = defaultdict(list)
    famplex_relations_df = pd.read_csv(FAMPLEX_RELATIONS_URL)
    for source_db, source_name, rel, target_db, target_name in famplex_relations_df.values:
        if source_db == 'HGNC' and rel == 'isa' and target_db == 'FPLX':
            try:
                hgnc_symbol = hgnc_symbol_to_id[source_name]
            except KeyError:
                logger.warning(f'Could not find {source_name} for fplx:{target_name}')
                continue

            famplex_id_to_members[target_name].append((hgnc_symbol, source_name))

    xrefs = defaultdict(list)
    for chebi_id, _, modulation, entity_type, db, db_id, db_name in xrefs_df.values:
        xrefs[chebi_id].append((modulation, entity_type, db, db_id, db_name))
        if db == 'fplx':
            for hgnc_id, hgnc_symbol in famplex_id_to_members.get(db_id, []):
                xrefs[chebi_id].append((modulation, 'protein', 'hgnc', hgnc_id, hgnc_symbol))

    rv = [
        (child_chebi_id, child_name, *xref)
        for child_chebi_id, child_name, role_chebi_id in _iterate_roles(graph, xrefs_df.chebi_id)
        for xref in xrefs.get(role_chebi_id, [])
    ]
    return pd.DataFrame(rv, columns=XREFS_COLUMNS)


def _iterate_roles(graph, chebi_ids):
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
    _get_curated_xrefs_df().sort_values(['chebi_name', 'modulation']).to_csv(XREFS_PATH, index=False, sep='\t')

    # Get the graph
    graph = get_graph()

    if suggest:
        suggest_activator_curation(graph)
        suggest_pathway_inhibitor_curation(graph)
        suggest_inhibitor_curation(graph)
        suggest_agonist_curation(graph)
        suggest_antagonist_curation(graph)
        suggest_inverse_agonist_curation(graph)

    relations_df = get_relations_df(graph)
    enzyme_inhibitor_df = get_enzyme_inhibitor_df(graph)
    df = pd.concat([
        relations_df,
        enzyme_inhibitor_df,
    ]).drop_duplicates()

    columns = ['modulation', 'entity_type', 'chebi_id', 'chebi_name', 'db', 'db_id', 'db_name']
    df[columns].sort_values(columns).to_csv(RELATIONS_OUTPUT_PATH, sep='\t', index=False)

    slim_columns = ['chebi_id', 'modulation', 'entity_type', 'db', 'db_id']
    df[slim_columns].sort_values(slim_columns).to_csv(RELATIONS_SLIM_OUTPUT_PATH, sep='\t', index=False)

    summary_df = df.groupby(['modulation', 'entity_type', 'db']).size().reset_index()
    summary_df.to_csv(SUMMARY_OUTPUT_PATH, sep='\t', index=False)

    try:
        from tabulate import tabulate
    except ImportError:
        pass
    else:
        print(tabulate(summary_df.values, ['relation', 'type', 'db', 'count'], tablefmt='github'))


if __name__ == '__main__':
    main()

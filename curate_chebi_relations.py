import json
import os
import pickle
from collections import defaultdict
from typing import Collection, Iterable, Tuple
from urllib.request import urlretrieve

import obonet
import pandas as pd
import requests
from networkx import MultiDiGraph, ancestors

HERE = os.path.abspath(os.path.dirname(__file__))
RESOURCES_DIRECTORY = os.path.join(HERE, 'resources')
EXPORT_DIRECTORY = os.path.join(HERE, 'export')

GILDA_URL = 'http://34.201.164.108:8001'

CHEBI_OBO_URL = 'ftp://ftp.ebi.ac.uk/pub/databases/chebi/ontology/chebi.obo'
CHEBI_OBO_PATH = os.path.join(RESOURCES_DIRECTORY, 'chebi.obo')
CHEBI_OBO_PICKLE_PATH = os.path.join(RESOURCES_DIRECTORY, 'chebi.obo.pickle')

RELATIONS_OUTPUT_PATH = os.path.join(EXPORT_DIRECTORY, 'relations.tsv')
XREFS_COLUMNS = ['chebi_id', 'chebi_name', 'modulation', 'entity_type', 'db', 'db_id', 'db_name']

XREFS_PATH = os.path.join(RESOURCES_DIRECTORY, 'xrefs.tsv')
RECLASSIFICATION_PATH = os.path.join(RESOURCES_DIRECTORY, 'reclassification.tsv')

PATHWAY_INHIBITOR_CHEBI_ID = 'CHEBI:76932'
ENZYME_INHIBITOR_CHEBI_ID = 'CHEBI:23924'
AGONIST_CHEBI_ID = 'CHEBI:48705'
INVERSE_AGONIST_CHEBI_ID = 'CHEBI:90847'
INHIBITOR_CHEBI_ID = 'CHEBI:35222'
ANTAGONIST_CHEBI_ID = 'CHEBI:48706'
BLACKLIST = [
    'CHEBI:48001',  # protein synthesis inhibitor
]


def _get_curated_xrefs_df():
    return pd.read_csv(XREFS_PATH, sep='\t', comment='#')


def _get_inhibitors_reclassification():
    return pd.read_csv(RECLASSIFICATION_PATH, sep='\t', comment='#')


def post_gilda(text: str, url: str = GILDA_URL) -> requests.Response:
    """Send text to GILDA."""
    return requests.post(f'{url}/ground', json={'text': text})


def get_graph(
        path: str = CHEBI_OBO_PATH,
        pickle_path: str = CHEBI_OBO_PICKLE_PATH,
) -> MultiDiGraph:
    if os.path.exists(pickle_path):
        with open(pickle_path, 'rb') as file:
            return pickle.load(file)
    if not os.path.exists(path):
        urlretrieve(CHEBI_OBO_URL, path)
    graph = obonet.read_obo(path)
    with open(pickle_path, 'wb') as file:
        pickle.dump(graph, file)
    return graph


def get_enzyme_inhibitor_df(graph: MultiDiGraph) -> pd.DataFrame:
    rv = []
    for chebi_id, data in graph.nodes(data=True):
        chebi_name = data['name']
        if chebi_name.startswith('EC ') and chebi_name.endswith('inhibitor'):
            ec_code = chebi_name[len('EC '):].split()[0]
            if ec_code.endswith('.*'):
                ec_code = ec_code[:-len('.*')]
            rv.append((chebi_id, chebi_name, 'inhibitor', 'enzyme', 'ec-code', ec_code, ec_code))

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


def _single_suggest(graph, chebi_id, modulation, file=None) -> None:
    it = set(ancestors(graph, chebi_id))
    for t in _suggest_xrefs_curation(graph, chebi_id, it, modulation):
        print(*t, sep='\t', file=file)


def _suggest_xrefs_curation(
        graph: MultiDiGraph,
        top_chebi_id: str,
        it: Collection[str],
        suffix: str,
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
        if not results:
            yield node, name, suffix, '?', '?', '?', '?'
        else:
            for result in results:
                term = result["term"]
                yield node, name, suffix, term['db'].lower(), term['id'], term['entry_name']


def get_relations_df(graph: MultiDiGraph) -> pd.DataFrame:
    xrefs_df = _get_curated_xrefs_df()
    xrefs = {
        chebi_id: (modulation, entity_type, db, db_id, db_name)
        for chebi_id, _, modulation, entity_type, db, db_id, db_name in xrefs_df.values
    }

    rv = [
        (child_chebi_id, child_name, *xrefs[role_chebi_id])
        for child_chebi_id, child_name, role_chebi_id in _iterate_roles(graph, xrefs_df.chebi_id)
        if role_chebi_id in xrefs
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


def main() -> None:
    graph = get_graph()

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
    ]).sort_values(['modulation', 'entity_type', 'chebi_id'])
    df.to_csv(RELATIONS_OUTPUT_PATH, sep='\t', index=False)


if __name__ == '__main__':
    main()

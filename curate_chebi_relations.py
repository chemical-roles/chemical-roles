import json
import os
import pickle
from collections import defaultdict
from urllib.request import urlretrieve

import obonet
import pandas as pd
import requests
from networkx import ancestors

HERE = os.path.abspath(os.path.dirname(__file__))
RESOURCES_DIRECTORY = os.path.join(HERE, 'resources')
EXPORT_DIRECTORY =  os.path.join(HERE, 'export')

GILDA_URL = 'http://34.201.164.108:8001'

CHEBI_OBO_URL = 'ftp://ftp.ebi.ac.uk/pub/databases/chebi/ontology/chebi.obo'
CHEBI_OBO_PATH = os.path.join(RESOURCES_DIRECTORY, 'chebi.obo')
CHEBI_OBO_PICKLE_PATH = os.path.join(RESOURCES_DIRECTORY, 'chebi.obo.pickle')


PATHWAY_OUTPUT_PATH = os.path.join(EXPORT_DIRECTORY, 'pathways.tsv')
ENZYME_OUTPUT_PATH = os.path.join(EXPORT_DIRECTORY, 'enzymes.tsv')
INHIBITORS_PATH = os.path.join(RESOURCES_DIRECTORY, 'inhibitors.tsv')
INHIBITORS_RECLASSIFICATION_PATH = os.path.join(RESOURCES_DIRECTORY, 'inhibitors_reclassification.tsv')

PATHWAY_INHIBITOR_CHEBI_ID = 'CHEBI:76932'
ENZYME_INHIBITOR_CHEBI_ID = 'CHEBI:23924'
INHIBITOR_CHEBI_ID = 'CHEBI:35222'
BLACKLIST = [
    'CHEBI:48001',  # protein synthesis inhibitor
]

COLUMNS = ['chebi_id', 'chebi_name', 'mode', 'db', 'db_id', 'db_name']


def _get_inhibitors():
    return pd.read_csv(INHIBITORS_PATH, sep='\t', comment='#')


def _get_pathway_inhibitors_reclassification():
    return pd.read_csv(INHIBITORS_RECLASSIFICATION_PATH, sep='\t', comment='#')


def post_gilda(text: str, url: str = GILDA_URL) -> requests.Response:
    """Send text to GILDA."""
    return requests.post(f'{url}/ground', json={'text': text})


def get_graph(path: str = CHEBI_OBO_PATH, pickle_path: str = CHEBI_OBO_PICKLE_PATH):
    if os.path.exists(pickle_path):
        with open(pickle_path, 'rb') as file:
            return pickle.load(file)
    if not os.path.exists(path):
        urlretrieve(CHEBI_OBO_URL, path)
    graph = obonet.read_obo(path)
    with open(pickle_path, 'wb') as file:
        pickle.dump(graph, file)
    return graph


def get_enzyme_inhibitor_df(graph):
    rv = []
    for chebi_id, data in graph.nodes(data=True):
        chebi_name = data['name']
        if chebi_name.startswith('EC ') and chebi_name.endswith('inhibitor'):
            ec_code = chebi_name[len('EC '):].split()[0]
            if ec_code.endswith('.*'):
                ec_code = ec_code[:-len('.*')]
            rv.append((chebi_id, chebi_name, 'enzyme', 'ec-code', ec_code, ec_code))

    return pd.DataFrame(rv, columns=COLUMNS).sort_values('chebi_id')


def suggest_pathway_inhibitor_curation(graph):
    inhibitors = _get_inhibitors()
    curated_chebi_ids = set(inhibitors.chebi_id)

    reclassify_df = _get_pathway_inhibitors_reclassification()
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


def suggest_inhibitor_curation(graph):
    inhibitors = _get_inhibitors()
    curated_chebi_ids = set(inhibitors.chebi_id)

    print(f'Children of {INHIBITOR_CHEBI_ID} ({graph.nodes[INHIBITOR_CHEBI_ID]["name"]})')
    it = (
            set(ancestors(graph, INHIBITOR_CHEBI_ID))
            - set(ancestors(graph, PATHWAY_INHIBITOR_CHEBI_ID))
            - set(ancestors(graph, ENZYME_INHIBITOR_CHEBI_ID))
    )
    for node in it:
        if node in curated_chebi_ids:
            continue
        name = graph.nodes[node]['name']
        if name.endswith('inhibitor'):
            results = post_gilda(name[:-len(' inhibitor')]).json()
            if results:
                for result in results:
                    term = result["term"]
                    print(node, name, term['db'].lower(), term['id'], term['entry_name'], sep='\t')
            else:
                print(node, name, '?', '?', '?', '?', sep='\t')


def get_inhibitor_df(graph):
    pathway_inhibitors_df = _get_inhibitors()

    d = {
        chebi_id: (mode, db, db_id, db_name)
        for mode, chebi_id, _, db, db_id, db_name in pathway_inhibitors_df.values
    }

    rv = []
    for chebi_id in pathway_inhibitors_df.chebi_id:
        for child_chebi_id, _ in graph.in_edges(chebi_id):
            child_data = graph.nodes[child_chebi_id]
            child_name = child_data['name']
            relationships = defaultdict(list)
            for r in child_data.get('relationship', []):
                role, role_chebi_id = r.split()
                relationships[role].append(role_chebi_id)
            relationships = dict(relationships)
            for role_chebi_id in relationships.get('has_role', []):
                incident = d.get(role_chebi_id, )
                if incident is not None:
                    rv.append((child_chebi_id, child_name, *incident))

    return pd.DataFrame(rv, columns=COLUMNS).sort_values(['mode', 'chebi_id'])


def main():
    graph = get_graph()

    suggest_pathway_inhibitor_curation(graph)
    suggest_inhibitor_curation(graph)

    enzyme_inhibitor_df = get_enzyme_inhibitor_df(graph)
    enzyme_inhibitor_df.to_csv(ENZYME_OUTPUT_PATH, sep='\t', index=False)

    inhibitor_df = get_inhibitor_df(graph)
    inhibitor_df.to_csv(PATHWAY_OUTPUT_PATH, sep='\t', index=False)


if __name__ == '__main__':
    main()

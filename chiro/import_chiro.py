import logging

from pyobo import get_id_name_mapping, get_obo_graph

RELATIONSHIPS = [
    'activator_of',
    'agonist_of',
    'antagonist_of',
    'destabilizer_of',
    'inducer_of',
    'inhibitor_of',
    'modulator_of',
    'sensitizer_of',
    'stabilizier_of',
]

MAPPING_PREFIXES = ['ncbitaxon', 'go', 'pr', 'hp', 'mp']


def main():
    graph = get_obo_graph('chiro')
    chebi_mapping = get_id_name_mapping('chebi')
    mappings = {
        prefix: get_id_name_mapping(prefix)
        for prefix in MAPPING_PREFIXES
    }

    triples = []
    for h, data in graph.nodes(data=True):
        if not data:
            continue
        r, t = data['relationship'][0].split()
        r = r[:-len('_of')]

        h_name = chebi_mapping.get(h)
        if h_name is None:
            print(f'Could not find name for chemical {h}')
            continue

        t_namespace = t.split(':')[0].lower()
        t_mapping = mappings[t_namespace]
        t_name = t_mapping.get(t)
        if t_name is None:
            print(f'Could not find name for target {t}')
            continue

        triples.append(('chebi', h, h_name, r, t_namespace, t, t_name))

    with open('chiro_import.tsv', 'w') as file:
        print('source_db	source_id	source_name	modulation	type	target_db	target_id	target_name',
              file=file)
        for t in sorted(triples):
            print(*t, sep='\t', file=file)


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    main()

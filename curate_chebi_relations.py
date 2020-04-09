# -*- coding: utf-8 -*-

"""A script to help curate ChEBI relations and infer new ones."""

import json
import logging
import os
import sys
from collections import defaultdict
from typing import Iterable, List, Mapping, Tuple

import click
import networkx as nx
import pandas as pd
import pyobo
from protmapper import uniprot_client
from protmapper.api import hgnc_id_to_up
from pyobo import Obo
from pyobo.sources.expasy import get_obo as get_expasy_obo
from tabulate import tabulate

from utils import (
    EXPORT_DIRECTORY, RESOURCES_DIRECTORY, XREFS_COLUMNS, get_xrefs_df, post_gilda,
    sort_xrefs_df,
)

logger = logging.getLogger(__name__)

FAMPLEX_RELATIONS_URL = 'https://raw.githubusercontent.com/sorgerlab/famplex/master/relations.csv'
FAMPLEX_EQUIVALENCES_URL = 'https://raw.githubusercontent.com/sorgerlab/famplex/master/equivalences.csv'
FAMPLEX_HGNC_SYMBOL_MAP_URL = 'https://raw.githubusercontent.com/sorgerlab/famplex/master/export/hgnc_symbol_map.csv'

RELATIONS_OUTPUT_PATH = os.path.join(EXPORT_DIRECTORY, 'relations.tsv')
RELATIONS_SLIM_OUTPUT_PATH = os.path.join(EXPORT_DIRECTORY, 'relations_slim.tsv')
SUMMARY_OUTPUT_PATH = os.path.join(EXPORT_DIRECTORY, 'relations_summary.tsv')

RECLASSIFICATION_PATH = os.path.join(RESOURCES_DIRECTORY, 'reclassification.tsv')

BIOCHEMICAL_ROLE_CHEBI_ID = '52206'
PATHWAY_INHIBITOR_CHEBI_ID = '76932'
ENZYME_INHIBITOR_CHEBI_ID = '23924'
AGONIST_CHEBI_ID = '48705'
INVERSE_AGONIST_CHEBI_ID = '90847'
INHIBITOR_CHEBI_ID = '35222'
ANTAGONIST_CHEBI_ID = '48706'
BLACKLIST = [
    '48001',  # protein synthesis inhibitor
    '64106',  # protein kinase agonist
]

expasy_obo = get_expasy_obo()
chebi_obo = pyobo.get('chebi')
chebi_id_to_name = pyobo.get_id_name_mapping('chebi')


def _get_inhibitors_reclassification() -> pd.DataFrame:
    return pd.read_csv(RECLASSIFICATION_PATH, sep='\t', comment='#')


def get_expasy_closure() -> Tuple[nx.DiGraph, Mapping[str, List[str]]]:
    """Get the ExPASy closure map."""
    _graph = nx.DiGraph()
    for term in expasy_obo:
        for parent_term in term.parents:
            _graph.add_edge(
                (term.prefix, term.identifier, term.identifier),
                (parent_term.prefix, parent_term.identifier, parent_term.identifier),
            )
        for xref in term.xrefs:
            if xref.prefix in {'uniprot', 'prosite'}:
                _graph.add_edge(
                    (xref.prefix, xref.identifier, xref.name),
                    (term.prefix, term.identifier, term.identifier),
                )

    rv = {
        identifier: list(nx.ancestors(_graph, (db, identifier, name)))
        for db, identifier, name in _graph
        if db == expasy_obo.ontology
    }
    return _graph, rv


def propose_enzyme_modulators() -> pd.DataFrame:
    """Suggest enzyme inhibitors for curation."""
    rv = []

    source_db = 'chebi'
    for term in chebi_obo:
        # Do this as a loop since there is at least one entry that corresponds to several EC codes
        ec_codes = []

        # Requested fix in https://github.com/ebi-chebi/ChEBI/issues/3651
        if term.name == 'EC 1.22* (oxidoreductase acting on halogen in donors) inhibitor':
            modulation = 'inhibitor'
            ec_codes.append('1.22.-.-')

        # Not sure why this has two
        elif term.name == 'EC 1.1.1.34/EC 1.1.1.88 (hydroxymethylglutaryl-CoA reductase) inhibitor':
            modulation = 'inhibitor'
            ec_codes.append('1.1.1.34')
            ec_codes.append('1.1.1.88')

        # Requested rename in https://github.com/ebi-chebi/ChEBI/issues/3653
        elif term.name == 'EC 1.11.1.11 (L-ascorbate peroxidase) inhibitors':
            modulation = 'inhibitor'
            ec_codes.append('1.11.1.11')

        # Requested typo fix in https://github.com/ebi-chebi/ChEBI/issues/3652
        elif term.name == 'EC 3.5.5.1 (nitrilase) inhhibitor':
            modulation = 'inhibitor'
            ec_codes.append('3.5.5.1')

        # All other normal cases
        elif term.name.startswith('EC '):
            if term.name.endswith('inhibitor'):
                modulation = 'inhibitor'
            elif term.name.endswith('activator'):
                modulation = 'activator'
            else:
                logger.warning(f'Unhandled suffix: {term}')
                continue

            ec_code = term.name[len('EC '):].split()[0].replace('*', '-').rstrip('.')

            # Add any remaining dashes
            for _ in range(3 - ec_code.count('.')):
                ec_code += '.-'

            ec_codes.append(ec_code)

        else:
            continue

        rv.extend(
            (
                source_db, term.identifier, term.name, modulation, 'protein family', 'ec-code', ec_code, ec_code,
            )
            for ec_code in ec_codes
        )

    return pd.DataFrame(rv, columns=XREFS_COLUMNS)


def get_enzyme_inhibitor_df(obo: Obo, add_proteins: bool = False) -> pd.DataFrame:
    """Suggest enzyme inhibitors for curation."""
    expasy_graph, ec_code_to_children = get_expasy_closure()
    rv = []

    source_db = 'chebi'
    for term in obo:
        # Do this as a loop since there is at least one entry that corresponds to several EC codes
        ec_codes = []

        # Requested fix in https://github.com/ebi-chebi/ChEBI/issues/3651
        if term.name == 'EC 1.22* (oxidoreductase acting on halogen in donors) inhibitor':
            modulation = 'inhibitor'
            ec_codes.append('1.22.-.-')

        # Not sure why this has two
        elif term.name == 'EC 1.1.1.34/EC 1.1.1.88 (hydroxymethylglutaryl-CoA reductase) inhibitor':
            modulation = 'inhibitor'
            ec_codes.append('1.1.1.34')
            ec_codes.append('1.1.1.88')

        # Requested rename in https://github.com/ebi-chebi/ChEBI/issues/3653
        elif term.name == 'EC 1.11.1.11 (L-ascorbate peroxidase) inhibitors':
            modulation = 'inhibitor'
            ec_codes.append('1.11.1.11')

        # Requested typo fix in https://github.com/ebi-chebi/ChEBI/issues/3652
        elif term.name == 'EC 3.5.5.1 (nitrilase) inhhibitor':
            modulation = 'inhibitor'
            ec_codes.append('3.5.5.1')

        # All other normal cases
        elif term.name.startswith('EC '):
            if term.name.endswith('inhibitor'):
                modulation = 'inhibitor'
            elif term.name.endswith('activator'):
                modulation = 'activator'
            else:
                logger.warning(f'Unhandled suffix: {term}')
                continue

            ec_code = term.name[len('EC '):].split()[0].replace('*', '-').rstrip('.')

            # Add any remaining dashes
            for _ in range(3 - ec_code.count('.')):
                ec_code += '.-'

            ec_codes.append(ec_code)

        else:
            continue

        for ec_code in ec_codes:
            rv.append((
                source_db, term.identifier, term.name, modulation, 'protein family', 'ec-code', ec_code, ec_code,
            ))
            if not add_proteins:
                continue

            children_ec_codes = ec_code_to_children.get(ec_code)
            if children_ec_codes is None:
                print(f'could not find {ec_code} (for {term})')
                continue

            for target_db, target_id, target_name in children_ec_codes:
                target_type = DB_TO_TYPE[target_db]
                rv.append((
                    source_db, term.identifier, term.name, modulation, target_type, target_db, target_id, target_name,
                ))

    return pd.DataFrame(rv, columns=XREFS_COLUMNS)


DB_TO_TYPE = {
    'ec-code': 'protein family',
    'uniprot': 'protein',
    'prosite': 'protein',
}


def suggest_pathway_inhibitor_curation() -> None:
    """Suggest pathway inhibitors for curation."""
    curated_chebi_ids = _get_curated_chebi_ids()

    reclassify_df = _get_inhibitors_reclassification()
    reclassify_chebi_ids = set(reclassify_df.chebi_id)

    print(f'Children of {PATHWAY_INHIBITOR_CHEBI_ID} ({chebi_id_to_name[PATHWAY_INHIBITOR_CHEBI_ID]})')
    for chebi_id in chebi_obo.descendants(PATHWAY_INHIBITOR_CHEBI_ID):
        if any(chebi_id in group for group in (curated_chebi_ids, reclassify_chebi_ids, BLACKLIST)):
            continue  # we already curated this!
        name = chebi_id_to_name[chebi_id]
        if name.endswith('inhibitor') and not name.startswith('EC '):
            print('chebi', f'CHEBI:{chebi_id}', name, 'inhibitor', '?', '?', '?', '?', sep='\t')
            results = post_gilda(name[:-len(' inhibitor')]).json()
            if results:
                print(json.dumps(results, indent=2))


def suggest_inhibitor_curation() -> None:
    """Suggest inhibitors for curation."""
    chebi_ids = (
        chebi_obo.descendants(INHIBITOR_CHEBI_ID)
        - chebi_obo.descendants(PATHWAY_INHIBITOR_CHEBI_ID)
        - chebi_obo.descendants(ENZYME_INHIBITOR_CHEBI_ID)
    )
    for t in _suggest_xrefs_curation(chebi_ids=chebi_ids, suffix='inhibitor'):
        print(*t, sep='\t')


def suggest_agonist_curation() -> None:
    """Suggest agonists for curation."""
    _single_suggest(AGONIST_CHEBI_ID, 'agonist')


def suggest_antagonist_curation() -> None:
    """Suggest antagonists for curation."""
    _single_suggest(ANTAGONIST_CHEBI_ID, 'antagonist')


def suggest_inverse_agonist_curation() -> None:
    """Suggest inverse agonists for curation."""
    _single_suggest(INVERSE_AGONIST_CHEBI_ID, 'inverse agonist')


def suggest_activator_curation() -> None:
    """Suggest activators for curation."""
    _single_suggest(BIOCHEMICAL_ROLE_CHEBI_ID, 'activator')


def _single_suggest(chebi_id: str, suffix, file=None) -> None:
    for t in _suggest_xrefs_curation(suffix=suffix, chebi_ids=chebi_obo.descendants(chebi_id)):
        print(*t, sep='\t', file=file)


def _get_curated_chebi_ids():
    xrefs_df = get_xrefs_df()
    return {
        source_id[len('CHEBI:'):]
        for source_db, source_id in xrefs_df[['source_db', 'source_id']].values
        if source_db == 'chebi'
    }


def _suggest_xrefs_curation(
    *,
    suffix: str,
    chebi_ids: Iterable[str],
    show_missing: bool = False
) -> Iterable[Tuple[str, str, str, str, str, str]]:
    """Suggest curation.

    :param suffix: If the term's name doesn't end with this, skip it
    """
    curated_chebi_ids = _get_curated_chebi_ids()

    for chebi_id in chebi_ids:
        if chebi_id in curated_chebi_ids:
            continue
        name = chebi_id_to_name[chebi_id]
        if not name.casefold().endswith(suffix.casefold()):
            continue
        results = post_gilda(name[:-len(suffix)].rstrip()).json()
        if results:
            for result in results:
                term = result["term"]
                yield (
                    'chebi', f'CHEBI:{chebi_id}', name,
                    suffix,
                    term['db'].lower(), term['id'], term['entry_name'],
                )
        elif show_missing:
            yield 'chebi', chebi_id, name, suffix, '?', '?', '?', '?'


def get_relations_df(graph: nx.MultiDiGraph) -> pd.DataFrame:
    """Assemble the relations dataframe."""
    xrefs_df = get_xrefs_df()

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

    xrefs = defaultdict(list)
    for source_db, source_id, _, modulation, target_type, target_db, target_id, target_name in xrefs_df.values:
        if target_db == 'hgnc':
            for uniprot_id, uniprot_name in get_uniprot_id_names(target_id):
                xrefs[source_id].append((modulation, 'protein', 'uniprot', uniprot_id, uniprot_name))

        elif target_db == 'fplx':
            xrefs[source_id].append((modulation, target_type, target_db, target_id, target_name))

            for hgnc_id, hgnc_symbol in famplex_id_to_members.get(target_id, []):
                xrefs[source_id].append((modulation, 'protein', 'hgnc', hgnc_id, hgnc_symbol))
                for uniprot_id, uniprot_name in get_uniprot_id_names(hgnc_id):
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


def get_uniprot_id_names(hgnc_id: str) -> Iterable[Tuple[str, str]]:
    """Get all of the UniProt identifiers for a given gene."""
    try:
        r = hgnc_id_to_up[str(hgnc_id)]
    except KeyError:
        _k, _v = list(hgnc_id_to_up.items())[0]
        print(f'could not find {hgnc_id} ({type(hgnc_id)} in dict. Example: {_k} ({type(_k)}), {_v} ({type(_v)})')
        raise

    for _uniprot_id in r.split(', '):
        yield _uniprot_id, uniprot_client.get_mnemonic(_uniprot_id)


@click.command()
@click.option('-s', '--suggest', is_flag=True)
@click.option('-b', '--debug', is_flag=True)
def main(suggest: bool, debug: bool) -> None:
    """Run the ChEBI curation pipeline."""
    level = logging.DEBUG if debug else logging.INFO
    logger.setLevel(level)
    logging.basicConfig(level=level)

    # Sort the curated xrefs file
    sort_xrefs_df()

    if True:
        suggest_activator_curation()
        suggest_pathway_inhibitor_curation()
        suggest_inhibitor_curation()
        suggest_agonist_curation()
        suggest_antagonist_curation()
        suggest_inverse_agonist_curation()
        propose_enzyme_modulators()
        return sys.exit(0)

    df = get_relations_df().drop_duplicates()

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

# -*- coding: utf-8 -*-

"""A script to run inference and generate more relationships."""

import itertools as itt
import logging
import os
import time
from collections import defaultdict
from typing import Iterable, List, Mapping, Tuple

import matplotlib.pyplot as plt
import networkx as nx
import pandas as pd
import pybel
import pyobo
import seaborn as sns
from protmapper import uniprot_client
from protmapper.api import hgnc_id_to_up, hgnc_name_to_id
from pybel import BELGraph, dsl
from pyobo.sources.chebi import get_chebi_role_to_children
from pyobo.sources.expasy import get_ec2go, get_obo as get_expasy_obo
from pyobo.struct.typedef import has_member
from tabulate import tabulate
from tqdm import tqdm

from utils import EXPORT_DIRECTORY, HERE, XREFS_COLUMNS, get_xrefs_df

logger = logging.getLogger(__name__)

RELATIONS_OUTPUT_PATH = os.path.join(EXPORT_DIRECTORY, 'relations.tsv')
RELATIONS_SLIM_OUTPUT_PATH = os.path.join(EXPORT_DIRECTORY, 'relations_slim.tsv')

FAMPLEX_RELATIONS_URL = 'https://raw.githubusercontent.com/sorgerlab/famplex/master/relations.csv'
FAMPLEX_EQUIVALENCES_URL = 'https://raw.githubusercontent.com/sorgerlab/famplex/master/equivalences.csv'
FAMPLEX_HGNC_SYMBOL_MAP_URL = 'https://raw.githubusercontent.com/sorgerlab/famplex/master/export/hgnc_symbol_map.csv'

EXPASY_OBO = get_expasy_obo()

DB_TO_TYPE = {
    'ec-code': 'protein family',
    'uniprot': 'protein',
    'prosite': 'protein',
}

sns.set(font_scale=1.3, style='whitegrid')


def get_expasy_closure() -> Tuple[nx.DiGraph, Mapping[str, List[str]]]:
    """Get the ExPASy closure map."""
    _graph = nx.DiGraph()
    for term in EXPASY_OBO:
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
        if db == EXPASY_OBO.ontology
    }
    return _graph, rv


def get_relations_df(use_sub_roles=False) -> pd.DataFrame:
    """Assemble the relations dataframe."""
    xrefs_df = get_xrefs_df()

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

    logger.info('getting enzyme classes')
    expasy_graph, ec_code_to_children = get_expasy_closure()
    logger.info('getting ec2go')
    ec2go = get_ec2go()

    x = defaultdict(list)
    it = tqdm(
        xrefs_df.values,
        total=len(xrefs_df.index),
        desc='inferring over target hierarchies',
    )
    for source_db, source_id, _, modulation, target_type, target_db, target_id, target_name in it:
        if source_db != 'chebi':
            continue

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

        elif target_db == 'ec-code':
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

    logger.info('inferring over role hiearchies')
    db_to_role_to_chemical_curies = {
        'chebi': get_chebi_role_to_children(),
    }

    rows = []
    for (role_db, role_id), entries in tqdm(x.items(), desc='inferring over role hierarchies'):
        sub_role_curies = {(role_db, role_id)}

        if role_db == 'chebi' and use_sub_roles:
            sub_role_curies |= {
                pyobo.normalize_curie(c)
                for c in pyobo.get_subhierarchy(role_db, role_id)
            }

        for modulation, target_type, target_db, target_id, target_name in entries:
            chemical_curies = set(itt.chain.from_iterable(
                db_to_role_to_chemical_curies[sub_role_db].get(sub_role_id, [])
                for sub_role_db, sub_role_id in sub_role_curies
            ))
            if not chemical_curies:
                logger.debug('no inference for %s:%s', role_db, role_id)
                continue
            for chemical_db, chemical_id in chemical_curies:
                rows.append((
                    chemical_db, chemical_id, pyobo.get_name(chemical_db, chemical_id),
                    modulation, target_type, target_db, target_id, target_name,
                ))
    return pd.DataFrame(rows, columns=XREFS_COLUMNS)


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


def rewrite_repo_readme():
    """Rewrite the summary of curated content in the repository's readme, automatically."""
    df = get_xrefs_df()

    summary_df = df.groupby(['source_db', 'modulation', 'type', 'target_db']).size().reset_index()
    summary_df.columns = ['Source Database', 'Modulation', 'Target Type', 'Target Database', 'Count']
    summary_df.to_csv(os.path.join(EXPORT_DIRECTORY, 'curated_summary.tsv'), sep='\t', index=False)

    modulation_summary_df = df.groupby('modulation').size().reset_index()
    modulation_summary_df.columns = ['Modulation', 'Count']
    modulation_summary_df.to_csv(
        os.path.join(EXPORT_DIRECTORY, 'curated_summary_by_modulation.tsv'), sep='\t', index=False,
    )
    type_summary_df = df.groupby('type').size().reset_index()
    type_summary_df.columns = ['Target Type', 'Count']
    type_summary_df.to_csv(
        os.path.join(EXPORT_DIRECTORY, 'curated_summary_by_type.tsv'), sep='\t', index=False,
    )
    namespace_summary_df = df.groupby('target_db').size().reset_index()
    namespace_summary_df.columns = ['Target Database', 'Count']
    namespace_summary_df.to_csv(
        os.path.join(EXPORT_DIRECTORY, 'curated_summary_by_namespace.tsv'), sep='\t', index=False,
    )

    logger.info('Plotting modulation and target type summary')
    fig, (lax, rax) = plt.subplots(nrows=1, ncols=2, figsize=(12, 5))

    modulation_summary_df['Modulation'] = modulation_summary_df['Modulation'].map(str.title)
    g = sns.barplot(y="Modulation", x='Count', data=modulation_summary_df, ax=lax)
    g.set_xscale("log")
    lax.set_title(
        f'Modulation ({len(df.index)} in {len(modulation_summary_df.index)} relations)',
        fontdict={'fontweight': 'bold'},
    )
    lax.set_ylabel('')

    type_summary_df['Target Type'] = type_summary_df['Target Type'].map(str.title)
    g = sns.barplot(y="Target Type", x='Count', data=type_summary_df, ax=rax)
    g.set_xscale("log")
    rax.set_title(
        f'Target Type ({len(df.index)} in {len(type_summary_df.index)} types)',
        fontdict={'fontweight': 'bold'},
    )
    rax.set_ylabel('')

    plt.tight_layout()
    plt.savefig(os.path.join(EXPORT_DIRECTORY, 'curated_summary.png'), dpi=300)

    text = f'There are {len(df.index)} curated roles as of export on {time.asctime()}\n\n'
    text += tabulate(modulation_summary_df.values, ['Modulation', 'Count'], tablefmt='rst')
    text += '\n\n'
    text += tabulate(type_summary_df.values, ['Target Entity Type', 'Count'], tablefmt='rst')
    text += '\n\n'
    text += tabulate(namespace_summary_df.values, ['Target Database', 'Count'], tablefmt='rst')
    text += '\n'

    path = os.path.join(HERE, 'README.rst')
    with open(path) as file:
        readme = [line.rstrip() for line in file]

    for i, line in enumerate(readme):
        if line == 'Summary':
            start = i + 2
            break
    else:
        raise ValueError('could not find summary block')

    for i, line in enumerate(readme):
        if line == 'Axioms':
            end = i
            break
    else:
        raise ValueError('could not find end block')

    with open(path, 'w') as file:
        for line in readme[:start]:
            print(line, file=file)
        print(text, file=file)
        for line in readme[end:]:
            print(line, file=file)


def write_export():
    """Generate export TSVs.

    1. Full TSV at ``export/relations.tsv``
    2. Slim TSV at ``export/relations_slim.tsv``, appropriate for machine learning
    """
    df = get_relations_df().drop_duplicates()
    logger.info('got relations df with %s rows', len(df.index))

    columns = [
        'modulation', 'target_type', 'source_db', 'source_id', 'source_name', 'target_db', 'target_id', 'target_name',
    ]
    df[columns].sort_values(columns).to_csv(RELATIONS_OUTPUT_PATH, sep='\t', index=False)

    logger.info('outputting slim df to %s', RELATIONS_SLIM_OUTPUT_PATH)
    slim_columns = ['source_db', 'source_id', 'modulation', 'target_db', 'target_id']
    df[slim_columns].sort_values(slim_columns).to_csv(RELATIONS_SLIM_OUTPUT_PATH, sep='\t', index=False)

    logger.info('making summary df')
    summary_df = df.groupby(['source_db', 'modulation', 'target_type', 'target_db']).size().reset_index()
    summary_df.columns = ['Source Database', 'Modulation', 'Target Type', 'Target Database', 'Count']
    summary_df.to_csv(os.path.join(EXPORT_DIRECTORY, 'inferred_summary.tsv'), sep='\t', index=False)
    summary_df_str = tabulate(
        summary_df.values,
        ['source_db', 'relation', 'target_type', 'target_db', 'count'],
        tablefmt='github',
    )

    logger.info('making modulation summary df')
    modulation_summary_df = df.groupby(['modulation']).size().reset_index()
    modulation_summary_df.columns = ['Modulation', 'Count']
    modulation_summary_df.to_csv(
        os.path.join(EXPORT_DIRECTORY, 'inferred_summary_by_modulation.tsv'), sep='\t', index=False,
    )
    modulation_str = tabulate(
        modulation_summary_df.values,
        ['modulation', 'count'],
        tablefmt='github',
    )

    logger.info('making namespace summary df')
    namespace_summary_df = df.groupby(['target_db']).size().reset_index()
    namespace_summary_df.columns = ['Target Database', 'Count']
    namespace_summary_df.to_csv(
        os.path.join(EXPORT_DIRECTORY, 'inferred_summary_by_namespace.tsv'), sep='\t', index=False,
    )
    ns_str = tabulate(
        namespace_summary_df.values,
        ['namespace', 'count'],
        tablefmt='github',
    )

    logger.info('making type summary df')
    type_summary_df = df.groupby(['target_type']).size().reset_index()
    type_summary_df.columns = ['Target Type', 'Count']
    type_summary_df.to_csv(os.path.join(EXPORT_DIRECTORY, 'inferred_summary_by_type.tsv'), sep='\t', index=False)
    type_str = tabulate(
        type_summary_df.values,
        ['type', 'count'],
        tablefmt='github',
    )

    logger.info('Plotting modulation and target type summary')
    fig, (lax, rax) = plt.subplots(nrows=1, ncols=2, figsize=(12, 5))

    modulation_summary_df['Modulation'] = modulation_summary_df['Modulation'].map(str.title)
    g = sns.barplot(y="Modulation", x='Count', data=modulation_summary_df, ax=lax)
    g.set_xscale("log")
    lax.set_title(
        f'Modulation ({len(df.index)} in {len(modulation_summary_df.index)} relations)',
        fontdict={'fontweight': 'bold'},
    )
    lax.set_ylabel('')

    type_summary_df['Target Type'] = type_summary_df['Target Type'].map(str.title)
    g = sns.barplot(y="Target Type", x='Count', data=type_summary_df, ax=rax)
    g.set_xscale("log")
    rax.set_title(
        f'Target Type ({len(df.index)} in {len(type_summary_df.index)} types)',
        fontdict={'fontweight': 'bold'},
    )
    rax.set_ylabel('')

    plt.tight_layout()
    plt.savefig(os.path.join(EXPORT_DIRECTORY, 'inferred_summary.png'), dpi=300)

    with open(os.path.join(EXPORT_DIRECTORY, 'README.md'), 'w') as file:
        print('# Export Summary\n', file=file)
        print(f'Exported {len(df.index)} relations on {time.asctime()}\n', file=file)
        print('\n## Summary by Modulation\n', file=file)
        print(modulation_str, file=file)
        print('\n## Summary by Type\n', file=file)
        print(type_str, file=file)
        print('\n## Summary by Namespace\n', file=file)
        print(ns_str, file=file)
        print('\n## Relation Summary\n', file=file)
        print(summary_df_str, file=file)

    write_bel_export(df)


_type_map = {
    'biological process': dsl.BiologicalProcess,
    'chemical': dsl.Abundance,
    'organism': dsl.Population,
    'phenotype': dsl.Pathology,
    'protein': dsl.Protein,
    'protein family': dsl.Protein,
    'protein complex': dsl.NamedComplexAbundance,
}

_adders = {
    'activator': BELGraph.add_activates,
    'agonist': BELGraph.add_activates,
    'antagonist': BELGraph.add_inhibits,
    'inhibitor': BELGraph.add_inhibits,
    'inverse agonist': BELGraph.add_activates,
    'modulator': BELGraph.add_regulates,
}


def write_bel_export(df: pd.DataFrame):
    """Write BEL export."""
    graph = BELGraph(name='Chemical Roles')
    it = tqdm(df.dropna().values, total=len(df.index), desc='mapping to BEL', unit_scale=True)
    for source_db, source_id, source_name, modulation, target_type, target_db, target_id, target_name in it:
        if target_type == 'molecular function':
            continue
        source = pybel.dsl.Abundance(
            namespace=source_db,
            identifier=source_id,
            name=source_name,
        )
        target = _type_map[target_type](
            namespace=target_db,
            identifier=target_id,
            name=target_name,
        )
        adder = _adders[modulation]
        adder(graph, source, target, citation='', evidence='')

    pybel.dump(graph, os.path.join(EXPORT_DIRECTORY, 'export.bel.nodelink.json.gz'), indent=2)


def main():
    """Rewrite readme and generate new export."""
    rewrite_repo_readme()
    write_export()


if __name__ == '__main__':
    main()

# -*- coding: utf-8 -*-

"""A script to run inference and generate more relationships."""

import logging
import os
import time

import matplotlib.pyplot as plt
import seaborn as sns
from tabulate import tabulate

from chemical_roles.constants import DATA, DOCS, IMG, RELATIONS_OUTPUT_PATH, RELATIONS_SLIM_OUTPUT_PATH, ROOT
from chemical_roles.export.utils import get_relations_df
from chemical_roles.resources import get_xrefs_df

logger = logging.getLogger(__name__)

FAMPLEX_EQUIVALENCES_URL = 'https://raw.githubusercontent.com/sorgerlab/famplex/master/equivalences.csv'
FAMPLEX_HGNC_SYMBOL_MAP_URL = 'https://raw.githubusercontent.com/sorgerlab/famplex/master/export/hgnc_symbol_map.csv'


def rewrite_repo_readme():
    """Rewrite the summary of curated content in the repository's readme, automatically."""
    df = get_xrefs_df()

    summary_df = df.groupby(['source_db', 'modulation', 'type', 'target_db']).size().reset_index()
    summary_df.columns = ['Source Database', 'Modulation', 'Target Type', 'Target Database', 'Count']
    summary_df.to_csv(os.path.join(DATA, 'curated_summary.tsv'), sep='\t', index=False)

    modulation_summary_df = df.groupby('modulation').size().reset_index()
    modulation_summary_df.columns = ['Modulation', 'Count']
    modulation_summary_df.to_csv(
        os.path.join(DATA, 'curated_summary_by_modulation.tsv'), sep='\t', index=False,
    )
    type_summary_df = df.groupby('type').size().reset_index()
    type_summary_df.columns = ['Target Type', 'Count']
    type_summary_df.to_csv(
        os.path.join(DATA, 'curated_summary_by_type.tsv'), sep='\t', index=False,
    )
    namespace_summary_df = df.groupby('target_db').size().reset_index()
    namespace_summary_df.columns = ['Target Database', 'Count']
    namespace_summary_df.to_csv(
        os.path.join(DATA, 'curated_summary_by_namespace.tsv'), sep='\t', index=False,
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
    plt.savefig(os.path.join(IMG, 'curated_summary.png'), dpi=300)
    plt.savefig(os.path.join(IMG, 'curated_summary.svg'))

    text = f'There are {len(df.index)} curated roles as of export on {time.asctime()}\n\n'
    text += tabulate(modulation_summary_df.values, ['Modulation', 'Count'], tablefmt='rst')
    text += '\n\n'
    text += tabulate(type_summary_df.values, ['Target Entity Type', 'Count'], tablefmt='rst')
    text += '\n\n'
    text += tabulate(namespace_summary_df.values, ['Target Database', 'Count'], tablefmt='rst')
    text += '\n'

    readme_path = os.path.join(ROOT, 'README.rst')
    with open(readme_path) as file:
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

    with open(readme_path, 'w') as file:
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
    summary_df.to_csv(os.path.join(DATA, 'inferred_summary.tsv'), sep='\t', index=False)
    summary_df_str = tabulate(
        summary_df.values,
        ['source_db', 'relation', 'target_type', 'target_db', 'count'],
        tablefmt='github',
    )

    logger.info('making modulation summary df')
    modulation_summary_df = df.groupby(['modulation']).size().reset_index()
    modulation_summary_df.columns = ['Modulation', 'Count']
    modulation_summary_df.to_csv(
        os.path.join(DATA, 'inferred_summary_by_modulation.tsv'), sep='\t', index=False,
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
        os.path.join(DATA, 'inferred_summary_by_namespace.tsv'), sep='\t', index=False,
    )
    ns_str = tabulate(
        namespace_summary_df.values,
        ['namespace', 'count'],
        tablefmt='github',
    )

    logger.info('making type summary df')
    type_summary_df = df.groupby(['target_type']).size().reset_index()
    type_summary_df.columns = ['Target Type', 'Count']
    type_summary_df.to_csv(os.path.join(DATA, 'inferred_summary_by_type.tsv'), sep='\t', index=False)
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
    plt.savefig(os.path.join(IMG, 'inferred_summary.png'), dpi=300)
    plt.savefig(os.path.join(IMG, 'inferred_summary.svg'))

    with open(os.path.join(DOCS, 'index.md'), 'w') as file:
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

# -*- coding: utf-8 -*-

"""Export to BEL."""

import click
import pybel
from pybel import BELGraph, dsl
from tqdm import tqdm

from .utils import get_relations_df

__all__ = [
    'get_bel',
]

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


def get_bel() -> BELGraph:
    """Get Chemical Roles as BEL."""
    df = get_relations_df().drop_duplicates()
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
    return graph


@click.command()
@click.argument('path')
def bel(path):
    """Write BEL export."""
    graph = get_bel()
    pybel.dump(graph, path)


if __name__ == '__main__':
    bel()

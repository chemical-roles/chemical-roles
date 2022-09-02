# -*- coding: utf-8 -*-

"""Export to OBO."""

from typing import Iterable, Mapping, Optional, Tuple

from pyobo import Obo, Reference, Term, TypeDef
from tqdm import tqdm

from .utils import get_relations_df

__all__ = [
    "get_obo",
]


def get_obo() -> Obo:
    """Get Chemical Roles as OBO."""
    return Obo(
        name="Chemical Roles Graph",
        ontology="crog",
        iter_terms=iter_terms,
    )


def iter_terms() -> Iterable[Term]:
    df = get_relations_df()
    it = tqdm(df.dropna().values, total=len(df.index), desc="mapping to OBO", unit_scale=True)
    ref_term = {}
    for (
        source_db,
        source_id,
        source_name,
        modulation,
        target_type,
        target_db,
        target_id,
        target_name,
    ) in it:
        source = Reference(source_db.upper(), source_id, source_name)
        term = ref_term.get(source)
        if term is None:
            term = ref_term[source] = Term(reference=source)

        typedef = _get_typedef(target_db=target_db, target_type=target_type, modulation=modulation)
        if typedef is not None:
            term.append_relationship(typedef, Reference(target_db.upper(), target_id, target_name))

    for term in ref_term.values():
        if len(term.relationships) > 0:
            yield term


# TODO fill out rest
_typedefs: Mapping[Tuple[str, str, str], TypeDef] = {
    ("go", "biological process", "activator"): TypeDef(
        Reference("RO", "0002213", "positively regulates")
    ),
    ("go", "biological process", "inhibitor"): TypeDef(
        Reference("RO", "0002212", "negatively regulates")
    ),
}

_logged = set()


def _get_typedef(target_db, target_type, modulation) -> Optional[TypeDef]:
    t = (target_db, target_type, modulation)
    rv = _typedefs.get(t)
    if rv is not None:
        return rv
    if t not in _logged:
        _logged.add(t)
        tqdm.write(f"no strategy for: {target_db} {target_type} {modulation}")

# -*- coding: utf-8 -*-

"""Export to OBO."""

from typing import Iterable, Optional

from pyobo import Annotation, Obo, Reference, Term, TypeDef, build_ontology
from pyobo.struct.typedef import (
    directly_regulates_activity_of,
    has_contributor,
    is_agonist_of,
    is_antagonist_of,
    is_inverse_agonist_of,
    negatively_regulates,
    positively_regulates,
)
from pyobo.struct.vocabulary import CHARLIE
from tqdm import tqdm

from .utils import get_relations_df

__all__ = [
    "get_obo",
]


KEY_TO_TYPEDEF: dict[str, TypeDef] = {
    "activator": positively_regulates,
    "inhibitor": negatively_regulates,
    "agonist": is_agonist_of,
    "inverse agonist": is_inverse_agonist_of,
    "antagonist": is_antagonist_of,
    "modulator": directly_regulates_activity_of,
}


def get_obo(*, use_inferred: bool, use_sub_roles: bool) -> Obo:
    """Get Chemical Roles as OBO."""
    return build_ontology(
        prefix="crog",
        name="Chemical Roles Graph",
        terms=list(iter_terms(use_inferred=use_inferred, use_sub_roles=use_sub_roles)),
        typedefs=[has_contributor, *KEY_TO_TYPEDEF.values()],
    )


def iter_terms(*, use_inferred: bool, use_sub_roles: bool) -> Iterable[Term]:
    df = get_relations_df(use_inferred=use_inferred, use_sub_roles=use_sub_roles)
    it = tqdm(df.dropna().values, total=len(df.index), desc="mapping to OBO", unit_scale=True)
    ref_term: dict[Reference, Term] = {}
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
        source = Reference(prefix=source_db.upper(), identifier=source_id)
        term = ref_term.get(source)
        if term is None:
            term = ref_term[source] = Term(reference=source)

        predicate = _predicate_from_modulation(modulation)
        if predicate is not None:
            term.annotate_object(
                predicate,
                Reference(prefix=target_db.upper(), identifier=target_id, name=target_name),
                annotations=[Annotation(predicate=has_contributor.reference, value=CHARLIE)],
            )

    for term in ref_term.values():
        if len(term.relationships) > 0:
            yield term


_logged: set[str] = set()


def _predicate_from_modulation(modulation: str) -> Optional[TypeDef]:
    if tt := KEY_TO_TYPEDEF.get(modulation):
        return tt
    if modulation not in _logged:
        _logged.add(modulation)
        tqdm.write(f"no strategy for: {modulation}")
    return None

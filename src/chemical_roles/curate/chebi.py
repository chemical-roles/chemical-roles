# -*- coding: utf-8 -*-

"""A script to help curate ChEBI relations and infer new ones."""

import itertools as itt
import json
import logging
from functools import lru_cache
from typing import Iterable, Optional, Set, TextIO

import click
import pandas as pd
import pyobo
from curies import Reference
from more_click import verbose_option
from pyobo.sources.expasy import get_ec2go
from tqdm import tqdm

from ..resources import (
    RECLASSIFICATION_PATH,
    UNCURATED_CHEBI_PATH,
    get_blacklist_roles_df,
    get_irrelevant_roles_df,
    get_xrefs_df,
)
from ..utils import (
    SUFFIXES,
    XREFS_COLUMNS,
    GildaTuple,
    post_gilda,
    sort_xrefs_df,
    yield_gilda,
)

logger = logging.getLogger(__name__)

BIOLOGICAL_ROLE_ID = "24432"
APPLICATION_ROLE_ID = "33232"
BIOCHEMICAL_ROLE_CHEBI_ID = "52206"
PATHWAY_INHIBITOR_CHEBI_ID = "76932"
ENZYME_INHIBITOR_CHEBI_ID = "23924"
AGONIST_CHEBI_ID = "48705"
INVERSE_AGONIST_CHEBI_ID = "90847"
INHIBITOR_CHEBI_ID = "35222"
ANTAGONIST_CHEBI_ID = "48706"


@lru_cache(maxsize=1)
def get_curated_role_chebi_ids() -> Set[str]:
    """Get the curated role chebi id set."""
    xrefs_df = get_xrefs_df()
    rv = {
        source_id
        for source_db, source_id in xrefs_df[["source_db", "source_id"]].values
        if source_db == "chebi"
    }
    logger.info("%d pre-curated ChEBI role identifiers", len(rv))
    return rv


@lru_cache(maxsize=1)
def _get_irrelevant_role_chebi_references() -> Set[Reference]:
    references = set(
        itt.chain.from_iterable(
            pyobo.get_descendants("chebi", chebi_id)
            for chebi_id in get_irrelevant_roles_df().identifier
        )
    )
    logger.info("%d irrelevant ChEBI role identifiers", len(references))
    return references


@lru_cache(maxsize=1)
def _get_blacklist_chebi_ids() -> Set[str]:
    rv = set(get_blacklist_roles_df().identifier)
    logger.info("%d blacklisted ChEBI role identifiers", len(rv))
    return rv


def _get_inhibitors_reclassification() -> pd.DataFrame:
    return pd.read_csv(RECLASSIFICATION_PATH, sep="\t", comment="#")


def propose_enzyme_modulators() -> pd.DataFrame:
    """Suggest enzyme inhibitors for curation."""
    ec2go = get_ec2go()
    rv = []

    for identifier, name in pyobo.get_id_name_mapping("chebi").items():
        # Do this as a loop since there is at least one entry that corresponds to several EC codes
        ec_codes = []

        # Requested fix in https://github.com/ebi-chebi/ChEBI/issues/3651
        if name == "EC 1.22* (oxidoreductase acting on halogen in donors) inhibitor":
            modulation = "inhibitor"
            ec_codes.append("1.22.-.-")

        # Not sure why this has two
        elif name == "EC 1.1.1.34/EC 1.1.1.88 (hydroxymethylglutaryl-CoA reductase) inhibitor":
            modulation = "inhibitor"
            ec_codes.append("1.1.1.34")
            ec_codes.append("1.1.1.88")

        # Requested rename in https://github.com/ebi-chebi/ChEBI/issues/3653
        elif name == "EC 1.11.1.11 (L-ascorbate peroxidase) inhibitors":
            modulation = "inhibitor"
            ec_codes.append("1.11.1.11")

        # Requested typo fix in https://github.com/ebi-chebi/ChEBI/issues/3652
        elif name == "EC 3.5.5.1 (nitrilase) inhhibitor":
            modulation = "inhibitor"
            ec_codes.append("3.5.5.1")

        # All other normal cases
        elif name.startswith("EC "):
            if name.endswith("inhibitor"):
                modulation = "inhibitor"
            elif name.endswith("activator"):
                modulation = "activator"
            else:
                logger.warning(f"Unhandled suffix: {identifier} ! {name}")
                continue

            ec_code = name[len("EC ") :].split()[0].replace("*", "-").rstrip(".")

            # Add any remaining dashes
            for _ in range(3 - ec_code.count(".")):
                ec_code += ".-"

            ec_codes.append(ec_code)

        else:
            continue

        for ec_code in ec_codes:
            rv.append(
                (
                    "chebi",
                    identifier,
                    name,
                    modulation,
                    "protein family",
                    "ec-code",
                    ec_code,
                    pyobo.get_name("eccode", ec_code) or ec_code,
                )
            )

            for go_activity_id, go_activity_name in ec2go.get(ec_code, []):
                rv.append(
                    (
                        "chebi",
                        identifier,
                        name,
                        modulation,
                        "activity",
                        "go",
                        go_activity_id,
                        go_activity_name,
                    )
                )

    return pd.DataFrame(rv, columns=XREFS_COLUMNS)


def suggest_pathway_inhibitor_curation() -> None:
    """Suggest pathway inhibitors for curation."""
    reclassify_df = _get_inhibitors_reclassification()
    reclassify_chebi_ids = set(reclassify_df.chebi_id)

    print(
        f"Children of {PATHWAY_INHIBITOR_CHEBI_ID} ({pyobo.get_name('chebi', PATHWAY_INHIBITOR_CHEBI_ID)})"
    )
    chebi_referenes = pyobo.get_descendants("chebi", PATHWAY_INHIBITOR_CHEBI_ID)
    for chebi_reference in chebi_referenes:
        if any(
            chebi_reference.identifier in group
            for group in (
                get_curated_role_chebi_ids(),
                reclassify_chebi_ids,
                _get_blacklist_chebi_ids(),
            )
        ):
            continue  # we already curated this!
        name = pyobo.get_name("chebi", chebi_reference)
        if name is None:
            logger.warning("could not find chebi:%s", chebi_reference)
            raise KeyError(chebi_reference.curie)
        if name.endswith("inhibitor") and not name.startswith("EC "):
            print(
                "chebi", chebi_reference.identifier, name, "inhibitor", "?", "?", "?", "?", sep="\t"
            )
            results = post_gilda(name[: -len(" inhibitor")]).json()
            if results:
                print(json.dumps(results, indent=2))


def suggest_inhibitor_curation() -> None:
    """Suggest inhibitors for curation."""
    chebi_references: set[Reference] = (
        pyobo.get_descendants("chebi", INHIBITOR_CHEBI_ID)
        - pyobo.get_descendants("chebi", PATHWAY_INHIBITOR_CHEBI_ID)
        - pyobo.get_descendants("chebi", ENZYME_INHIBITOR_CHEBI_ID)
    )
    for t in _suggest_xrefs_curation(chebi_references=chebi_references, suffix="inhibitor"):
        print(*t, sep="\t")


def suggest_agonist_curation() -> None:
    """Suggest agonists for curation."""
    _single_suggest(AGONIST_CHEBI_ID, "agonist")


def suggest_antagonist_curation() -> None:
    """Suggest antagonists for curation."""
    _single_suggest(ANTAGONIST_CHEBI_ID, "antagonist")


def suggest_inverse_agonist_curation() -> None:
    """Suggest inverse agonists for curation."""
    _single_suggest(INVERSE_AGONIST_CHEBI_ID, "inverse agonist")


def suggest_activator_curation() -> None:
    """Suggest activators for curation."""
    _single_suggest(BIOCHEMICAL_ROLE_CHEBI_ID, "activator")


def suggest_all_roles(show_ungrounded: bool = False, file: Optional[TextIO] = None) -> None:
    """Suggest all roles."""
    logger.info(
        "Getting descendants of chebi:%s and chebi:%s",
        BIOLOGICAL_ROLE_ID,
        APPLICATION_ROLE_ID,
    )
    chebi_references = pyobo.get_descendants("chebi", BIOLOGICAL_ROLE_ID) | pyobo.get_descendants(
        "chebi", APPLICATION_ROLE_ID
    )
    print(*get_xrefs_df().columns, sep="\t", file=file)
    for row in _iter_gilda(chebi_references=chebi_references, show_missing=show_ungrounded):
        print(*row, sep="\t", file=file)


def _single_suggest(
    chebi_reference: str, suffix: str, file: TextIO | None = None, show_missing: bool = False
) -> None:
    descendant_references = pyobo.get_descendants(chebi_reference)
    logger.info(
        "Suggesting for %d descendants of %s ! %s",
        len(descendant_references),
        chebi_reference.curie,
        pyobo.get_name(chebi_reference),
    )
    for t in _suggest_xrefs_curation(
        suffix=suffix, chebi_references=descendant_references, show_missing=show_missing
    ):
        print(*t, sep="\t", file=file)


def _suggest_xrefs_curation(
    *,
    suffix: str,
    chebi_references: Iterable[Reference],
    show_missing: bool = False,
) -> Iterable[GildaTuple]:
    """Suggest curation.

    :param suffix: If the term's name doesn't end with this, skip it
    """
    chebi_references = (
        chebi_reference
        for chebi_reference in chebi_references
        if pyobo.get_name(chebi_reference).casefold().endswith(suffix.casefold())
    )
    yield from _iter_gilda(chebi_references, suffix=suffix, show_missing=show_missing)


def _iter_gilda(
    chebi_references: Iterable[Reference],
    show_missing: bool,
    suffix: Optional[str] = None,
    use_tqdm: bool = True,
) -> Iterable[GildaTuple]:
    for chebi_reference in tqdm(list(chebi_references), desc="making ChEBI curation sheet"):
        if (
            chebi_reference in get_curated_role_chebi_ids()
            or chebi_reference in _get_irrelevant_role_chebi_references()
        ):
            continue  # already curated, skip
        name = pyobo.get_name(chebi_reference)
        if name is None:
            logger.warning("could not look up %s", chebi_reference.curie)
            continue
            # raise KeyError(f'chebi:{chebi_id}')

        if suffix is not None:
            search_text = name[: -len(suffix)].rstrip()
            yield from yield_gilda(
                "chebi", chebi_reference.identifier, name, suffix, search_text, show_missing
            )
        else:
            for _suffix in SUFFIXES:
                if name.endswith(_suffix):
                    search_text = name[: -len(_suffix)].rstrip()
                    yield from yield_gilda(
                        "chebi",
                        chebi_reference.identifier,
                        name,
                        _suffix,
                        search_text,
                        show_missing,
                    )
                    break


@click.command(name="chebi")
@verbose_option
@click.option("--show-ungrounded", is_flag=True)
@click.option("--output", type=click.File("w"), default=UNCURATED_CHEBI_PATH)
def curate_chebi(show_ungrounded: bool, output: Optional[TextIO]) -> None:
    """Run the ChEBI curation pipeline."""
    sort_xrefs_df()
    # suggest_activator_curation()
    # suggest_pathway_inhibitor_curation()
    # suggest_inhibitor_curation()
    # suggest_agonist_curation()
    # suggest_antagonist_curation()
    # suggest_inverse_agonist_curation()
    # propose_enzyme_modulators()
    suggest_all_roles(show_ungrounded=show_ungrounded or output is not None, file=output)


if __name__ == "__main__":
    curate_chebi()

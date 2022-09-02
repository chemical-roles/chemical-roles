# -*- coding: utf-8 -*-

"""Get wikidata mappings."""

# pip install sparqlwrapper
# https://rdflib.github.io/sparqlwrapper/

from SPARQLWrapper import JSON, SPARQLWrapper

WIKIDATA_ENDPOINT_URL = "https://query.wikidata.org/sparql"

protein_query = """SELECT ?gene ?hgnc_id ?hgnc_symbol ?protein ?uniprot_id
WHERE {
?gene wdt:P354 ?hgnc_id ;
      wdt:P353 ?hgnc_symbol ;
      wdt:P688 ?protein .
?protein wdt:P352 ?uniprot_id .
}
ORDER BY ?hgnc_id
"""

chemical_query = """SELECT ?chebi_id ?chemical
WHERE {
?chemical wdt:P683 ?chebi_id .
}
ORDER BY ?chebi_id
"""


def _get_results(endpoint_url, query):
    sparql = SPARQLWrapper(endpoint_url)
    sparql.setQuery(query)
    sparql.setReturnFormat(JSON)
    return sparql.query().convert()


def main():
    """Get wikidata mappings."""
    with open("resources/wd_proteins.tsv", "w") as file:
        print(
            "hgnc_id",
            "hgnc_symbol",
            "gene_wikidata_id",
            "uniprot_id",
            "protein_wikidata_id",
            sep="\t",
            file=file,
        )
        protein_results = _get_results(WIKIDATA_ENDPOINT_URL, protein_query)
        protein_results = (
            (
                int(result["hgnc_id"]["value"]),
                result["hgnc_symbol"]["value"],
                result["gene"]["value"].split("/")[-1],
                result["uniprot_id"]["value"],
                result["protein"]["value"].split("/")[-1],
            )
            for result in protein_results["results"]["bindings"]
        )
        for result in sorted(protein_results):
            print(*result, sep="\t", file=file)

    with open("resources/wd_chemicals.tsv", "w") as file:
        print(
            "chebi_id",
            "wikidata_id",
            sep="\t",
            file=file,
        )
        chemical_results = _get_results(WIKIDATA_ENDPOINT_URL, chemical_query)
        chemical_results = (
            (
                int(result["chebi_id"]["value"]),
                result["chebi_id"]["value"],
                result["chemical"]["value"].split("/")[-1],
            )
            for result in chemical_results["results"]["bindings"]
        )
        for _, chebi_id, wikidata_id in sorted(chemical_results):
            print(chebi_id, wikidata_id, sep="\t", file=file)


if __name__ == "__main__":
    main()

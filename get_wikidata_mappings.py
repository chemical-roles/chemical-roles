# pip install sparqlwrapper
# https://rdflib.github.io/sparqlwrapper/

from SPARQLWrapper import SPARQLWrapper, JSON

endpoint_url = "https://query.wikidata.org/sparql"

protein_query = """SELECT ?gene ?hgnc_id ?hgnc_symbol ?protein ?uniprot_id
WHERE
{
	?gene wdt:P354 ?hgnc_id ;
          wdt:P353 ?hgnc_symbol ;
          wdt:P688 ?protein .
    ?protein wdt:P352 ?uniprot_id .   
}
ORDER BY ?hgnc_id
"""

chemical_query = """SELECT ?chemical ?chebi_id
WHERE
{
	?chemical wdt:P683 ?chebi_id .
}
ORDER BY ?chebi_id
"""


def get_results(endpoint_url, query):
    sparql = SPARQLWrapper(endpoint_url)
    sparql.setQuery(query)
    sparql.setReturnFormat(JSON)
    return sparql.query().convert()


def main():
    with open('resources/wd_proteins.tsv', 'w') as file:
        print(
            'gene_wikidata_id',
            'hgnc_id',
            'hgnc_symbol',
            'protein_wikidata_id',
            'uniprot_id',
            sep='\t',
            file=file,
        )
        protein_results = get_results(endpoint_url, protein_query)
        for result in protein_results["results"]["bindings"]:
            print(
                result['gene']['value'].split('/')[-1],
                result['hgnc_id']['value'],
                result['hgnc_symbol']['value'],
                result['protein']['value'].split('/')[-1],
                result['uniprot_id']['value'],
                sep='\t',
                file=file,
            )
            
    with open('resources/wd_chemicals.tsv', 'w') as file:
        print(
            'wikidata_id',
            'chebi_id',
            sep='\t',
            file=file,
        )
        chemical_results = get_results(endpoint_url, chemical_query)
        for result in chemical_results["results"]["bindings"]:
            print(
                result['chemical']['value'].split('/')[-1],
                f"CHEBI:{result['chebi_id']['value']}",
                sep='\t',
                file=file,
            )


if __name__ == '__main__':
    main()

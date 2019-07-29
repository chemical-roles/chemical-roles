ChEBI Relations
===============
This repository is meant to house links between ChEBI terms describing various
types of inhibitors and their targets in other vocabularies. For example,
there are several enzyme inhibitors that refer to ExPASy enzyme class codes,
but there are no explicit links in ChEBI to support this.

See also: https://github.com/ebi-chebi/ChEBI/issues/3429

Repository Structure
--------------------
- ``resources/``: Manually curated resources
  - `resources/xrefs.tsv <https://github.com/cthoyt/chebi-relations/blob/master/resources/xrefs.tsv>`_:
    Structured annotations of ChEBI entries describing the relationship to an external vocabulary. For example,
    the term CETP inhibitor (CHEBI:49205) describes the ``inhibition`` relationship to CETP (hgnc:1869).
  - `resources/reclassification.tsv <https://github.com/cthoyt/chebi-relations/blob/master/resources/reclassification.tsv>`_:
    Suggestions on which entities inherit from the wrong type of ChEBI inhibitor role.
- ``export/``: Final results exported by combining the ChEBI OBO information
  with the manually curated resources in ``resources/``
- ``curate_chebi_relations.py``: The python script that suggests new curation
  and generates the derived resources in ``export/``

License
-------
- Code in this repository is under the MIT License.
- Data in this repository is under CC BY 4.0.

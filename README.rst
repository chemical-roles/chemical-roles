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
- ``export/``: Final results exported by combining the ChEBI OBO information
  with the manually curated resources in ``resources/``
- ``curate_chebi_relations.py``: The python script that suggests new curation
  and generates the derived resources in ``export/``

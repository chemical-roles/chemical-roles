Chemical Relations |zenodo|
===========================
This repository is meant to house links between ChEBI terms describing various
types of inhibitors and their targets in other vocabularies. For example,
there are several enzyme inhibitors that refer to ExPASy enzyme class codes,
but there are no explicit links in ChEBI to support this.

See also: https://github.com/ebi-chebi/ChEBI/issues/3429

License
-------
- Code in this repository is under the MIT License.
- Data in this repository is under CC0.

Repository Structure
--------------------
- ``resources/``: Manually curated resources

  - `resources/xrefs.tsv </resources/xrefs.tsv>`_:
    Structured annotations of ChEBI entries describing the relationship to an external vocabulary. For example,
    the term CETP inhibitor (CHEBI:49205) describes the ``inhibition`` relationship to CETP (hgnc:1869).
  - `resources/reclassification.tsv </resources/reclassification.tsv>`_:
    Suggestions on which entities inherit from the wrong type of ChEBI inhibitor role.
- ``export/``: Final results exported by combining the ChEBI OBO information
  with the manually curated resources in ``resources/``

  - `export/relations.tsv </export/relations.tsv>`_:
    A combination of the manually curated resources and reasoning (TODO) over FamPlex, InterPro, GO, HGNC
    Gene Familes, and ExPASy to generate new relationships between instances of the roles curated here.
- ``curate_chebi_relations.py``: The python script that suggests new curation using `Gilda <https://github.com/indralab/gilda>`_
  and generates the derived resources in ``export/``

Contributing
------------
1. Fork the repository and clone it with:

.. code-block:: sh

    git clone https://github.com/<your github username>/chemical-relations.git
    cd chemical-relations
    pip install -r requirements.txt

2. Run either ``python curate_chebi_relations.py`` or ``tox -e mesh`` to get a list of
   possible curatable rows. For example:

.. code-block::

   mesh	D064449	Sequestering Agents	agent	?	?	?	?	?
   mesh	D015842	Serine Proteinase Inhibitors	inhibitor	?	?	?	?	?
   mesh	D058825	Serotonin 5-HT1 Receptor Agonists	agonist	?	?	?	?	?
   mesh	D058829	Serotonin 5-HT1 Receptor Antagonists	antagonist	?	?	?	?	?

3. Copy/paste the rows you want to curate into the end of ``resources/xrefs.tsv``. Fill out the last 5
   columns, which are:
   1. modulation (the relationship)
   2. type (the type of the target)
   3. target_db
   4. target_id
   5. target_name

4. If you find rows that should be blacklisted in MeSH, add them to the blacklist in ``curate_mesh_relations.py``. If
   You find irrelevant rows while curating ChEBI, add them to ``resoures/irrelevant_roles.tsv``. This will be
   reorganized to be more consistent for MeSH soon.
5. Run ``tox -e lint`` to make sure you didn't make a mess
6. Make a PR back to https://github.com/chemical-roles/chemical-roles.git

Data Sources
------------
Automatically Retrieved Data Sources
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
- ChEBI (https://www.ebi.ac.uk/chebi)
- HUGO Gene Nomenclature Committee at the European Bioinformatics Institute. See https://www.genenames.org/ for more information
- UniProt (https://www.uniprot.org)
- FamPlex (https://github.com/sorgerlab/famplex)
- ExPASy (https://www.expasy.org/)

Manually Referenced Data Sources
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
- Gene Ontology (http://geneontology.org)
- EFO (https://www.ebi.ac.uk/efo)
- NCIT (https://ncit.nci.nih.gov)

Summary
-------
There are 829 curated roles as of export on Sun Jun  7 13:37:55 2020

===============  =======
Modulation         Count
===============  =======
activator             15
agonist               87
antagonist            94
inhibitor            625
inverse agonist        2
modulator              6
===============  =======

====================  =======
Target Entity Type      Count
====================  =======
biological process         42
chemical                   14
molecular function         39
organism                    9
phenotype                  15
protein                   124
protein complex            14
protein family            572
====================  =======

=================  =======
Target Database      Count
=================  =======
chebi                   13
ec-code                499
efo                     11
fplx                    39
go                      88
hgnc                   109
hgnc.genefamily          9
hp                       6
interpro                 1
mesh                    27
ncbitaxon                9
ncit                     1
pr                      15
uniprot                  2
=================  =======

Axioms
------
One of the main goals of this repository is to provide a framework for reasoning over roles (or families)
in ChEBI that don't have enough metadata.

Chemical-Physical Entity and Chemical-Process
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
This repo annotates relationships between chemical role ``X`` and protein ``Z`` such that:

- X isA chebi:agonist
- Y hasRole X
- Y agonistOf Z

It can also be reasoned over the hierarchy of children of chemical ``Y`` such that:

- X isA chebi:agonist
- Y hasRole X
- Y agonistOf Z
- y isA* Y
- y agonistOf Z

This repo annotates relationships between chemical role ``X`` and protein family ``Z`` such that:

- X isA chebi:agonist
- Y hasRole X
- Y agonistOf Z
- z isA* Z
- Y agonistOf z

And a combination of both the hierarchy of children of chemical ``Y`` and the children of protein family ``Z`` such
that:

- X isA chebi:agonist
- Y hasRole X
- Y agonistOf Z
- y isA* Y
- z isA* Z
- y agonistOf z

In general, this repository maps many ChEBI roles ``R`` to relationships ``r`` such that:

- X isA R
- Y hasRole X
- R roleHasRelation r
- Y r Z
- y isA* Y
- z isA* Z
- y r z

Chemical and Activity
~~~~~~~~~~~~~~~~~~~~~
This repo annotates relationships between chemical role ``X`` and activity ``A`` such

- X hasRole R
- R roleHasActivityRelation ar
- X ar A

When this is true, we can further infer the action of chemical role ``X`` on protein ``P``
that has activity ``A``:

- R roleHasEntityRelation er
- P isA protein
- P hasActivity A
- X er P

.. |zenodo| image:: https://zenodo.org/badge/199155107.svg
   :target: https://zenodo.org/badge/latestdoi/199155107

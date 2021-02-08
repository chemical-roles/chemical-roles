.. image:: art/CRoG-logotype-1024.png
   :alt: CRoG Logotype

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

Citation
--------
`Extension of Roles in the ChEBI Ontology <https://doi.org/10.26434/chemrxiv.12591221>`_.
Hoyt, C. T., Mungall, C., Vasilevsky, N., Domingo-Fernández, D., Healy, M., and Colluru, V. (2020).
*chemRxiv*, 12591221.

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

2. Run either ``tox -e chebi`` to curate ChEBI terms that Gilda was able to look up,
   ``tox -e chebi-ungrounded`` to curate all ChEBI terms, or ``tox -e mesh`` to curate MeSH.
   For example, the output ``tox -e mesh`` will include something like:

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
There are 836 curated roles as of export on Mon Feb  8 01:13:55 2021

===============  =======
Modulation         Count
===============  =======
Activator             15
Agonist               87
Antagonist            94
Inhibitor            629
Inverse Agonist        2
Modulator              9
===============  =======

====================  =======
Target Entity Type      Count
====================  =======
Biological Process         42
Chemical                   14
Molecular Function         39
Organism                    9
Phenotype                  15
Protein                   126
Protein Complex            15
Protein Family            576
====================  =======

=================  =======
Target Database      Count
=================  =======
chebi                   13
eccode                 499
efo                     11
fplx                    43
go                      88
hgnc                   110
hgnc.genefamily          9
hp                       6
interpro                 1
mesh                    29
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

.. raw:: html

   <p align="center">
    <img src="art/CRoG-king-100.png" alt="Keep Calm and CRoG On">
   </p>

.. |zenodo| image:: https://zenodo.org/badge/199155107.svg
   :target: https://zenodo.org/badge/latestdoi/199155107

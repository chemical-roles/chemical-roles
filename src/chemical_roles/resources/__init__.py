# -*- coding: utf-8 -*-

"""Chemical Roles resources."""

import os

HERE = os.path.abspath(os.path.dirname(__file__))

RECLASSIFICATION_PATH = os.path.join(HERE, 'reclassification.tsv')
UNCURATED_CHEBI_PATH = os.path.join(HERE, 'uncurated_chebi.tsv')
UNCURATED_MESH_PATH = os.path.join(HERE, 'uncurated_mesh.tsv')
XREFS_PATH = os.path.join(HERE, 'xrefs.tsv')
IRRELEVANT_ROLES_PATH = os.path.join(HERE, 'irrelevant_roles.tsv')
BLACKLIST_ROLES_PATH = os.path.join(HERE, 'blacklist.tsv')

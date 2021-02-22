# -*- coding: utf-8 -*-

"""Constants for Chemical Roles."""

import os

HERE = os.path.abspath(os.path.dirname(__file__))
ROOT = os.path.abspath(os.path.join(HERE, os.pardir, os.pardir))
DOCS = os.path.join(ROOT, 'docs')
DATA = os.path.join(DOCS, '_data')
IMG = os.path.join(DOCS, 'img')

RELATIONS_OUTPUT_PATH = os.path.join(DATA, 'relations.tsv')
RELATIONS_SLIM_OUTPUT_PATH = os.path.join(DATA, 'relations_slim.tsv')
EXPORT_BEL_PATH = os.path.join(DATA, 'export.bel.nodelink.json.gz')

# -*- coding: utf-8 -*-

"""Constants for Chemical Roles."""

import os

HERE = os.path.abspath(os.path.dirname(__file__))
EXPORT_DIRECTORY = os.path.join(HERE, 'export')
RELATIONS_OUTPUT_PATH = os.path.join(EXPORT_DIRECTORY, 'relations.tsv')
RELATIONS_SLIM_OUTPUT_PATH = os.path.join(EXPORT_DIRECTORY, 'relations_slim.tsv')
EXPORT_BEL_PATH = os.path.join(EXPORT_DIRECTORY, 'export.bel.nodelink.json.gz')

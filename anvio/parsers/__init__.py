#!/usr/bin/env python
# -*- coding: utf-8

import anvio
import anvio.terminal as terminal

from anvio.errors import ConfigError

from anvio.parsers.kraken_hll import KrakenHLL
from anvio.parsers.defaultmatrix import DefaultMatrix
from anvio.parsers.centrifuge import Centrifuge
from anvio.parsers.kaiju import Kaiju
from anvio.parsers.hmmscan import HMMScan
from anvio.parsers.concoct import CONCOCT
from anvio.parsers.interproscan import InterProScan


__author__ = "Developers of anvi'o (see AUTHORS.txt)"
__copyright__ = "Copyleft 2015-2018, the Meren Lab (http://merenlab.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "A. Murat Eren"
__email__ = "a.murat.eren@gmail.com"


parser_modules = {}
parser_modules['taxonomy_genes']  = {"default_matrix": DefaultMatrix, "centrifuge": Centrifuge, 'kaiju': Kaiju}
parser_modules['taxonomy_layers'] = {"kraken_hll": KrakenHLL}
parser_modules['functions']       = {"interproscan": InterProScan}
parser_modules['search']          = {"hmmscan": HMMScan}
parser_modules['collections']     = {"concoct": CONCOCT}

run = terminal.Run()

def get_parser_names(module):
    parser_module = get_parser_module(module)

    return list(parser_module.keys())


def get_parser_module(module):
    if module not in parser_modules:
        raise ConfigError("Anvi'o parser modules do not recognize any module called '%s'. But it has\
                           these in case if your honour would change their minds: '%s'." % \
                                (module, ', '.join(list(parser_modules.keys()))))

    return parser_modules[module]


def get_parser_obj(module, parser):
    parser_module = get_parser_module(module)

    if parser not in parser_module:
        raise ConfigError("Parser modules speaking: the parser module '%s' does exist (yay),\
                           but there is no parser '%s' in it (boo). It has these instad: '%s'."\
                                % (module, parser, ', '.join(list(parser_module.keys()))))

    return parser_modules[module][parser]

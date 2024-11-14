#!/usr/bin/env python
# -*- coding: utf-8

import anvio
import anvio.terminal as terminal

from anvio.errors import ConfigError

from anvio.parsers.krakenuniq import KrakenUniq
from anvio.parsers.defaultmatrix import DefaultMatrix
from anvio.parsers.centrifuge import Centrifuge
from anvio.parsers.kaiju import Kaiju
from anvio.parsers.hmmer import HMMERTableOutput, HMMERStandardOutput
from anvio.parsers.concoct import CONCOCT
from anvio.parsers.interproscan import InterProScan
from anvio.parsers.agnostos import AGNOSTOS


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "A. Murat Eren"
__email__ = "a.murat.eren@gmail.com"


parser_modules = {}
parser_modules["taxonomy_genes"] = {
    "default_matrix": DefaultMatrix,
    "centrifuge": Centrifuge,
    "kaiju": Kaiju,
}
parser_modules["taxonomy_layers"] = {"krakenuniq": KrakenUniq}
parser_modules["functions"] = {"interproscan": InterProScan, "AGNOSTOS": AGNOSTOS}
parser_modules["search"] = {
    "hmmer_table_output": HMMERTableOutput,
    "hmmer_std_output": HMMERStandardOutput,
}
parser_modules["collections"] = {"concoct": CONCOCT}

run = terminal.Run()


def get_parser_names(module):
    parser_module = get_parser_module(module)

    return list(parser_module.keys())


def get_parser_module(module):
    if module not in parser_modules:
        raise ConfigError(
            "Anvi'o parser modules do not recognize any module called '%s'. But it has "
            "these in case your honour would change their minds: '%s'."
            % (module, ", ".join(list(parser_modules.keys())))
        )

    return parser_modules[module]


def get_parser_obj(module, parser):
    parser_module = get_parser_module(module)

    if parser not in parser_module:
        raise ConfigError(
            "Parser modules speaking: the parser module '%s' does exist (yay),\
                           but there is no parser '%s' in it (boo). It has these instad: '%s'."
            % (module, parser, ", ".join(list(parser_module.keys())))
        )

    return parser_modules[module][parser]

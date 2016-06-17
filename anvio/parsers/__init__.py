#!/usr/bin/env python
# -*- coding: utf-8

# Copyright (C) 2015, anvio Developers
#
# This program is free software; you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free
# Software Foundation; either version 2 of the License, or (at your option)
# any later version.
#
# Please read the COPYING file.

from anvio.parsers.defaultmatrix import DefaultMatrix
from anvio.parsers.centrifuge import Centrifuge
from anvio.parsers.hmmscan import HMMScan
from anvio.parsers.concoct import CONCOCT
from anvio.parsers.interproscan import InterProScan

parser_modules = {}
parser_modules['taxonomy'] = {"default_matrix": DefaultMatrix,
                              "centrifuge": Centrifuge}
parser_modules['functions'] = {"interproscan": InterProScan}
parser_modules['search'] = {"hmmscan": HMMScan}
parser_modules['collections'] = {"concoct": CONCOCT}

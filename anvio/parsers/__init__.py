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
from anvio.parsers.myrastgui import MyRastGUI
from anvio.parsers.myrastcmdline import MyRastCMDLine 
from anvio.parsers.myrastcmdline_do_not_use import MyRastCMDLine_DO_NOT_USE
from anvio.parsers.hmmscan import HMMScan
from anvio.parsers.concoct import CONCOCT

parser_modules = {}
parser_modules['genes']       = {"default_matrix": DefaultMatrix,
                                 "myrast_gui": MyRastGUI,
                                 "myrast_cmdline": MyRastCMDLine,
                                 "myrast_cmdline_dont_use": MyRastCMDLine_DO_NOT_USE}
parser_modules['search']      = {"hmmscan": HMMScan}
parser_modules['collections'] = {"concoct": CONCOCT}

#!/usr/bin/env python
# -*- coding: utf-8

# Copyright (C) 2015, PaPi Developers
#
# This program is free software; you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free
# Software Foundation; either version 2 of the License, or (at your option)
# any later version.
#
# Please read the COPYING file.

from PaPi.parsers.defaultmatrix import DefaultMatrix
from PaPi.parsers.myrastgui import MyRastGUI
from PaPi.parsers.myrastcmdline import MyRastCMDLine 
from PaPi.parsers.myrastcmdline_do_not_use import MyRastCMDLine_DO_NOT_USE
from PaPi.parsers.hmmscan import HMMScan
from PaPi.parsers.concoct import CONCOCT

parser_modules = {}
parser_modules['annotation'] = {"default_matrix": DefaultMatrix,
                                "myrast_gui": MyRastGUI,
                                "myrast_cmdline": MyRastCMDLine,
                                "myrast_cmdline_dont_use": MyRastCMDLine_DO_NOT_USE}
parser_modules['search']     = {'hmmscan': HMMScan}
parser_modules['binning']    = {'concoct': CONCOCT}
# -*- coding: utf-8
# pylint: disable=line-too-long
"""A module to find Diversity Generating Retroelements"""

import subprocess
import argparse
import json
import re
import xml.etree.ElementTree as ET

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

import anvio.terminal as terminal

__author__ = "Developers of anvi'o (see AUTHORS.txt)"
__copyright__ = "Copyleft 2015-2024, the Meren Lab (http://merenlab.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "Katy Lambert-Slosarska"
__email__ = "klambertslosarska@gmail.com"
__status__ = "Development"

class DGR_Finder:
    def __init__(self, args, run=terminal.Run(), progress=terminal.Progress()):
        self.args = args
        self.run = run
        self.progress = progress
        

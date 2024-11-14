"""
Positional arguments:
1. INPUT - file path to FASTA file
2. OUTPUT - file path of output PIR file
"""

import sys

FASTA = sys.argv[1]
PIR = sys.argv[2]

from modeller import *

e = environ()
a = alignment(e, file=FASTA, alignment_format="FASTA")
a.write(file=PIR, alignment_format="PIR")

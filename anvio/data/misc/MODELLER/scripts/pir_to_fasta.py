"""
Positional arguments:
1. INPUT - file path to FASTA file
2. OUTPUT - file path of output PIR file
"""

import sys

PIR = sys.argv[1]
FASTA = sys.argv[2]

from modeller import *

e = environ()
a = alignment(e, file=PIR, alignment_format="PIR")
a.write(file=FASTA, alignment_format="FASTA")

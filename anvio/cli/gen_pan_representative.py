#!/usr/bin/env python
"""Generates a pangenome-supplemented representative genome from a pangenome analysis"""

import sys

import anvio
import anvio.terminal as terminal

from anvio.panrep import PanRepresenter
from anvio.argparse import ArgumentParser
from anvio.errors import ConfigError, FilesNPathsError


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ["Med1Bel", "sarilog"]
__requires__ = ["pan-db", "genomes-storage-db", "external-genomes"]
__provides__ = ["contigs-db"]
__description__ = "Generates a pangenome-supplemented representative genome from a pangenome"


run = terminal.Run()
progress = terminal.Progress()

def main():
    args = get_args()

    try:
        pan_representative = PanRepresenter(args, run, progress)
        pan_representative.process()
    except ConfigError as e:
        print(e)
        sys.exit(-1)
    except FilesNPathsError as e:
        print(e)
        sys.exit(-2)

def get_args():
    parser = ArgumentParser(description=__description__)

    groupA = parser.add_argument_group("INPUT FILES", "Input files from the pangenome analysis.")
    groupA.add_argument(*anvio.A("external-genomes"), **anvio.K("external-genomes", {'required': True}))
    groupA.add_argument(*anvio.A("genomes-storage"), **anvio.K("genomes-storage", {'required': True}))
    groupA.add_argument(*anvio.A("pan-db"), **anvio.K("pan-db"))

    groupB = parser.add_argument_group("SUPPLEMENTARY-CONTIG-OPTIONS", "All these options are related to the supplementary contig that will be added to the representative genome.")
    groupB.add_argument("--gap-size", metavar='INT', default=20)
    groupB.add_argument("--alpha", metavar='FLOAT', default=0.8)
    groupB.add_argument("--representative", metavar='GENOME-NAME')
    groupB.add_argument("--max-num-contigs", metavar='INT', default=99999)
    groupB.add_argument("--keep-synteny", action="store_true")
    groupB.add_argument("--keep-promoter", action="store_true")

    groupC = parser.add_argument_group("OUTPUT", "All these options are related to the contigs_db that will be generated.")
    groupC.add_argument(*anvio.A("output-file"), required=True, **anvio.K("output-file"))
    groupC.add_argument(*anvio.A("project-name"), **anvio.K("project-name"))
    groupC.add_argument(*anvio.A("description"), **anvio.K("description"))
    groupC.add_argument(*anvio.A("kmer-size"), **anvio.K("kmer-size"))
    groupC.add_argument(*anvio.A("split-length"), **anvio.K("split-length"))
    groupC.add_argument(*anvio.A("skip-mindful-splitting"), **anvio.K("skip-mindful-splitting"))

    return parser.get_args(parser)


if __name__ == "__main__":
    main()

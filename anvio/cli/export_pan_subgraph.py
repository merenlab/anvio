#!/usr/bin/env python
# coding: utf-8

import sys

import anvio
import anvio.terminal as terminal
import anvio.pangenome_graph_preprocess as preprocess

from anvio.argparse import ArgumentParser
from anvio.errors import ConfigError, FilesNPathsError

__copyright__ = "Copyleft 2015-2026, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['meren', 'ahenoch']
__requires__ = ['external-genomes', 'pan-graph-db']
__provides__ = ['contigs-db']
__description__ = "Export genomic loci that is in between two nodes in a given pangenome graph."


@terminal.time_program
def main():
    args = get_args()

    try:
        prep = preprocess.external_genomes_preprocess(args)
        prep.process()
    except ConfigError as e:
        print(e)
        sys.exit(-1)
    except FilesNPathsError as e:
        print(e)
        sys.exit(-2)


def get_args():
    parser = ArgumentParser(description=__description__)
    groupA = parser.add_argument_group('INPUT FILES')
    groupA.add_argument(*anvio.A('pan-graph-db'), **anvio.K('pan-graph-db'))
    groupA.add_argument(*anvio.A('external-genomes'), **anvio.K('external-genomes', {'required': True}))

    groupB = parser.add_argument_group('SUBGRAPH COORDINATES', description="This is where you describe the loci that you wish anvi'o to export. "
                                       "Users will typically learn which nodes they are interested in by visualizing the pangenome graph in anvi'o "
                                       "using the program `anvi-display-pan-graph`.")
    groupB.add_argument(*anvio.A('graph-nodes'), **anvio.K('graph-nodes'))

    groupC = parser.add_argument_group('OUTPUT', description="Provide a directory path, and anvi'o will generate all your contigs databases for "
                                "exported loci in it.")
    groupC.add_argument(*anvio.A('output-dir'), **anvio.K('output-dir', {'required': True}))

    return parser.get_args(parser)


if __name__ == '__main__':
    main()

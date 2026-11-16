#!/usr/bin/env python
"""Get sequences for genes with hits in a functional annotation source"""

import sys

import anvio
import anvio.terminal as terminal
import anvio.functionops as functionops

from anvio.argparse import ArgumentParser
from anvio.errors import ConfigError, FilesNPathsError


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['meren']
__requires__ = ['contigs-db', 'functions']
__can_use__ = ['profile-db', 'collection', 'bin', 'external-genomes', 'internal-genomes']
__provides__ = ['genes-fasta']
__description__ = ("Get sequences for genes with hits in a single functional annotation source, optionally only for "
                   "those that match to specific accession ids, from one or more contigs-dbs, bins, or genomes. Similar "
                   "to anvi-get-sequences-for-hmm-hits, but without any gene concatenation / phylogenomics functionality.")


@terminal.time_program
def main():
    try:
        run_program()
    except ConfigError as e:
        print(e)
        sys.exit(-1)
    except FilesNPathsError as e:
        print(e)
        sys.exit(-2)


def run_program():
    args = get_args()
    run = terminal.Run()
    progress = terminal.Progress()

    if args.list_defline_variables:
        functionops.SequencesForFunctionHits.list_defline_variables(run=run)
        sys.exit()

    s = functionops.SequencesForFunctionHits(args, run=run, progress=progress)

    if args.list_annotation_sources:
        s.list_annotation_sources()
        sys.exit()

    s.process()


def get_args():
    parser = ArgumentParser(description=__description__)

    groupA = parser.add_argument_group('INPUT OPTION #1: CONTIGS DB', "Get sequences for all genes with hits in a single contigs-db.")
    groupA.add_argument(*anvio.A('contigs-db'), **anvio.K('contigs-db', {'required': False}))

    groupB = parser.add_argument_group('INPUT OPTION #2: CONTIGS DB + PROFILE DB', "Get sequences for genes in bins of a collection. "
                                       "Without a bin id or a bin ids file, every bin in the collection will be used.")
    groupB.add_argument(*anvio.A('profile-db'), **anvio.K('profile-db', {'required': False}))
    groupB.add_argument(*anvio.A('collection-name'), **anvio.K('collection-name'))
    groupB.add_argument(*anvio.A('bin-id'), **anvio.K('bin-id'))
    groupB.add_argument(*anvio.A('bin-ids-file'), **anvio.K('bin-ids-file'))

    groupC = parser.add_argument_group('INPUT OPTION #3: EXTERNAL / INTERNAL GENOMES', "Get sequences for genes in many genomes. "
                                       "You can use an external genomes file, an internal genomes file, or both.")
    groupC.add_argument(*anvio.A('external-genomes'), **anvio.K('external-genomes'))
    groupC.add_argument(*anvio.A('internal-genomes'), **anvio.K('internal-genomes'))

    groupD = parser.add_argument_group('FUNCTIONS', "Declare exactly one annotation source, and optionally limit the output "
                                       "to genes with hits to specific accession ids in it. If you are recovering sequences "
                                       "from multiple genomes you better have the same annotation source in all of them.")
    groupD.add_argument(*anvio.A('annotation-source'), **anvio.K('annotation-source', {'required': False}))
    groupD.add_argument(*anvio.A('accession-ids'), **anvio.K('accession-ids'))
    groupD.add_argument(*anvio.A('list-annotation-sources'), **anvio.K('list-annotation-sources'))

    groupE = parser.add_argument_group('OUTPUT', "A FASTA file with one entry per gene.")
    groupE.add_argument(*anvio.A('output-file'), **anvio.K('output-file', {'required': False}))
    groupE.add_argument(*anvio.A('get-aa-sequences'), **anvio.K('get-aa-sequences'))
    groupE.add_argument(*anvio.A('defline-format'), **anvio.K('defline-format'))
    groupE.add_argument(*anvio.A('list-defline-variables'), **anvio.K('list-defline-variables'))
    groupE.add_argument(*anvio.A('no-wrap'), **anvio.K('no-wrap'))

    return parser.get_args(parser)


if __name__ == '__main__':
    main()

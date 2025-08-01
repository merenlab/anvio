#!/usr/bin/env python
# -*- coding: utf-8

import sys

import anvio
import anvio.dbops as dbops
import anvio.terminal as terminal

from anvio.errors import ConfigError, FilesNPathsError


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['meren', 'ekiefl', 'ozcan']
__requires__ = ['contigs-fasta', 'external-gene-calls']
__provides__ = ['contigs-db']
__anvio_workflows__ = ['metagenomics']
__description__ = "Generate a new anvi'o contigs database"


@terminal.time_program
def main():
    args = get_args()
    run = terminal.Run()
    progress = terminal.Progress()

    try:
        a = dbops.ContigsDatabase(args.output_db_path, run, progress, quiet=False, skip_init = True)
        a.create(args)
    except ConfigError as e:
        print(e)
        sys.exit(-1)
    except FilesNPathsError as e:
        print(e)
        sys.exit(-2)


def get_args():
    from anvio.argparse import ArgumentParser
    parser = ArgumentParser(description=__description__)

    groupA = parser.add_argument_group('MANDATORY INPUTS', 'Things you really need to provide to be in business.')
    groupA.add_argument(*anvio.A('contigs-fasta'), **anvio.K('contigs-fasta'))
    groupA.add_argument(*anvio.A('project-name'), **anvio.K('project-name'))

    groupB = parser.add_argument_group('PERFORMANCE', 'You have multiple cores? WELL, USE THEM MAYBE?.')
    groupB.add_argument(*anvio.A('num-threads'), **anvio.K('num-threads'))

    groupC = parser.add_argument_group('OPTIONAL INPUTS', 'Things you may want to tweak.')
    groupC.add_argument(*anvio.A('output-db-path'), **anvio.K('output-db-path', {'default': 'CONTIGS.db'}))
    groupC.add_argument(*anvio.A('db-variant'), **anvio.K('db-variant'))
    groupC.add_argument(*anvio.A('description'), **anvio.K('description'))
    groupC.add_argument(*anvio.A('split-length'), **anvio.K('split-length'))
    groupC.add_argument(*anvio.A('skip-mindful-splitting'), **anvio.K('skip-mindful-splitting'))
    groupC.add_argument(*anvio.A('kmer-size'), **anvio.K('kmer-size'))

    groupD = parser.add_argument_group('GENES IN CONTIGS', 'Expert thingies.')
    groupD.add_argument(*anvio.A('skip-gene-calling'), **anvio.K('skip-gene-calling'))
    groupD.add_argument(*anvio.A('prodigal-single-mode'), **anvio.K('prodigal-single-mode'))
    groupD.add_argument(*anvio.A('prodigal-translation-table'), **anvio.K('prodigal-translation-table'))
    groupD.add_argument(*anvio.A('full-gene-calling-report'), **anvio.K('full-gene-calling-report'))
    groupD.add_argument(*anvio.A('external-gene-calls'), **anvio.K('external-gene-calls'))
    groupD.add_argument(*anvio.A('ignore-internal-stop-codons'), **anvio.K('ignore-internal-stop-codons'))
    groupD.add_argument(*anvio.A('skip-predict-frame'), **anvio.K('skip-predict-frame'))

    return parser.get_args(parser)


if __name__ == '__main__':
    main()

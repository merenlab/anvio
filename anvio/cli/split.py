#!/usr/bin/env python
# -*- coding: utf-8
"""Split an anvi'o profile into smaller pieces."""

import sys
from anvio.argparse import ArgumentParser

import anvio
import anvio.splitter as splitter
import anvio.constants as constants

from anvio.errors import ConfigError, FilesNPathsError


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['meren']
__requires__ = ['profile-db', 'contigs-db', 'genomes-storage-db', 'pan-db', 'collection',]
__provides__ = ['split-bins',]
__resources__ = [("Anvi-split in action in the pangenomics tutorial", "http://merenlab.org/2016/11/08/pangenomics-v2/#splitting-the-pangenome")]
__description__ = ("Split an anvi'o pan or profile database into smaller, self-contained projects. Black magic.")


def main():
    args = get_args()

    try:
        splitter.DBSplitter(args).get()(args).process()
    except ConfigError as e:
        print(e)
        sys.exit(-1)
    except FilesNPathsError as e:
        print(e)
        sys.exit(-1)


def get_args():
    parser = ArgumentParser(description=__description__)

    groupA = parser.add_argument_group('DATABASES', "You will either provide a PROFILE/CONTIGS or a PAN/GENOMES STORAGE pair here.")
    groupA.add_argument(*anvio.A('pan-or-profile-db'), **anvio.K('pan-or-profile-db'))
    groupA.add_argument(*anvio.A('contigs-db'), **anvio.K('contigs-db', {'required': False}))
    groupA.add_argument(*anvio.A('genomes-storage'), **anvio.K('genomes-storage', {'required': False}))


    groupB = parser.add_argument_group('PROFILE/CONTIGS OPTIONS', "Some options that are specific to this only.")
    groupB.add_argument(*anvio.A('skip-variability-tables'), **anvio.K('skip-variability-tables'))
    groupB.add_argument(*anvio.A('compress-auxiliary-data'), **anvio.K('compress-auxiliary-data'))

    groupC = parser.add_argument_group('COLLECTION', 'You should provide a valid collection name. If you do not provide\
                                                      bin names, the program will generate an output for each bin in your\
                                                      collection separately.')
    groupC.add_argument(*anvio.A('collection-name'), **anvio.K('collection-name'))
    groupC.add_argument(*anvio.A('bin-id'), **anvio.K('bin-id'))

    groupD = parser.add_argument_group('OUTPUT', 'Where do we want the resulting split profiles to be stored.')
    groupD.add_argument(*anvio.A('output-dir'), **anvio.K('output-dir'))

    groupE = parser.add_argument_group('EXTRAS', "Stuff that you rarely need, but you really really need when the time comes.\
                                                  Following parameters will aply to each of the resulting anvi'o profile that\
                                                  will be split from the mother anvi'o profile.")
    groupE.add_argument(*anvio.A('list-collections'), **anvio.K('list-collections'))
    groupE.add_argument(*anvio.A('skip-hierarchical-clustering'), **anvio.K('skip-hierarchical-clustering'))
    groupE.add_argument(*anvio.A('enforce-hierarchical-clustering'), **anvio.K('enforce-hierarchical-clustering'))
    groupE.add_argument(*anvio.A('distance'), **anvio.K('distance', {'default': None, 'help':
                      'The distance metric for the hierarchical clustering. If you do not use this flag,\
                       the default distance metric will be used for each clustering configuration\
                       which is "%s".' % constants.distance_metric_default}))
    groupE.add_argument(*anvio.A('linkage'), **anvio.K('linkage', {'default': None, 'help':
                      'The same story with the `--distance`, except, the system default for this one\
                       is %s.' % constants.linkage_method_default}))

    return parser.get_args(parser)


if __name__ == '__main__':
    main()

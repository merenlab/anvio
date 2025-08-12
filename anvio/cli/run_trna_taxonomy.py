#!/usr/bin/env python
# -*- coding: utf-8

import sys

import anvio
import anvio.terminal as terminal
import anvio.taxonomyops.trna as trnataxonomyops

from anvio.errors import ConfigError, FilesNPathsError


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['meren']
__requires__ = ['contigs-db', 'trna-taxonomy-db']
__provides__ = ['trna-taxonomy']
__resources__ = []
__description__ = ("The purpose of this program is to affiliate tRNA gene sequences in an anvi'o contigs database with "
                   "taxonomic names. A properly setup local tRNA taxonomy database is required for this program to perform properly. "
                   "After its successful run, `anvi-estimate-trna-taxonomy` will be useful to estimate taxonomy at genome-, collection-, or metagenome-level).")


@terminal.time_program
def main():
    args = get_args()

    try:
        t = trnataxonomyops.PopulateContigsDatabaseWithTRNATaxonomy(args)
        t.populate_contigs_database()
    except ConfigError as e:
        print(e)
        sys.exit(-1)
    except FilesNPathsError as e:
        print(e)
        sys.exit(-2)


def get_args():
    from anvio.argparse import ArgumentParser

    parser = ArgumentParser(description=__description__)

    groupA = parser.add_argument_group('INPUT DATABASE', "An anvi'o contigs databaes to search for and  store the taxonomic\
                                        affiliations of tRNA genes.")
    groupA.add_argument(*anvio.A('contigs-db'), **anvio.K('contigs-db', {'required': True}))

    groupA = parser.add_argument_group("ADVANCED STUFF")
    groupA.add_argument(*anvio.A('trna-taxonomy-data-dir'), **anvio.K('trna-taxonomy-data-dir'))
    groupA.add_argument(*anvio.A('min-percent-identity'), **anvio.K('min-percent-identity', {'help': "The defualt value for this is \
                        %(default).1f%%, and in an ideal world you sholdn't really change it. Lowering this value will probably give \
                        you too many hits from neighboring genomes, which may ruin your consensus taxonomy (imagine, at 90%% identity \
                        you may match to a single species, but at 70%% identity you may match to every species in a genus and your \
                        consensus assignment may be influenced by that). But once in a while you will have a genome that doesn't have any\
                        close match in GTDB, and you will be curious to find out what it could be. So, when you are getting no tRNA hits\
                        whatsoever, only then you may want to play with this value. In those cases you can run anvi-estimate-trna-taxonomy\
                        with a `--debug` flag to see what is really going on. We strongly advice you to do this only with single genomes,\
                        and never with metagenomes.", 'default': 90}))
    groupA.add_argument(*anvio.A('max-num-target-sequences'), **anvio.K('max-num-target-sequences', {'help': "This parameter is used to determine \
                        how many hits from the database that has a reasonable match to the query sequence should be taken into consideration \
                        to make a final decision about the consensus taxonomy for each individual transfer RNA gene sequence. The default \
                        is %(default)d, which has been quite reasonable in our tests, however, you may need to increase this number to get more \
                        accurate results for your own data. In cases where you think this is what you need, the best way to test the parameter \
                        space for `--max-num-target-sequences` is to run the program multiple times on the same database with `--debug` \
                        and compare results.", 'default': 100}))

    groupH = parser.add_argument_group("PERFORMANCE")
    groupH.add_argument(*anvio.A('num-parallel-processes'), **anvio.K('num-parallel-processes'))
    groupH.add_argument(*anvio.A('num-threads'), **anvio.K('num-threads'))
    groupH.add_argument(*anvio.A('write-buffer-size'), **anvio.K('write-buffer-size'))

    groupI = parser.add_argument_group("OUTPUT", "By default, this program does not generate an output and instead simply store taxonomy \
                                                  information into the contigs database. But if the user wants more, they get more.")
    groupI.add_argument(*anvio.A('all-hits-output-file'), **anvio.K('all-hits-output-file'))

    return parser.get_args(parser)

if __name__ == '__main__':
    main()

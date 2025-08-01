#!/usr/bin/env python
# -*- coding: utf-8

import sys

import anvio
import anvio.taxonomyops.trna as trnataxonomyops

from anvio.errors import ConfigError, FilesNPathsError

__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['meren']
__provides__ = ['trna-taxonomy-db']
__resources__ = []
__description__ = ("The purpose of this program is to setup necessary databases for tRNA genes collected from GTDB (https://gtdb.ecogenomic.org/), "
                   "genomes in your local anvi'o installation so taxonomy information for a given set of tRNA sequences can be identified "
                   "using `anvi-run-trna-taxonomy` and made sense of via `anvi-estimate-trna-taxonomy`)")


def main():
    args = get_args()

    try:
        s = trnataxonomyops.SetupLocalTRNATaxonomyData(args)
        s.setup()
    except ConfigError as e:
        print(e)
        sys.exit(-1)
    except FilesNPathsError as e:
        print(e)
        sys.exit(-2)


def get_args():
    from anvio.argparse import ArgumentParser

    parser = ArgumentParser(description=__description__)
    parser.add_argument(*anvio.A('num-threads'), **anvio.K('num-threads'))
    parser.add_argument(*anvio.A('trna-taxonomy-data-dir'), **anvio.K('trna-taxonomy-data-dir'))
    parser.add_argument(*anvio.A('reset'), **anvio.K('reset'))
    parser.add_argument(*anvio.A('redo-databases'), **anvio.K('redo-databases'))

    return parser.get_args(parser)


if __name__ == '__main__':
    main()

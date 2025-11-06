#!/usr/bin/env python
# -*- coding: utf-8

import sys

import anvio
import anvio.taxonomyops.scg as scgtaxonomyops

from anvio.errors import ConfigError, FilesNPathsError

__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['meren', 'qclayssen']
__provides__ = ['scgs-taxonomy-db']
__resources__ = [("Usage examples and warnings", "http://merenlab.org/scg-taxonomy")]
__description__ = ("The purpose of this program is to download necessary information from GTDB (https://gtdb.ecogenomic.org/), "
                   "and set it up in such a way that your anvi'o installation is able to assign taxonomy to single-copy core "
                   "genes using `anvi-run-scg-taxonomy` and estimate taxonomy for genomes or metagenomes using "
                   "`anvi-estimate-scg-taxonomy`)")


def main():
    args = get_args()

    try:
        s = scgtaxonomyops.SetupLocalSCGTaxonomyData(args)
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
    parser.add_argument(*anvio.A('scgs-taxonomy-data-dir'), **anvio.K('scgs-taxonomy-data-dir'))
    parser.add_argument(*anvio.A('gtdb-release'), **anvio.K('gtdb-release'))

    return parser.get_args(parser)


if __name__ == '__main__':
    main()

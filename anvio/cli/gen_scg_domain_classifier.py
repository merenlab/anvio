#!/usr/bin/env python
# -*- coding: utf-8
# pylint: disable=line-too-long

"""
    Generate a random forest classifier to predict SCG domain
"""

import sys

import anvio

from anvio.scgdomainclassifier import Train
from anvio.errors import ConfigError, FilesNPathsError


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['meren']
__description__ = "Train a classifier for SCG domain prediction"


def main():
    args = get_args()

    try:
        c = Train(args)
        c.train()
    except ConfigError as e:
        print(e)
        sys.exit(-1)
    except FilesNPathsError as e:
        print(e)
        sys.exit(-1)


def get_args():
    from anvio.argparse import ArgumentParser

    parser = ArgumentParser(description=__description__)

    parser.add_argument('--genomes-dir', help = "This should be a directory that contains a directory per domain for \
                                                 single-copy core gene collections a given version of anvi'o knows about.\
                                                 For instance, if there are collections for archaea, bacteria, and\
                                                 eukarya, then this directory should contain subdirectories with these\
                                                 names. Contents of which should be contigs databases that belong to\
                                                 those domains. These genomes will be used to generate the classifier.")
    parser.add_argument('-o', '--output', default = "domain.classifier", help = "Output file name for the classifier.",
                        metavar='PATH')

    return parser.get_args(parser)


if __name__ == '__main__':
    main()

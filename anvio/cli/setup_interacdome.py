#!/usr/bin/env python
# -*- coding: utf-8

import sys

import anvio
import anvio.interacdome as interacdome

from anvio.errors import ConfigError, FilesNPathsError


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['ekiefl']
__provides__ = ['interacdome-data']
__description__ = "Setup InteracDome data"
__resources__ = [("The setup step in the InteracDome technical blogpost", "http://merenlab.org/2020/07/22/interacdome/#anvi-setup-interacdome")]


def main():
    args = get_args()

    try:
        setup = interacdome.InteracDomeSetup(args)
        setup.setup()
    except ConfigError as e:
        print(e)
        sys.exit(-1)
    except FilesNPathsError as e:
        print(e)
        sys.exit(-1)


def get_args():
    from anvio.argparse import ArgumentParser

    parser = ArgumentParser(description=__description__)
    parser.add_argument(*anvio.A('interacdome-data-dir'), **anvio.K('interacdome-data-dir'))
    parser.add_argument(*anvio.A('reset'), **anvio.K('reset'))

    args = parser.get_args(parser)


if __name__ == '__main__':
    main()

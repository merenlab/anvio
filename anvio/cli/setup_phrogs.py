#!/usr/bin/env python

import sys

import anvio
import anvio.phrogs as phrogs

from anvio.errors import ConfigError, FilesNPathsError


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['copilot']
__provides__ = ['phrogs-data']
__description__ = "Download and setup PHROGs data"


def main():
    args = get_args()

    try:
        setup = phrogs.PHROGsSetup(args)
        setup.download()
    except ConfigError as e:
        print(e)
        sys.exit(-1)
    except FilesNPathsError as e:
        print(e)
        sys.exit(-1)


def get_args():
    from anvio.argparse import ArgumentParser

    parser = ArgumentParser(description=__description__)
    parser.add_argument('--phrogs-data-dir', default=None, type=str, help="The directory for PHROGs data to be stored. If you leave it "
                        "empty, anvi'o uses the default location. If you choose a custom path, you should provide the same path "
                        "later when running `anvi-run-phrogs`.")
    parser.add_argument('--reset', default=False, action="store_true", help="Delete existing PHROGs setup and download everything again.")
    parser.add_argument('--phrogs-version', default='v4', help="PHROGs annotation version suffix to download for annotations. "
                        "For instance, use 'v4' to download phrog_annot_v4.tsv.")

    return parser.get_args(parser)


if __name__ == '__main__':
    main()

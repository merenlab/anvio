#!/usr/bin/env python
# -*- coding: utf-8

import sys

import anvio
import anvio.cazymes as cazymes

from anvio.errors import ConfigError, FilesNPathsError


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['mschecht']
__provides__ = ['cazyme-data']
__description__ = "Download and setup Pfam data from the EBI"


def main():
    args = get_args()

    try:
        setup = cazymes.CAZymeSetup(args)
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
    parser.add_argument('--cazyme-data-dir', default=None, type=str, help="The directory for CAZyme data to be stored. If you leave it\
                        as is without specifying anything, the default destination for the data directory will be used to set things\
                        up. The advantage of it is that everyone will be using a single data directory, but then you may need\
                        superuser privileges to do it. Using this parameter you can choose the location of the data directory somewhere\
                        you like. However, when it is time to run CAZyme, you will need to remember that path and provide it to the program.")
    parser.add_argument('--reset', default=False, action="store_true", help="This program by default attempts to use previously\
                        downloaded files in your Pfam data directory if there are any. If something is wrong for some reason you\
                        can use this to tell anvi'o to remove everything, and start over.")
    parser.add_argument('--cazyme-version', default=None, help="By default, the most current version available (V11) will be downloaded.\
                        If you have specific tastes for a different version, you can provide it here. For example, `V10`. Here are\
                        all possible versions: https://bcb.unl.edu/dbCAN2/download/Databases/")

    return parser.get_args(parser)


if __name__ == '__main__':
    main()

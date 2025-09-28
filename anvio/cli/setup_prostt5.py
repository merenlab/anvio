#!/usr/bin/env python
# -*- coding: utf-8

import sys
import anvio

from anvio.drivers.foldseek import Prostt5SetupWeight
from anvio.errors import ConfigError, FilesNPathsError


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "Metehan Sever"
__email__ = "metehaansever@gmail.com"
__provides__ = ['prostt5-weights']
__description__ = "Setup PROSTT5 weights"


def main():
    args = get_args()

    try:
        setup = Prostt5SetupWeight(args)
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
    parser.add_argument(*anvio.A('prostt5-data-dir'), **anvio.K('prostt5-data-dir'))
    parser.add_argument(*anvio.A('reset'), **anvio.K('reset'))

    return parser.get_args(parser)

if __name__ == '__main__':
    main()

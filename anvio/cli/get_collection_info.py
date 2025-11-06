#!/usr/bin/env python
# -*- coding: utf-8

import sys

from anvio.errors import ConfigError

__status__ = "Deprecated"
__description__ = "Replaced by 'anvi-estimate-genome-completeness'"

def main():
    try:
        raise ConfigError("This program is now known as 'anvi-estimate-genome-completeness'. At least "
                          "we made sure you can re-run the same command line above by only replacing "
                          "the program name.")
    except ConfigError as e:
        print(e)
        sys.exit(1)


if __name__ == '__main__':
    main()

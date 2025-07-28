#!/usr/bin/env python
# -*- coding: utf-8

import sys

from anvio.errors import ConfigError


__description__ = ("This program has been superseded by the beefier `anvi-compute-genome-similarity`")


def main():
    raise ConfigError(__description__)


if __name__ == '__main__':
    try:
        main()
    except ConfigError as e:
        print(e)
        sys.exit(-1)

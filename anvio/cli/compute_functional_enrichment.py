#!/usr/bin/env python
# -*- coding: utf-8
import sys

from anvio.errors import ConfigError
from anvio.terminal import Run

__description__ = ("A program that is no more. Read on")


def main():
    try:
        run = Run()
        run.warning("This program is broken into multiple programs. Our VERY FAVORITE one is "
                    "the metabolic enrichment one, but you feel free to choose the one from "
                    "the list of minions below based on what you wish to do:",
                    header="THIS PROGRAM IS NO MORE", lc="yellow")

        run.info_single('anvi-compute-metabolic-enrichment')
        run.info_single('anvi-compute-functional-enrichment-in-pan')
        run.info_single('anvi-compute-functional-enrichment-across-genomes')

        raise ConfigError("You've been served (but in a good way).")
    except ConfigError as e:
        print(e)
        sys.exit(-1)


if __name__ == '__main__':
    main()

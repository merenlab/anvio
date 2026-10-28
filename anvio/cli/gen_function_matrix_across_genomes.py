#!/usr/bin/env python

import sys

from anvio.errors import ConfigError

__status__ = "Deprecated"
__description__ = "Replaced by 'anvi-gen-function-matrix'"

def main():
    try:
        raise ConfigError("This program is now known as 'anvi-gen-function-matrix'. We renamed it "
                          "because then additional features made applicable to metagenomes and "
                          "other kinds of data in contigs-db files. You can re-run the same "
                          "commands by simply replacing the program name.")
    except ConfigError as e:
        print(e)
        sys.exit(1)


if __name__ == '__main__':
    main()

#!/usr/bin/env python

import sys

from anvio.errors import ConfigError

__status__ = "Deprecated"
__description__ = "Replaced by 'anvi-gen-hmm-hits-matrix'"

def main():
    try:
        raise ConfigError("This program is now known as 'anvi-gen-hmm-hits-matrix'. The previous "
                          "name was way too long for it. Everything else remained the same, so you "
                          "can simply replace the program name, and run the same command.")
    except ConfigError as e:
        print(e)
        sys.exit(1)


if __name__ == '__main__':
    main()

#!/usr/bin/env python
# -*- coding: utf-8

import sys
import argparse

import anvio.db as db
import anvio.tables as t
import anvio.utils as utils
import anvio.terminal as terminal

from anvio.errors import ConfigError


run = terminal.Run()
progress = terminal.Progress()


def migrate(db_path):
    pass

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='A simple script to upgrade KEGG Modules database from version 1 to version 2')
    parser.add_argument('modules_db', metavar = 'MODULES_DB', help = 'KEGG Modules database')
    args, unknown = parser.parse_known_args()

    try:
        migrate(args.modules_db)
    except ConfigError as e:
        print(e)
        sys.exit(-1)

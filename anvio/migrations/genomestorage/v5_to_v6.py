#!/usr/bin/env python
# -*- coding: utf-8

import sys
import argparse
import anvio.db as db
import anvio.terminal as terminal

from anvio.errors import ConfigError

current_version = "5"
next_version = "6"

run = terminal.Run()
progress = terminal.Progress()


def migrate(db_path):
    if db_path is None:
        raise ConfigError("No database path is given.")

    genomes_db = db.DB(db_path, current_version)
    gene_function_sources = ",".join(
        genomes_db.get_single_column_from_table(
            "gene_function_calls", "source", unique=True
        )
    )
    genomes_db.set_meta_value("gene_function_sources", gene_function_sources)
    genomes_db.set_version(next_version)
    genomes_db.disconnect()

    run.info_single(
        "Your genomes storage is now at version %s. " % (next_version),
        nl_after=1,
        nl_before=1,
        mc="green",
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="A simple script to upgrade genomes storage from version %s to version %s"
        % (current_version, next_version)
    )
    parser.add_argument(
        "genomes_storage",
        metavar="GENOMES_STORAGE",
        help="An anvi'o genomes storage of version %s" % current_version,
    )
    args, unknown = parser.parse_known_args()

    try:
        migrate(args.genomes_storage)
    except ConfigError as e:
        print(e)
        sys.exit(-1)

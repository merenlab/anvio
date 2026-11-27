#!/usr/bin/env python

import sys
import argparse

import anvio.db as db
import anvio.utils as utils
import anvio.terminal as terminal

from anvio.errors import ConfigError

run = terminal.Run()
progress = terminal.Progress()

current_version, next_version = [x[1:] for x in __name__.split('_to_')]

# these are the names as they were at the time this script was written. they are hard-coded
# here on purpose, so this script does not depend on table definitions in the codebase.
inseq_table_name = 'gene_level_inseq_stats'
standard_table_name = 'gene_level_coverage_stats'
gene_level_stats_parameter_keys = ['min_cov_for_detection', 'outliers_threshold', 'zeros_are_outliers']


def migrate(db_path):
    if db_path is None:
        raise ConfigError("No database path is given.")

    # make sure someone is not being funny
    utils.is_genes_db(db_path)

    # make sure the version is accurate
    genes_db = db.DB(db_path, None, ignore_version = True)
    if str(genes_db.get_version()) != current_version:
        genes_db.disconnect()
        raise ConfigError(f"Version of this genes database is not {current_version} (hence, this script cannot really do anything).")

    progress.new("Migrating")
    progress.update("Removing INSEQ/Tn-Seq remnants")

    # INSEQ/Tn-Seq support is no longer a part of anvi'o, so the table that used to hold
    # INSEQ stats goes away regardless of whether it contains anything
    genes_db._exec(f'''DROP TABLE IF EXISTS {inseq_table_name}''')

    # if this database was generated in INSEQ mode, its gene-level stats are not usable anymore.
    # here we reset it to the state of a freshly created genes database, so anvi'o will compute
    # standard gene-level stats the next time they are needed
    mode = genes_db.get_meta_value('mode', try_as_type_int=False, return_none_if_not_in_table=True)
    was_inseq = mode is not None and str(mode).upper() == 'INSEQ'

    if was_inseq:
        if standard_table_name in genes_db.get_table_names():
            genes_db._exec(f'''DELETE FROM {standard_table_name}''')

        for key in gene_level_stats_parameter_keys:
            genes_db.remove_meta_key_value_pair(key)

        genes_db.update_meta_value('gene_level_coverages_stored', False)
        genes_db.update_meta_value('mode', 'STANDARD')

    # set the version
    progress.update("Setting the version")
    genes_db.remove_meta_key_value_pair('version')
    genes_db.set_version(next_version)

    # bye
    genes_db.disconnect()
    progress.end()

    if was_inseq:
        run.warning("This genes database was generated in INSEQ/Tn-Seq mode. Since anvi'o no longer supports "
                    "INSEQ/Tn-Seq, the INSEQ data stored in it has been removed, and the database is now set to "
                    "the 'STANDARD' mode. Anvi'o will recompute standard gene-level coverage stats for this bin "
                    "the next time they are needed (i.e., when you run `anvi-interactive` with `--gene-mode`, "
                    "or `anvi-gen-gene-level-stats-databases`). Collections, states, and other data you may have "
                    "stored in this database are not affected.", header="INSEQ DATA REMOVED")

    run.info_single(f"Your genes db is now version {next_version}.", nl_after=1, nl_before=1, mc='green')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='A simple script to upgrade genes database from version %s to version %s' % (current_version, next_version))
    parser.add_argument('genes_db', metavar = 'GENES_DB', help = "An anvi'o genes database of version %s" % current_version)
    args, unknown = parser.parse_known_args()

    try:
        migrate(args.genes_db)
    except ConfigError as e:
        print(e)
        sys.exit(-1)

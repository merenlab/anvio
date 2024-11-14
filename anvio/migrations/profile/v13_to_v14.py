#!/usr/bin/env python
# -*- coding: utf-8

import sys
import argparse

import anvio.db as db
import anvio.utils as utils
import anvio.terminal as terminal

from anvio.errors import ConfigError


run = terminal.Run()
progress = terminal.Progress()


def migrate(db_path):
    if db_path is None:
        raise ConfigError("No profile database is given.")

    # make sure someone is not being funny
    utils.is_profile_db(db_path)

    # make sure the version is 5
    profile_db = db.DB(db_path, None, ignore_version=True)
    if str(profile_db.get_version()) != "13":
        raise ConfigError(
            "Version of this profile database is not 13 (hence, this script cannot really do anything)."
        )

    is_merged = profile_db.get_meta_value("merged")

    progress.new(
        "Trying to upgrade the %s profile database" % "merged"
        if is_merged
        else "single"
    )

    # serious stuff ###################################################################################################
    progress.update("working on the variable nts table ...")
    profile_db.cursor.execute(
        "ALTER TABLE variable_nucleotide_positions RENAME TO variable_nucleotide_positions_TEMP;"
    )
    profile_db.cursor.execute(
        "CREATE TABLE variable_nucleotide_positions (entry_id numeric, sample_id text, split_name text, pos numeric, pos_in_contig numeric, corresponding_gene_call numeric, in_partial_gene_call numeric, in_complete_gene_call numeric, base_pos_in_codon numeric, codon_order_in_gene numeric, coverage numeric, cov_outlier_in_split bool, cov_outlier_in_contig bool, departure_from_reference numeric, competing_nts text, reference text, A numeric, T numeric, C numeric, G numeric, N numeric);"
    )
    profile_db.cursor.execute(
        "INSERT INTO variable_nucleotide_positions(entry_id, sample_id, split_name, pos, pos_in_contig, corresponding_gene_call, in_partial_gene_call, in_complete_gene_call, base_pos_in_codon, codon_order_in_gene, coverage, cov_outlier_in_split, cov_outlier_in_contig, departure_from_reference, competing_nts, reference, A, T, C, G, N) SELECT entry_id, sample_id, split_name, pos, pos_in_contig, corresponding_gene_call, in_partial_gene_call, in_complete_gene_call, base_pos_in_codon, codon_order_in_gene, coverage, cov_outlier_in_split, cov_outlier_in_contig, departure_from_consensus, competing_nts, consensus, A, T, C, G, N FROM variable_nucleotide_positions_TEMP;"
    )
    profile_db.cursor.execute("DROP TABLE variable_nucleotide_positions_TEMP;")

    progress.update("working on the aa freqs table ...")
    profile_db.cursor.execute(
        "ALTER TABLE variable_amino_acid_frequencies RENAME TO variable_amino_acid_frequencies_TEMP;"
    )
    profile_db.cursor.execute(
        "CREATE TABLE variable_amino_acid_frequencies (entry_id numeric, sample_id text, corresponding_gene_call numeric, codon_order_in_gene numeric, reference text, departure_from_reference numeric, coverage numeric, Ala numeric, Arg numeric, Asn numeric, Asp numeric, Cys numeric, Gln numeric, Glu numeric, Gly numeric, His numeric, Ile numeric, Leu numeric, Lys numeric, Met numeric, Phe numeric, Pro numeric, STP numeric, Ser numeric, Thr numeric, Trp numeric, Tyr numeric, Val numeric);"
    )
    profile_db.cursor.execute(
        "INSERT INTO variable_amino_acid_frequencies(entry_id, sample_id, corresponding_gene_call, codon_order_in_gene, reference, departure_from_reference, coverage, Ala, Arg, Asn, Asp, Cys, Gln, Glu, Gly, His, Ile, Leu, Lys, Met, Phe, Pro, STP, Ser, Thr, Trp, Tyr, Val) SELECT entry_id, sample_id, corresponding_gene_call, codon_order_in_gene, consensus, departure_from_consensus, coverage, Ala, Arg, Asn, Asp, Cys, Gln, Glu, Gly, His, Ile, Leu, Lys, Met, Phe, Pro, STP, Ser, Thr, Trp, Tyr, Val FROM variable_amino_acid_frequencies_TEMP;"
    )
    profile_db.cursor.execute("DROP TABLE variable_amino_acid_frequencies_TEMP;")
    # /serious stuff ##################################################################################################

    # set the version
    profile_db.remove_meta_key_value_pair("version")
    profile_db.set_version("14")

    # bye
    profile_db.disconnect()
    progress.end()

    run.info_single(
        "Database successfully upgraded to version 14!",
        nl_after=1,
        nl_before=1,
        mc="green",
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="A simple script to upgrade profile database to from version 13 version 14"
    )
    parser.add_argument(
        "profile_db", metavar="PROFILE_DB", help="Profile database (of version 13)"
    )
    args, unknown = parser.parse_known_args()

    try:
        migrate(args.profile_db)
    except ConfigError as e:
        print(e)
        sys.exit(-1)

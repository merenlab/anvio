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

current_version, next_version = [x[1:] for x in __name__.split("_to_")]

# tables that are specific to profile databases that require updatin'
codons = [
    "AAA",
    "AAC",
    "AAG",
    "AAT",
    "ACA",
    "ACC",
    "ACG",
    "ACT",
    "AGA",
    "AGC",
    "AGG",
    "AGT",
    "ATA",
    "ATC",
    "ATG",
    "ATT",
    "CAA",
    "CAC",
    "CAG",
    "CAT",
    "CCA",
    "CCC",
    "CCG",
    "CCT",
    "CGA",
    "CGC",
    "CGG",
    "CGT",
    "CTA",
    "CTC",
    "CTG",
    "CTT",
    "GAA",
    "GAC",
    "GAG",
    "GAT",
    "GCA",
    "GCC",
    "GCG",
    "GCT",
    "GGA",
    "GGC",
    "GGG",
    "GGT",
    "GTA",
    "GTC",
    "GTG",
    "GTT",
    "TAA",
    "TAC",
    "TAG",
    "TAT",
    "TCA",
    "TCC",
    "TCG",
    "TCT",
    "TGA",
    "TGC",
    "TGG",
    "TGT",
    "TTA",
    "TTC",
    "TTG",
    "TTT",
]
nucleotides = list("ATCGN")

tables = {
    "item_additional_data": {
        "structure": ["item_name", "data_key", "data_value", "data_type", "data_group"],
        "types": ["text", "text", "text", "text", "text"],
    },
    "layer_additional_data": {
        "structure": ["item_name", "data_key", "data_value", "data_type", "data_group"],
        "types": ["text", "text", "text", "text", "text"],
    },
    "variable_codons": {
        "structure": [
            "sample_id",
            "corresponding_gene_call",
            "codon_order_in_gene",
            "reference",
            "departure_from_reference",
            "coverage",
        ]
        + codons,
        "types": ["text", "numeric", "numeric", "text", "numeric", "numeric"]
        + ["numeric"] * len(codons),
    },
    "variable_nucleotides": {
        "structure": [
            "sample_id",
            "split_name",
            "pos",
            "pos_in_contig",
            "corresponding_gene_call",
            "in_noncoding_gene_call",
            "in_coding_gene_call",
            "base_pos_in_codon",
            "codon_order_in_gene",
            "coverage",
            "cov_outlier_in_split",
            "cov_outlier_in_contig",
            "departure_from_reference",
            "competing_nts",
            "reference",
        ]
        + nucleotides,
        "types": [
            "text",
            "text",
            "numeric",
            "numeric",
            "numeric",
            "numeric",
            "numeric",
            "numeric",
            "numeric",
            "numeric",
            "bool",
            "bool",
            "numeric",
            "text",
            "text",
        ]
        + ["numeric"] * len(nucleotides),
    },
    "indels": {
        "structure": [
            "sample_id",
            "split_name",
            "type",
            "sequence",
            "start_in_contig",
            "start_in_split",
            "length",
            "coverage",
        ],
        "types": [
            "text",
            "text",
            "text",
            "text",
            "numeric",
            "numeric",
            "numeric",
            "numeric",
        ],
    },
    "collections_bins_info": {
        "structure": ["collection_name", "bin_name", "source", "html_color"],
        "types": ["text", "text", "text", "text"],
    },
    "collections_of_contigs": {
        "structure": ["collection_name", "contig", "bin_name"],
        "types": ["text", "text", "text"],
    },
    "collections_of_splits": {
        "structure": ["collection_name", "split", "bin_name"],
        "types": ["text", "text", "text"],
    },
}


def drop_entry_id_column_from_table(db_path, table_name, table_properties):
    progress.new("Modifying '%s'" % table_name)

    structure = table_properties["structure"]
    types = table_properties["types"]
    db_fields = ", ".join(["%s %s" % (t[0], t[1]) for t in zip(structure, types)])
    temp_table_name = table_name + "_TEMP"

    _db = db.DB(db_path, None, ignore_version=True)

    progress.update("Creating a temporary table")
    _db._exec("""CREATE TABLE %s (%s)""" % (temp_table_name, db_fields))

    progress.update("Copying data into the temporary table")
    _db._exec(
        """INSERT INTO %s SELECT %s FROM %s"""
        % (temp_table_name, ", ".join(structure), table_name)
    )

    progress.update("Dropping the original table")
    _db._exec("""DROP TABLE IF EXISTS %s""" % (table_name))

    progress.update("Renaming temporary table to become the original")
    _db._exec("""ALTER TABLE %s RENAME TO %s""" % (temp_table_name, table_name))

    progress.update("Committing changes")
    _db.disconnect()

    progress.end()


def migrate(db_path):
    if db_path is None:
        raise ConfigError("No database path is given.")

    # make sure someone is not being funny
    utils.is_profile_db(db_path)

    profile_db = db.DB(db_path, None, ignore_version=True)
    if str(profile_db.get_version()) != current_version:
        profile_db.disconnect()
        raise ConfigError(
            "Version of this profile database is not %s (hence, this script cannot really do anything)."
            % current_version
        )
    profile_db.disconnect()

    # drop entry ids one by one
    for table_name in tables:
        drop_entry_id_column_from_table(
            db_path, table_name, table_properties=tables[table_name]
        )

    # set the version
    profile_db = db.DB(db_path, None, ignore_version=True)
    profile_db.remove_meta_key_value_pair("version")
    profile_db.set_version(next_version)
    profile_db.disconnect()

    progress.end()

    run.info_single(
        "Your profile db is now version %s. %d of its tables were cleaned from a historical "
        "design artifact." % (next_version, len(tables)),
        nl_after=1,
        nl_before=1,
        mc="green",
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="A simple script to upgrade profile database from version %s to version %s"
        % (current_version, next_version)
    )
    parser.add_argument(
        "profile_db",
        metavar="PROFILE_DB",
        help="An anvi'o profile database of version %s" % current_version,
    )
    args, unknown = parser.parse_known_args()

    try:
        migrate(args.profile_db)
    except ConfigError as e:
        print(e)
        sys.exit(-1)

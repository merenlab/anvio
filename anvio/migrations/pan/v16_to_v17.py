#!/usr/bin/env python
# -*- coding: utf-8

import sys
import argparse

import anvio.dbinfo as dbinfo
import anvio.terminal as terminal

from anvio.errors import ConfigError

run = terminal.Run()
progress = terminal.Progress()

current_version, next_version = [x[1:] for x in __name__.split("_to_")]

pan_gene_cluster_function_reactions_table_name = "gene_cluster_function_reactions"
pan_gene_cluster_function_reactions_table_structure = [
    "modelseed_reaction_id",
    "modelseed_reaction_name",
    "ko_kegg_reaction_source",
    "ko_ec_number_source",
    "other_kegg_reaction_ids",
    "other_ec_numbers",
    "metabolite_modelseed_ids",
    "stoichiometry",
    "compartments",
    "reversibility",
]
pan_gene_cluster_function_reactions_table_types = [
    "text",
    "text",
    "text",
    "text",
    "text",
    "text",
    "text",
    "text",
    "text",
    "bool",
]

pan_gene_cluster_function_metabolites_table_name = "gene_cluster_function_metabolites"
pan_gene_cluster_function_metabolites_table_structure = [
    "modelseed_compound_id",
    "modelseed_compound_name",
    "kegg_aliases",
    "formula",
    "charge",
]
pan_gene_cluster_function_metabolites_table_types = [
    "text",
    "text",
    "text",
    "text",
    "numeric",
]


def migrate(db_path):
    if db_path is None:
        raise ConfigError("No database path is given.")

    pan_db_info = dbinfo.PanDBInfo(db_path)
    if str(pan_db_info.version) != current_version:
        raise ConfigError(
            f"The version of the provided pan database is {pan_db_info.version}, not the required version, "
            f"{current_version}, so this script cannot upgrade the database."
        )

    pan_db = pan_db_info.load_db()

    if str(pan_db.get_version()) != current_version:
        raise ConfigError(
            "Version of this contigs database is not %s (hence, this script cannot really do anything)."
            % current_version
        )

    progress.new("Migrating")
    progress.update("Creating two new tables for reactions and metabolites")

    # To be on the safe side, remove any reaction network tables and metadata that might exist.
    try:
        pan_db.drop_table(pan_gene_cluster_function_reactions_table_name)
        pan_db.drop_table(pan_gene_cluster_function_metabolites_table_name)
    except:
        pass

    try:
        pan_db.remove_meta_key_value_pair("reaction_network_ko_annotations_hash")
        pan_db.remove_meta_key_value_pair("reaction_network_kegg_database_release")
        pan_db.remove_meta_key_value_pair("reaction_network_modelseed_database_sha")
        pan_db.remove_meta_key_value_pair("reaction_network_consensus_threshold")
        pan_db.remove_meta_key_value_pair("reaction_network_discard_ties")
    except:
        pass

    pan_db.set_meta_value("reaction_network_ko_annotations_hash", None)
    pan_db.set_meta_value("reaction_network_kegg_database_release", None)
    pan_db.set_meta_value("reaction_network_modelseed_database_sha", None)
    pan_db.set_meta_value("reaction_network_consensus_threshold", None)
    pan_db.set_meta_value("reaction_network_discard_ties", None)
    pan_db.create_table(
        pan_gene_cluster_function_reactions_table_name,
        pan_gene_cluster_function_reactions_table_structure,
        pan_gene_cluster_function_reactions_table_types,
    )
    pan_db.create_table(
        pan_gene_cluster_function_metabolites_table_name,
        pan_gene_cluster_function_metabolites_table_structure,
        pan_gene_cluster_function_metabolites_table_types,
    )

    progress.update("Updating version")
    pan_db.remove_meta_key_value_pair("version")
    pan_db.set_version(next_version)

    progress.update("Committing changes")
    pan_db.disconnect()

    progress.end()

    message = (
        "Congratulations! Your pan database is now version 17, which means it now contains two new empty tables. "
        "These tables are not as boring as they first appear, because you can now run `anvi-reaction-network` to store "
        "a network of the metabolic reactions that may be encoded by gene clusters. A metabolic model representing the "
        "network can be exported from the database using `anvi-get-metabolic-model-file`."
    )
    run.info_single(message, nl_after=1, nl_before=1, mc="green")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="A simple script to upgrade the pan database from version %s to version %s"
        % (current_version, next_version)
    )
    parser.add_argument(
        "pan_db",
        metavar="PAN_DB",
        help="An anvi'o pan database of version %s" % current_version,
    )
    args, unknown = parser.parse_known_args()

    try:
        migrate(args.pan_db)
    except ConfigError as e:
        print(e)
        sys.exit(-1)

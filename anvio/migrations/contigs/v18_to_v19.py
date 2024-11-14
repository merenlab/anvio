#!/usr/bin/env python
# -*- coding: utf-8

import sys
import argparse

import anvio.db as db
import anvio.utils as utils


import anvio.terminal as terminal

from anvio.errors import ConfigError

current_version, next_version = [x[1:] for x in __name__.split("_to_")]

trna_taxonomy_table_name = "trna_taxonomy"
trna_taxonomy_table_structure = [
    "gene_callers_id",
    "amino_acid",
    "anticodon",
    "accession",
    "percent_identity",
    "t_domain",
    "t_phylum",
    "t_class",
    "t_order",
    "t_family",
    "t_genus",
    "t_species",
]
trna_taxonomy_table_types = [
    "numeric",
    "text",
    "text",
    "text",
    "text",
    "text",
    "text",
    "text",
    "text",
    "text",
    "text",
    "text",
]

run = terminal.Run()
progress = terminal.Progress()


def migrate(db_path):
    if db_path is None:
        raise ConfigError("No database path is given.")

    utils.is_contigs_db(db_path)

    contigs_db = db.DB(db_path, None, ignore_version=True)
    if str(contigs_db.get_version()) != current_version:
        raise ConfigError(
            "Version of this contigs database is not %s (hence, this script cannot really do anything)."
            % current_version
        )

    progress.new("Le migrateaoux")
    progress.update("Creating a new table for tRNA taxonomy")

    # just to be on the safe side.
    try:
        contigs_db.drop_table(trna_taxonomy_table_name)
    except:
        pass

    try:
        contigs_db.remove_meta_key_value_pair("trna_taxonomy_was_run")
        contigs_db.remove_meta_key_value_pair("trna_taxonomy_database_version")
    except:
        pass

    contigs_db.set_meta_value("trna_taxonomy_was_run", False)
    contigs_db.set_meta_value("trna_taxonomy_database_version", None)
    contigs_db.create_table(
        trna_taxonomy_table_name,
        trna_taxonomy_table_structure,
        trna_taxonomy_table_types,
    )

    progress.update("Removing tRNA hits")
    relevant_gene_calls = [
        g[0]
        for g in contigs_db._exec(
            "select gene_callers_id from hmm_hits where source='Transfer_RNAs'"
        ).fetchall()
    ]

    if len(relevant_gene_calls):
        for table_name in [
            "hmm_hits_info",
            "hmm_hits",
            "hmm_hits_in_splits",
            "genes_in_contigs",
            "gene_functions",
        ]:
            contigs_db.remove_some_rows_from_table(
                table_name, 'source IN ("Transfer_RNAs")'
            )

        CLAUSE = "gene_callers_id in (%s)" % (
            ",".join([str(x) for x in relevant_gene_calls])
        )
        for table_name in [
            "gene_amino_acid_sequences",
            "genes_taxonomy",
            "genes_in_splits",
        ]:
            contigs_db.remove_some_rows_from_table(table_name, CLAUSE)

        gene_function_sources = contigs_db.get_meta_value("gene_function_sources")
        if len(gene_function_sources):
            new_gene_function_sources = ",".join(
                [f for f in gene_function_sources.split(",") if f != "Transfer_RNAs"]
            )

            contigs_db.remove_meta_key_value_pair("gene_function_sources")
            contigs_db.set_meta_value(
                "gene_function_sources", new_gene_function_sources
            )

    progress.update("Updating version")
    contigs_db.remove_meta_key_value_pair("version")
    contigs_db.set_version(next_version)

    progress.update("Committing changes")
    contigs_db.disconnect()

    progress.end()
    message = (
        f"The contigs database is now {next_version}. This upgrade added an empty tRNA taxonomy table "
        f"in it. Probably you will never use it because who cares about tRNA taxonomy amirite? "
        f"Well, IT WILL FOREVER BE THERE ANYWAY. BOOM."
    )

    if relevant_gene_calls:
        message += (
            " This update also removed tRNA HMM hits from your contis database :/ It was really very "
            "necessary since one of the developers of anvi'o (*COUGH* meren *COUGH*) was very confused "
            "about anticodons. Now everything is fixed, but unfortunately you will have to re-run "
            "`anvi-scan-trnas` program on your contigs database :( "
        )

    run.info_single(message, nl_after=1, nl_before=1, mc="green")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="A simple script to upgrade CONTIGS.db from version %s to version %s"
        % (current_version, next_version)
    )
    parser.add_argument(
        "contigs_db",
        metavar="CONTIGS_DB",
        help="Contigs database at version %s" % current_version,
    )
    args, unknown = parser.parse_known_args()

    try:
        migrate(args.contigs_db)
    except ConfigError as e:
        print(e)
        sys.exit(-1)

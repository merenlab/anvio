#!/usr/bin/env python
# -*- coding: utf-8

import sys
import argparse

import anvio.db as db
import anvio.utils as utils


import anvio.terminal as terminal

from anvio.errors import ConfigError

current_version, next_version = [x[1:] for x in __name__.split('_to_')]


blast_hits_table_name                    = 'blast_hits'
blast_hits_table_structure               = ['match_id' , 'gene_callers_id', 'gene_name', 'taxon_id', 'pourcentage_identity', 'bitscore']
blast_hits_table_types                   = ['text'     ,       'text'    ,      'text'   ,   'text'   ,     'text'   ,         'text']

taxon_names_table_name                   = 'taxon_names'
taxon_names_table_structure              = ['taxon_id', 't_domain', "t_phylum", "t_class", "t_order", "t_family", "t_genus", "t_species"]
taxon_names_table_types                  = [ 'numeric',   'text'  ,   'text'  ,  'text'  ,  'text'  ,  'text'   ,  'text'  ,   'text'   ]

taxonomy_estimation_metagenome_name      = 'taxonomy_estimation_metagenome'
taxonomy_estimation_metagenome_structure = ['gene_caller_id',      'gene_name',  'source' , 't_domain', "t_phylum", "t_class", "t_order", "t_family", "t_genus", "t_species"]
taxonomy_estimation_metagenome_types     = [ 'numeric',             'text'    ,  'text',      'text'  ,   'text'  ,  'text'  ,  'text'   ,  'text'  ,  'text'  ,   'text'   ]



run = terminal.Run()
progress = terminal.Progress()

def migrate(db_path):
    if db_path is None:
        raise ConfigError("No database path is given.")

    utils.is_contigs_db(db_path)

    contigs_db = db.DB(db_path, None, ignore_version = True)
    if str(contigs_db.get_version()) != current_version:
        raise ConfigError("Version of this contigs database is not %s (hence, this script cannot really do anything)." % current_version)

    progress.new("Actualyse databases ")
    progress.update("...")

    contigs_db.create_table(blast_hits_table_name, blast_hits_table_structure, blast_hits_table_types)
    contigs_db.drop_table(taxon_names_table_name)
    #contigs_db._exec('''ALTER TABLE %s ADD %s text''' % (taxon_names_table_name, "t_domain"))
    contigs_db.commit()
    contigs_db.create_table(taxon_names_table_name, taxon_names_table_structure, taxon_names_table_types)
    contigs_db.create_table(taxonomy_estimation_metagenome_name, taxonomy_estimation_metagenome_structure, taxonomy_estimation_metagenome_types)

    progress.update("Updating version")
    contigs_db.remove_meta_key_value_pair('version')
    contigs_db.set_version(next_version)

    progress.update("Committing changes")
    contigs_db.disconnect()

    progress.end()
    run.info_single("The contigs database is now %s. Unfortunatly this update removed Taxonomy \
                     from your contigs database :( We are very sorry about this, but we only did it to be\
                     able to offer you  nicer things." % (next_version), nl_after=1, nl_before=1, mc='green')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='A simple script to upgrade CONTIGS.db from version %s to version %s' % (current_version, next_version))
    parser.add_argument('contigs_db', metavar = 'CONTIGS_DB', help = 'Contigs database at version %s' % current_version)
    args, unknown = parser.parse_known_args()

    try:
        migrate(args.contigs_db)
    except ConfigError as e:
        print(e)
        sys.exit(-1)

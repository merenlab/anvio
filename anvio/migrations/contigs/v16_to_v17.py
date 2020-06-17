#!/usr/bin/env python
# -*- coding: utf-8

import sys
import argparse

import anvio.db as db
import anvio.utils as utils

import anvio.terminal as terminal

from anvio.errors import ConfigError

current_version, next_version = [x[1:] for x in __name__.split('_to_')]

run = terminal.Run()
progress = terminal.Progress()

codons = ['AAA', 'AAC', 'AAG', 'AAT', 'ACA', 'ACC', 'ACG', 'ACT', 'AGA', 'AGC', 'AGG', 'AGT', 'ATA', 'ATC', 'ATG', 'ATT', 'CAA', 'CAC', 'CAG', 'CAT', 'CCA', 'CCC', 'CCG', 'CCT', 'CGA', 'CGC', 'CGG', 'CGT', 'CTA', 'CTC', 'CTG', 'CTT', 'GAA', 'GAC', 'GAG', 'GAT', 'GCA', 'GCC', 'GCG', 'GCT', 'GGA', 'GGC', 'GGG', 'GGT', 'GTA', 'GTC', 'GTG', 'GTT', 'TAA', 'TAC', 'TAG', 'TAT', 'TCA', 'TCC', 'TCG', 'TCT', 'TGA', 'TGC', 'TGG', 'TGT', 'TTA', 'TTC', 'TTG', 'TTT']
nucleotides = list('ATCGN')

tables = {
    'gene_clusters':              {'structure': ['gene_caller_id', 'gene_cluster_id', 'genome_name', 'alignment_summary'],
                                       'types': ['numeric', 'str', 'str', 'str']},
    'genes_in_splits':            {'structure': ['split', 'gene_callers_id', 'start_in_split', 'stop_in_split', 'percentage_in_split'],
                                       'types': ['text', 'numeric', 'numeric', 'numeric', 'numeric']},
    'gene_functions':             {'structure': ['gene_callers_id', 'source', 'accession', 'function', 'e_value'],
                                       'types': ['numeric', 'text', 'text', 'text', 'numeric']},
    'hmm_hits_in_splits':         {'structure': ['hmm_hit_entry_id', 'split', 'percentage_in_split', 'source'],
                                       'types': ['numeric', 'text', 'numeric', 'text']},
    'scg_taxonomy':               {'structure': ['gene_callers_id', 'gene_name', 'accession', 'percent_identity', 't_domain', "t_phylum", "t_class", "t_order", "t_family", "t_genus", "t_species"],
                                       'types': ['numeric', 'text', 'text', 'text', 'text', 'text', 'text', 'text', 'text', 'text', 'text']},
    'nucleotide_additional_data': {'structure': ['item_name', 'data_key', 'data_value', 'data_type', 'data_group'],
                                       'types': ['text', 'text', 'text', 'text', 'text']},
    'amino_acid_additional_data': {'structure': ['item_name', 'data_key', 'data_value', 'data_type', 'data_group'],
                                       'types': ['text', 'text', 'text', 'text', 'text']},
    'gene_level_coverage_stats':  {'structure': ['gene_callers_id', 'sample_name', 'mean_coverage', 'detection', 'non_outlier_mean_coverage', 'non_outlier_coverage_std', 'gene_coverage_values_per_nt', 'non_outlier_positions'],
                                       'types': ['numeric', 'text', 'numeric', 'numeric', 'numeric', 'numeric', 'blob', 'blob']},
    'gene_level_inseq_stats':     {'structure': ['gene_callers_id', 'sample_name', 'mean_coverage', 'insertions', 'insertions_normalized', 'mean_disruption', 'below_disruption', 'gene_coverage_values_per_nt'],
                                       'types': ['numeric', 'text', 'numeric', 'numeric', 'numeric', 'numeric', 'numeric', 'blob']},
    'item_additional_data':       {'structure': ['item_name', 'data_key', 'data_value', 'data_type', 'data_group'],
                                       'types': ['text', 'text', 'text', 'text', 'text']},
    'layer_additional_data':      {'structure': ['item_name', 'data_key', 'data_value', 'data_type', 'data_group'],
                                       'types': ['text', 'text', 'text', 'text', 'text']},
    'variable_codons':            {'structure': ['sample_id', 'corresponding_gene_call', 'codon_order_in_gene', 'reference', 'departure_from_reference', 'coverage']+ codons,
                                       'types': ['text', 'numeric', 'numeric', 'text', 'numeric', 'numeric']+ ['numeric']* len(codons)},
    'variable_nucleotides':       {'structure': ['sample_id', 'split_name', 'pos', 'pos_in_contig', 'corresponding_gene_call', 'in_noncoding_gene_call', 'in_coding_gene_call', 'base_pos_in_codon', 'codon_order_in_gene', 'coverage', 'cov_outlier_in_split', 'cov_outlier_in_contig', 'departure_from_reference', 'competing_nts', 'reference']+ nucleotides,
                                       'types': ['text', 'text', 'numeric', 'numeric', 'numeric', 'numeric', 'numeric', 'numeric', 'numeric', 'numeric', 'bool', 'bool', 'numeric', 'text', 'text']+ ['numeric']* len(nucleotides)},
    'indels':                     {'structure': ['sample_id', 'split_name', 'type', 'sequence', 'start_in_contig', 'start_in_split', 'length', 'coverage'],
                                       'types': ['text', 'text', 'text', 'text', 'numeric', 'numeric', 'numeric', 'numeric']},
    'collections_bins_info':      {'structure': ['collection_name', 'bin_name', 'source', 'html_color'],
                                       'types': ['text', 'text', 'text', 'text']},
    'collections_of_contigs':     {'structure': ['collection_name', 'contig', 'bin_name'],
                                       'types': ['text', 'text', 'text']},
    'collections_of_splits':      {'structure': ['collection_name', 'split', 'bin_name'],
                                       'types': ['text', 'text', 'text']},
    'templates':                  {'structure': ['corresponding_gene_call', 'pdb_id', 'chain_id', 'ppi'],
                                       'types': ['integer', 'text', 'text', 'real']},
    'models':                     {'structure': ['corresponding_gene_call', 'molpdf', 'GA341_score', 'DOPE_score', 'picked_as_best'],
                                       'types': ['integer', 'real', 'real', 'real', 'integer']},
    'residue_info':               {'structure': ['corresponding_gene_call', 'codon_order_in_gene', 'contact_numbers', 'codon', 'amino_acid', 'codon_number', 'sec_struct', 'rel_solvent_acc', 'phi', 'psi', 'NH_O_1_index', 'NH_O_1_energy', 'O_NH_1_index', 'O_NH_1_energy', 'NH_O_2_index', 'NH_O_2_energy', 'O_NH_2_index', 'O_NH_2_energy'],
                                       'types': ['integer', 'integer', 'text', 'text', 'text', 'integer', 'text', 'real', 'real', 'real', 'integer', 'real', 'integer', 'real', 'integer', 'real', 'integer', 'real']}
}


def drop_entry_id_column_from_table(db_path, table_name):
    progress.new("Modifying '%s'" % table_name)

    structure = tables[table_name]['structure']
    types = tables[table_name]['types']
    db_fields = ', '.join(['%s %s' % (t[0], t[1]) for t in zip(structure, types)])
    temp_table_name = table_name + '_TEMP'

    _db = db.DB(db_path, None, ignore_version = True)

    progress.update("Creating a temporary table")
    _db._exec('''CREATE TABLE %s (%s)''' % (temp_table_name, db_fields))

    progress.update("Copying data into the temporary table")
    _db._exec('''INSERT INTO %s SELECT %s FROM %s''' % (temp_table_name, ', '.join(structure), table_name))

    progress.update("Dropping the original table")
    _db._exec('''DROP TABLE IF EXISTS %s''' % (table_name))

    progress.update("Renaming temporary table to become the original")
    _db._exec('''ALTER TABLE %s RENAME TO %s''' % (temp_table_name, table_name))

    progress.update("Committing changes")
    _db.disconnect()

    progress.end()


def migrate(db_path):
    if db_path is None:
        raise ConfigError("No database path is given.")

    utils.is_contigs_db(db_path)

    contigs_db = db.DB(db_path, None, ignore_version = True)
    if str(contigs_db.get_version()) != current_version:
        raise ConfigError("Version of this contigs database is not %s (hence, this script cannot really do anything)." % current_version)

    tables_require_modification = [table_name for table_name in contigs_db.get_table_names() if table_name in tables]
    contigs_db.disconnect()

    # drop entry ids one by one
    if len(tables_require_modification):
        for table_name in tables_require_modification:
            drop_entry_id_column_from_table(db_path, table_name)

    contigs_db = db.DB(db_path, None, ignore_version = True)
    contigs_db.remove_meta_key_value_pair('version')
    contigs_db.set_version(next_version)
    contigs_db.disconnect()

    run.info_single("Your contigs db is now %s. This update carried one more issue into the graveyard "
                    "of bad design decisions we've made years ago by altering %d tables in your database." \
                            % (current_version, len(tables_require_modification)), nl_after=1, nl_before=1, mc='green')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='A simple script to upgrade CONTIGS.db from version %s to version %s' % (current_version, next_version))
    parser.add_argument('contigs_db', metavar = 'CONTIGS_DB', help = 'Contigs database at version %s' % current_version)
    args, unknown = parser.parse_known_args()

    try:
        migrate(args.contigs_db)
    except ConfigError as e:
        print(e)
        sys.exit(-1)

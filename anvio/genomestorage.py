# -*- coding: utf-8
# pylint: disable=line-too-long
"""
    A module to dealing with genome storages.

    Pangenomic workflow heavily uses this module.

    Ad hoc access to make sene of internal or external genome descriptions is also welcome.
"""

import time
import hashlib

import anvio
import anvio.db as db
import anvio.tables as t
import anvio.fastalib as fastalib
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths

from anvio.errors import ConfigError


__author__ = "Developers of anvi'o (see AUTHORS.txt)"
__copyright__ = "Copyleft 2015-2018, the Meren Lab (http://merenlab.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "A. Murat Eren"
__email__ = "a.murat.eren@gmail.com"


run = terminal.Run()
progress = terminal.Progress()
pp = terminal.pretty_print


class GenomeStorage(object):
    def __init__(self, storage_path, storage_hash=None, genome_names_to_focus=None, create_new=False, skip_init_functions=False, run=run, progress=progress):
        self.db_type = 'genomestorage'
        self.version = anvio.__genomes_storage_version__
        self.run = run
        self.storage_hash = storage_hash
        self.progress = progress
        self.storage_path = storage_path
        self.genome_names_to_focus = genome_names_to_focus
        self.skip_init_functions = skip_init_functions

        if create_new:
            self.check_storage_path_for_create_new()
        else:
            self.check_storage_path_for_load()

        self.db = db.DB(self.storage_path, self.version, new_database=create_new)

        self.genome_info_entries = []
        self.gene_info_entries = []
        self.gene_functions_entries = []
        self.gene_functions_entry_id = 0

        if create_new:
            self.create_tables()
        else:
            self.init()


    def check_storage_path_for_create_new(self):
        if not self.storage_path.endswith('GENOMES.db'):
            raise ConfigError("The genomes storage file must end with '-GENOMES.db'. Anvi'o developers do know how ridiculous\
                                this requirement sounds like, but if you have seen the things they did, you would totally\
                                understand why this is necessary.")

        filesnpaths.is_output_file_writable(self.storage_path)


    def check_storage_path_for_load(self):
        if not self.storage_path:
            raise ConfigError("Anvi'o genomes storage is speaking. Someone called the init function,\
                               yet there is nothing to initialize since genome storage path variable\
                               (args.genomes_storage) is None. If you are an end user, please make sure\
                               you provide the genomes storage paramater to whatever program you were\
                               running. If you are a developer, you probably already figured what is\
                               wrong. If you are a cat, you need to send us an e-mail immediately.")

        if self.storage_path.endswith('.h5'):
            raise ConfigError("We recenlty switched from HD5 files (.h5) to Sqlite (.db) files for the genome storage, \
                              you can upgrade your genome storage by running 'anvi-migrate-db %s'." % self.storage_path)

        filesnpaths.is_file_exists(self.storage_path)


    def init(self):
        if self.db_type != self.db.get_meta_value('db_type'):
            raise ConfigError('The database "%s" does not look like a genome storage :/' % self.storage_path)

        if self.storage_hash:
            if self.storage_hash != self.get_storage_hash():
                raise ConfigError("The requested genome storage hash ('%s') does not match with the one read from the database ('%s')." %
                    (self.storage_hash, self.get_storage_hash()))

        self.genome_names_in_db = self.get_all_genome_names()

        if self.genome_names_to_focus:
            genome_names_to_focus_missing_from_db = [g for g in self.genome_names_to_focus if g not in self.genome_names_in_db]

            # make sure the user knows what they're doing
            if genome_names_to_focus_missing_from_db:
                raise ConfigError("%d of %d genome names you wanted to focus are missing from the genomes sotrage.\
                                 Although this may not be a show-stopper, anvi'o likes to be explicit, so here we\
                                 are. Not going anywhere until you fix this. For instance this is one of the missing\
                                 genome names: '%s', and this is one random genome name from the database: '%s'" % \
                                         (len(genome_names_to_focus_missing_from_db), len(self.genome_names_to_focus),\
                                         genome_names_to_focus_missing_from_db[0], ', '.join(self.genome_names_in_db)))

            self.genome_names = self.genome_names_to_focus
        else:
            self.genome_names = self.genome_names_in_db

        self.progress.new('Recovering data from the db')
        self.progress.update('Initializing...')

        self.num_genomes = len(self.genome_names)
        self.functions_are_available = self.db.get_meta_value('functions_are_available')
        self.gene_functions_entry_id = self.db.get_max_value_in_column(t.genome_gene_function_calls_table_name, 'entry_id')

        ## load the data
        self.progress.update('Loading genomes basic info...')
        where_clause = """genome_name IN (%s)""" % ",".join('"' + item + '"' for item in self.genome_names)
        self.genomes_info = self.db.get_some_rows_from_table_as_dict(t.genome_info_table_name, where_clause)

        self.progress.end()

        self.gene_info = {}
        num_genes = self.db.get_row_counts_from_table(t.gene_info_table_name, where_clause)
        self.progress.new('Loading gene info', progress_total_items=num_genes)

        for gene_num, gene_info_tuple in enumerate(self.db.get_some_rows_from_table(t.gene_info_table_name, where_clause)):
            genome_name, gene_caller_id, aa_sequence, dna_sequence, partial, length = gene_info_tuple

            if gene_num % 100 == 0:
                self.progress.increment(increment_to = gene_num)
                self.progress.update('Completed %d/%d' % (gene_num, num_genes))

            if genome_name not in self.gene_info:
                self.gene_info[genome_name] = {}

            self.gene_info[genome_name][gene_caller_id] = {
                'aa_sequence': aa_sequence,
                'dna_sequence': dna_sequence,
                'partial': partial,
                'length': length,
                'functions': {}
            }


            if not self.skip_init_functions:
                functions = self.db.get_some_rows_from_table(t.genome_gene_function_calls_table_name,
                                                             'genome_name = "%s" and gene_callers_id = "%s"' % (genome_name, gene_caller_id))

                for row in functions:
                    self.gene_info[genome_name][gene_caller_id]['functions'][row[3]] = "%s|||%s" % (row[4], row[5])

        self.progress.end()

        self.run.info('Genomes storage', 'Initialized (storage hash: %s)' % (self.get_storage_hash()))
        self.run.info('Num genomes in storage', len(self.get_all_genome_names()))
        self.run.info('Num genomes will be used', len(self.genome_names), mc='green')


    def create_tables(self):
        self.db.create_table(t.genome_info_table_name, t.genome_info_table_structure, t.genome_info_table_types)
        self.db.create_table(t.gene_info_table_name, t.gene_info_table_structure, t.gene_info_table_types)
        self.db.create_table(t.genome_gene_function_calls_table_name, t.genome_gene_function_calls_table_structure, t.genome_gene_function_calls_table_types)
        self.db._exec("CREATE INDEX covering_index ON %s (gene_callers_id, genome_name);" % t.genome_gene_function_calls_table_name)

        self.db.set_meta_value('db_type', self.db_type)
        self.db.set_meta_value('creation_date', time.time())
        self.db.set_meta_value('functions_are_available', False)


    def get_genomes_dict(self):
        # we retrieve all table at once to avoid seperate sql queries
        all_genomes_dict = self.db.get_table_as_dict(t.genome_info_table_name)
        result = {}

        # copy genomes requested by user to result dictionary
        for genome_name in self.genome_names:
            result[genome_name] = all_genomes_dict[genome_name]

        return result


    def update_storage_hash(self):
        # here we create a signature for the storage itself by concatenating all hash values from all genomes. even if one
        # split is added or removed to any of these genomes will change this signature. since we will tie this information
        # to the profile database we will generate for the pangenome analysis, even if one split is added or removed from any
        # of the genomes will make sure that the profile databases from this storage and storage itself are not compatible:

        concatenated_genome_hashes = '_'.join(sorted(map(str, self.db.get_single_column_from_table(t.genome_info_table_name, 'genome_hash'))))
        new_hash = 'hash' + str(hashlib.sha224(concatenated_genome_hashes.encode('utf-8')).hexdigest()[0:8])

        self.db.set_meta_value('hash', new_hash)


    def get_storage_hash(self):
        return str(self.db.get_meta_value('hash'))


    def store_genomes(self, genome_descriptions):
        self.functions_are_available = genome_descriptions.functions_are_available
        self.db.set_meta_value('functions_are_available', self.functions_are_available)
        self.db.set_meta_value('gene_function_sources', ','.join(genome_descriptions.function_annotation_sources))

        num_gene_calls_added_total = 0
        num_partial_gene_calls_total = 0

        genome_names_to_go_through = sorted(genome_descriptions.genomes.keys())

        for genome_name in genome_names_to_go_through:
            self.progress.new('Initializing genomes')
            self.progress.update('%s ...' % genome_name)
            num_gene_calls_added = 0
            num_partial_gene_calls = 0

            genome = genome_descriptions.genomes[genome_name]

            self.add_genome(genome_name, genome)

            functions_dict, aa_sequences_dict, dna_sequences_dict = genome_descriptions.get_functions_and_sequences_dicts_from_contigs_db(genome_name)

            for gene_caller_id in genome['gene_caller_ids']:
                is_partial_gene_call = gene_caller_id in genome['partial_gene_calls']

                self.add_gene_call(genome_name,
                                    gene_caller_id,
                                    aa_sequence=aa_sequences_dict[gene_caller_id]['sequence'],
                                    dna_sequence=dna_sequences_dict[gene_caller_id]['sequence'],
                                    partial=is_partial_gene_call)

                if gene_caller_id in functions_dict:
                    for source in functions_dict[gene_caller_id]:
                        self.add_gene_function_annotation(genome_name, 
                                                          gene_caller_id, 
                                                          source, 
                                                          functions_dict[gene_caller_id][source])

                num_gene_calls_added += 1
                if is_partial_gene_call:
                    num_partial_gene_calls += 1


            self.progress.end()
            self.run.info_single('%s is stored with %s genes (%s of which were partial)' % (genome_name, pp(num_gene_calls_added), pp(num_partial_gene_calls)),
                          cut_after=120,
                          nl_before = 1 if genome_name == genome_names_to_go_through[0] else 0,
                          nl_after  = 1 if genome_name == genome_names_to_go_through[-1] else 0)

            num_gene_calls_added_total += num_gene_calls_added
            num_partial_gene_calls_total += num_partial_gene_calls

        # write entries to the database.
        self.db.insert_many(t.genome_info_table_name, entries=self.genome_info_entries)
        self.db.insert_many(t.gene_info_table_name, entries=self.gene_info_entries)
        self.db.insert_many(t.genome_gene_function_calls_table_name, entries=self.gene_functions_entries)
        self.update_storage_hash()

        self.run.info('The new genomes storage', '%s (v%s, signature: %s)' % (self.storage_path, self.version, self.get_storage_hash()))
        self.run.info('Number of genomes', '%s (internal: %s, external: %s)' % (pp(len(genome_descriptions.genomes)), 
                                                                                pp(len(genome_descriptions.internal_genome_names)), 
                                                                                pp(len(genome_descriptions.external_genome_names))))
        self.run.info('Number of gene calls', '%s' % pp(num_gene_calls_added_total))
        self.run.info('Number of partial gene calls', '%s' % pp(num_partial_gene_calls_total))

        self.close()


    def add_genome(self, genome_name, genome_info_dict):
        values = (genome_name, )

        for column_name in t.genome_info_table_structure[1:]:
            if genome_info_dict[column_name]:
                values += (genome_info_dict[column_name], )
            else:
                # the following line will add a -1 for any `key` that has the value of `None`. the reason
                # we added this was to be able to work with contigs databases without any hmm hits for SCGs
                # which is covered in https://github.com/merenlab/anvio/issues/573
                values += (-1, )

        self.genome_info_entries.append(values)


    def add_gene_call(self, genome_name, gene_caller_id, aa_sequence, dna_sequence, partial=0):
        self.gene_info_entries.append((genome_name, gene_caller_id, aa_sequence, dna_sequence,
                                                partial, len(aa_sequence),))


    def add_gene_function_annotation(self, genome_name, gene_caller_id, source, annotation):
        if not annotation or len(annotation) != 3:
            return

        accession, function, e_value = annotation
        values = (genome_name, self.gene_functions_entry_id, gene_caller_id, source, accession, function, e_value, )

        self.gene_functions_entry_id += 1
        self.gene_functions_entries.append(values)


    def is_known_genome(self, genome_name, throw_exception=True):
        if genome_name not in self.genomes_info:
            if throw_exception:
                raise ConfigError('The database at "%s" does not know anything about "%s" :(' % (self.storage_path, genome_name))
            else:
                return False
        else:
            return True


    def is_known_gene_call(self, genome_name, gene_caller_id):
        if genome_name not in self.gene_info and gene_caller_id not in self.gene_info[genome_name]:
            raise ConfigError('The database at "%s" does not know anything gene caller id "%d" in genome "%s" :(' % (self.storage_path, gene_caller_id, genome_name))


    def is_partial_gene_call(self, genome_name, gene_caller_id):
        self.is_known_genome(genome_name)
        self.is_known_gene_call(genome_name, gene_caller_id)

        return (self.gene_info[genome_name][gene_caller_id]['partial'] == 1)


    def get_gene_caller_ids(self, genome_name):
        self.is_known_genome(genome_name)
        return self.gene_info[genome_name].keys()


    def get_gene_sequence(self, genome_name, gene_caller_id, report_DNA_sequences=False):
        """Returns gene amino acid sequence unless `report_DNA_sequences` is True."""
        self.is_known_genome(genome_name)
        self.is_known_gene_call(genome_name, gene_caller_id)

        column_name = 'dna_sequence' if report_DNA_sequences else 'aa_sequence' 

        return self.gene_info[genome_name][gene_caller_id][column_name]


    def get_gene_functions(self, genome_name, gene_callers_id):
        if not self.functions_are_available:
            raise ConfigError("Functions are not available in this genome storage ('%s'). " % self.storage_path)

        if self.skip_init_functions:
            raise ConfigError("Initialization of functions were skipped when the GenomeStorage\
                              class was called for '%s'. " % self.storage_path)

        self.is_known_genome(genome_name)
        self.is_known_gene_call(genome_name, gene_callers_id)

        return self.gene_info[genome_name][gene_callers_id]['functions']


    def get_all_genome_names(self):
        return self.db.get_single_column_from_table(t.genome_info_table_name, 'genome_name')


    def gen_combined_aa_sequences_FASTA(self, output_file_path, exclude_partial_gene_calls=False):
        self.run.info('Exclude partial gene calls', exclude_partial_gene_calls, nl_after=1)

        total_num_aa_sequences = 0
        total_num_excluded_aa_sequences = 0

        fasta_output = fastalib.FastaOutput(output_file_path)

        genome_info_dict = self.get_genomes_dict()

        for genome_name in self.genome_names:
            self.progress.new('Storing aa sequences')
            self.progress.update('%s ...' % genome_name)

            gene_caller_ids = sorted([int(gene_caller_id) for gene_caller_id in self.get_gene_caller_ids(genome_name)])

            for gene_caller_id in gene_caller_ids:
                is_partial = self.is_partial_gene_call(genome_name, gene_caller_id)

                if exclude_partial_gene_calls and is_partial:
                    total_num_excluded_aa_sequences += 1
                    continue

                aa_sequence = self.get_gene_sequence(genome_name, gene_caller_id)

                fasta_output.write_id('%s_%d' % (genome_info_dict[genome_name]['genome_hash'], int(gene_caller_id)))
                fasta_output.write_seq(aa_sequence, split=False)

                total_num_aa_sequences += 1

            self.progress.end()

        fasta_output.close()

        self.run.info('AA sequences FASTA', output_file_path)
        self.run.info('Num AA sequences reported', '%s' % pp(total_num_aa_sequences), nl_before=1)
        self.run.info('Num excluded gene calls', '%s' % pp(total_num_excluded_aa_sequences))

        return total_num_aa_sequences, total_num_excluded_aa_sequences


    def close(self):
        self.db.disconnect()

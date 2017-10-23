# -*- coding: utf-8
# pylint: disable=line-too-long
"""
    A module to dealing with genome storages.

    Pangenomic workflow heavily uses this module.

    Ad hoc access to make sene of internal or external genome descriptions is also welcome.
"""

import time
import hashlib
import argparse

import anvio
import anvio.db as db
import anvio.tables as t
import anvio.utils as utils
import anvio.dbops as dbops
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths
import anvio.auxiliarydataops as auxiliarydataops

from anvio.errors import ConfigError
from anvio.genomedescriptions import GenomeDescriptions


__author__ = "A. Murat Eren"
__copyright__ = "Copyright 2017, The anvio Project"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "A. Murat Eren"
__email__ = "a.murat.eren@gmail.com"


run = terminal.Run()
progress = terminal.Progress()
pp = terminal.pretty_print


class GenomeStorage(object):
    def __init__(self, storage_path, create_new=False, run=run, progress=progress):
        self.db_type = 'genomes data storage'
        self.version = anvio.__genomes_storage_version__
        self.run = run
        self.progress = progress

        if create_new:
            self.storage_path = storage_path or 'GENOMES.db'

        self.check_storage_path()
        self.db = db.DB(self.storage_path, self.version, new_database=create_new)

        if create_new:
            self.create_tables()


    def check_storage_path(self):
        if not self.storage_path.endswith('GENOMES.db'):
            raise ConfigError("The genomes storage file must end with '-GENOMES.db'. Anvi'o developers do know how ridiculous\
                                this requirement sounds like, but if you have seen the things they did, you would totally\
                                understand why this is necessary.")

        filesnpaths.is_output_file_writable(self.storage_path)


    def create_tables(self):
        self.db.create_table(t.genome_info_table_name, t.genome_info_table_structure, t.genome_info_table_types)
        self.db.create_table(t.gene_info_table_name, t.gene_info_table_structure, t.gene_info_table_types)
        self.db.create_table(t.genome_gene_function_calls_table_name, t.genome_gene_function_calls_table_structure, t.genome_gene_function_calls_table_types)

        self.db.set_meta_value('db_type', self.db_type)
        self.db.set_meta_value('creation_date', time.time())


    def init_genomes_data_storage(self):
        """Initializes an existing genomes storage by reading everything about genomes of interest"""


        self.storage_path = A('genomes_storage')
        self.genome_names_to_focus = A('genome_names')

        if not self.storage_path:
            raise ConfigError("Anvi'o genomes storage is speaking. Someone called the init function,\
                               yet there is nothing to initialize since genome storage path variable\
                               (args.genomes_storage) is None. If you are an end user, please make sure\
                               you provide the genomes storage paramater to whatever program you were\
                               running. If you are a developer, you probably already figured what is\
                               wrong. If you are a cat, you need to send us an e-mail immediately.")

        # let's take care of the genome names to focus, if there are any, first. 
        if self.genome_names_to_focus:
            if filesnpaths.is_file_exists(self.genome_names_to_focus, dont_raise=True):
                self.genome_names_to_focus = utils.get_column_data_from_TAB_delim_file(self.genome_names_to_focus, column_indices=[0], expected_number_of_fields=1)[0]
            else:
                self.genome_names_to_focus = [g.strip() for g in self.genome_names_to_focus.split(',')]

            self.run.warning("A subset of genome names is found, and anvi'o will focus only on to those.")

        filesnpaths.is_proper_genomes_storage_file(self.storage_path)

        self.genomes_storage = auxiliarydataops.GenomesDataStorage(self.storage_path, db_hash = None, genome_names_to_focus=self.genome_names_to_focus, ignore_hash = True)
        self.genomes_storage_hash = self.genomes_storage.get_storage_hash()

        self.genomes = self.genomes_storage.get_genomes_dict()

        self.external_genome_names = [g for g in self.genomes if self.genomes[g]['external_genome']]
        self.internal_genome_names = [g for g in self.genomes if not self.genomes[g]['external_genome']]

        for genome_name in self.genomes:
            self.hash_to_genome_name[self.genomes[genome_name]['genome_hash']] = genome_name


    def gen_storage_hash(self):
        # here we create a signature for the storage itself by concatenating all hash values from all genomes. even if one
        # split is added or removed to any of these genomes will change this signature. since we will tie this information
        # to the profile database we will generate for the pangenome analysis, even if one split is added or removed from any
        # of the genomes will make sure that the profile databases from this storage and storage itself are not compatible:

        concatenated_genome_hashes = '_'.join(self.genomes[genome_name]['genome_hash'] for genome_name in self.genomes)
        return hashlib.sha224(concatenated_genome_hashes).encode('utf-8').hexdigest()[0:8]


    def store_genomes(self, genome_descriptions):
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

            functions_dict, aa_sequences_dict, dna_sequences_dict = self.get_functions_and_sequences_dicts_from_contigs_db(genome['contigs_db_path'], 
                                                                                                                           genome['gene_caller_ids'],
                                                                                                                           functions_are_available=genome_descriptions.functions_are_available,
                                                                                                                           function_annotation_sources=genome_descriptions.function_annotation_sources)

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

        self.run.info('The new genomes storage', '%s (v%s, signature: %s)' % (self.storage_path, self.version, 'storage_hash')) # TO DO
        self.run.info('Number of genomes', '%s (internal: %s, external: %s)' % (pp(len(genome_descriptions.genomes)), 
                                                                                pp(len(genome_descriptions.internal_genome_names)), 
                                                                                pp(len(genome_descriptions.external_genome_names))))
        self.run.info('Number of gene calls', '%s' % pp(num_gene_calls_added_total))
        self.run.info('Number of partial gene calls', '%s' % pp(num_partial_gene_calls_total))


    def get_functions_and_sequences_dicts_from_contigs_db(self, contigs_db_path, gene_caller_ids=None, functions_are_available=False, function_annotation_sources=set([])):
        """Returns function calls, dna and amino acid sequences for `gene_caller_ids`
           from a contigs database"""

        args = argparse.Namespace(contigs_db=contigs_db_path)
        contigs_super = dbops.ContigsSuperclass(args, r=anvio.terminal.Run(verbose=False))

        if functions_are_available:
            contigs_super.init_functions(requested_sources=function_annotation_sources)
            function_calls_dict = contigs_super.gene_function_calls_dict
        else:
            function_calls_dict = {}

        # get dna sequences
        gene_caller_ids_list, dna_sequences_dict = contigs_super.get_sequences_for_gene_callers_ids(gene_caller_ids_list=list(gene_caller_ids))

        # get amino acid sequences.
        # FIXME: this should be done in the contigs super.
        contigs_db = dbops.ContigsDatabase(contigs_db_path)
        aa_sequences_dict = contigs_db.db.get_table_as_dict(t.gene_protein_sequences_table_name)
        contigs_db.disconnect()

        return (function_calls_dict, aa_sequences_dict, dna_sequences_dict)


    def add_genome(self, genome_name, genome_info_dict):
        if self.is_known_genome(genome_name, throw_exception=False):
            raise ConfigError("Genome '%s' is already in this data storage :/" % genome_name)

        values = (genome_name, )

        for column_name in t.genome_info_table_structure[1:]:
            if genome_info_dict[column_name]:
                values += (genome_info_dict[column_name], )
            else:
                # the following line will add a -1 for any `key` that has the value of `None`. the reason
                # we added this was to be able to work with contigs databases without any hmm hits for SCGs
                # which is covered in https://github.com/merenlab/anvio/issues/573
                values += (-1, )

        self.db.insert(t.genome_info_table_name, values=values)


    def add_gene_call(self, genome_name, gene_caller_id, aa_sequence, dna_sequence, partial=0):
        self.db.insert(t.gene_info_table_name, (genome_name, gene_caller_id, aa_sequence, dna_sequence,
                                                partial, len(aa_sequence),))


    def add_gene_function_annotation(self, genome_name, gene_caller_id, source, annotation):
        if not annotation or len(annotation) != 3:
            return

        accession, function, e_value = annotation
        values = (genome_name, 0, gene_caller_id, source, accession, function, e_value, )

        self.db.insert(t.gene_function_calls_table_name, values=values)


    def is_known_genome(self, genome_name, throw_exception=True):
        cursor = self.db._exec('''SELECT genome_name FROM %s WHERE genome_name = ?''' % t.genome_info_table_name, (genome_name, ))
        rows = cursor.fetchall()

        if len(rows) == 0:
            if throw_exception:
                raise ConfigError('The database at "%s" does not know anything about "%s" :(' % (self.storage_path, genome_name))
            else:
                return False
        else:
            return True


    def close(self):
        self.db.disconnect()
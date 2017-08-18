# -*- coding: utf-8
# pylint: disable=line-too-long
"""
    A module to dealing with genome storages.

    Pangenomic workflow heavily uses this module.

    Ad hoc access to make sene of internal or external genome descriptions is also welcome.
"""

import hashlib
import argparse

import anvio
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


class GenomeStorage(GenomeDescriptions):
    def __init__(self, args=None, run=run, progress=progress):
        self.args = args
        self.run = run
        self.progress = progress

        GenomeDescriptions.__init__(self, args, self.run, self.progress)

        self.storage_path = None
        self.genomes_storage = None
        self.genome_names_to_focus = None

        self.hash_to_genome_name = {}

        self.functions_are_available = False
        self.function_annotation_sources = set([])


    def init_genomes_data_storage(self):
        """Initializes an existing genomes storage by reading everything about genomes of interest"""

        A = lambda x: self.args.__dict__[x] if x in self.args.__dict__ else None
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


    def get_functions_and_sequences_dicts_from_contigs_db(self, contigs_db_path, gene_caller_ids=None):
        """Returns function calls, dna and amino acid sequences for `gene_caller_ids`
           from a contigs database"""

        args = argparse.Namespace(contigs_db=contigs_db_path)
        contigs_super = dbops.ContigsSuperclass(args, r=anvio.terminal.Run(verbose=False))

        # get functions
        if self.functions_are_available:
            contigs_super.init_functions(requested_sources=self.function_annotation_sources)
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


    def create_genomes_data_storage(self):
        """Creates an HDF5 file storing all genome related information for later access."""

        # some simple checks
        A = lambda x: self.args.__dict__[x] if x in self.args.__dict__ else None
        self.storage_path = A('output_file')
        if not self.storage_path:
            self.storage_path = "GENOMES.h5"
        else:
            if not self.storage_path.endswith('-GENOMES.h5'):
                raise ConfigError("The genomes storage file must end with '-GENOMES.h5'. Anvi'o developers do know how ridiculous\
                                    this requirement sounds like, but if you have seen the things they did, you would totally\
                                    understand why this is necessary.")

        filesnpaths.is_output_file_writable(self.storage_path)

        # let's read those internal and external genome files the user sent. the following function is inherited from
        # the GenomeDescriptions class:
        self.load_genomes_descriptions()

        # here we create a signature for the storage itself by concatenating all hash values from all genomes. even if one
        # split is added or removed to any of these genomes will change this signature. since we will tie this information
        # to the profile database we will generate for the pangenome analysis, even if one split is added or removed from any
        # of the genomes will make sure that the profile databases from this storage and storage itself are not compatible:
        storage_hash = hashlib.sha224('_'.join(self.genomes[genome_name]['genome_hash'] for genome_name in self.genomes).encode('utf-8')).hexdigest()[0:8]

        self.genomes_storage = auxiliarydataops.GenomesDataStorage(self.storage_path, storage_hash, create_new=True)
        self.genomes_storage.fp.attrs['functions_are_available'] = self.functions_are_available

        # some silly stuff for later fun
        num_gene_calls_added_total = 0
        num_partial_gene_calls_total = 0

        # main loop
        genome_names_to_go_through = sorted(self.genomes.keys())
        for genome_name in genome_names_to_go_through:
            self.progress.new('Initializing genomes')
            self.progress.update('%s ...' % genome_name)
            g = self.genomes[genome_name]

            self.genomes_storage.add_genome(genome_name, g)

            num_gene_calls_added = 0
            num_partial_gene_calls = 0

            functions_dict, aa_sequences_dict, dna_sequences_dict = self.get_functions_and_sequences_dicts_from_contigs_db(g['contigs_db_path'], g['gene_caller_ids'])

            for gene_caller_id in g['gene_caller_ids']:
                partial_gene_call = gene_caller_id in g['partial_gene_calls']

                functions = []
                if gene_caller_id in functions_dict:
                    for annotation_source in self.function_annotation_sources:
                        if annotation_source in functions_dict[gene_caller_id]:
                            annotation_tuple = functions_dict[gene_caller_id][annotation_source]
                            if annotation_tuple:
                                functions.append((annotation_source, '%s|||%s' % (annotation_tuple[0], annotation_tuple[1])),)

                self.genomes_storage.add_gene_call_data(genome_name,
                                                        gene_caller_id,
                                                        aa_sequence=aa_sequences_dict[gene_caller_id]['sequence'],
                                                        dna_sequence=dna_sequences_dict[gene_caller_id]['sequence'],
                                                        partial=partial_gene_call,
                                                        functions=functions)

                num_gene_calls_added += 1
                if partial_gene_call:
                    num_partial_gene_calls += 1

            self.progress.end()

            if genome_names_to_go_through.index(genome_name) == 0:
                self.run.info_single('%s is stored with %s genes (%s of which were partial)' % (genome_name, pp(num_gene_calls_added), pp(num_partial_gene_calls)), cut_after=120, nl_before = 1)
            elif genome_names_to_go_through.index(genome_name) == len(genome_names_to_go_through) -1:
                self.run.info_single('%s is stored with %s genes (%s of which were partial)' % (genome_name, pp(num_gene_calls_added), pp(num_partial_gene_calls)), cut_after=120, nl_after = 1)
            else:
                self.run.info_single('%s is stored with %s genes (%s of which were partial)' % (genome_name, pp(num_gene_calls_added), pp(num_partial_gene_calls)), cut_after=120)

            num_gene_calls_added_total += num_gene_calls_added
            num_partial_gene_calls_total += num_partial_gene_calls

        self.run.info('The new genomes storage', '%s (v%s, signature: %s)' % (self.storage_path, self.genomes_storage.version, storage_hash))
        self.run.info('Number of genomes', '%s (internal: %s, external: %s)' % (pp(len(self.genomes)), pp(len(self.internal_genome_names)), pp(len(self.external_genome_names))))
        self.run.info('Number of gene calls', '%s' % pp(num_gene_calls_added_total))
        self.run.info('Number of partial gene calls', '%s' % pp(num_partial_gene_calls_total))

        self.genomes_storage.close()

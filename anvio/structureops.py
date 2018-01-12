# -*- coding: utf-8
# pylint: disable=line-too-long

"""Classes to make sense of genes and variability within the context of protein structure"""

import os
import shutil

import pandas as pd
import anvio.db as db
import anvio.dbops as dbops
import anvio.fastalib as u
import anvio.utils as utils
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths
import anvio.drivers.MODELLER as MODELLER

from anvio.errors import ConfigError

class StructureDatabase(object):
    def __init__(self, file_path, db_hash, create_new=False, ignore_hash=False, run=terminal.Run(), progress=terminal.Progress(), quiet=False):
        self.db_type = 'structure'
        self.db_hash = str(db_hash)
        self.version = anvio.__auxiliary_data_version__
        self.file_path = file_path
        self.quiet = quiet
        self.run = run
        self.progress = progress
        self.entries = []

        self.db = db.DB(self.file_path, self.version, new_database=create_new)

        if create_new:
            self.create_tables()

        if not ignore_hash:
            self.check_hash()

    def create_tables(self):
        self.db.set_meta_value('db_type', self.db_type)
        self.db.set_meta_value('profile_db_hash', self.db_hash)
        self.db.set_meta_value('creation_date', time.time())

        self.db.create_table(t.structure_table_name, t.structure_table_structure, t.structure_table_types)


    def check_hash(self):
        actual_db_hash = str(self.db.get_meta_value('profile_db_hash'))
        if self.db_hash != actual_db_hash:
            raise ConfigError('The hash value inside Structure Database "%s" does not match with Profile Database hash "%s",\
                                      this files probaby belong to different projects.' % (actual_db_hash, self.db_hash))


    def append(self, gene_id, pdb_path):
        pdb_content = open(pdb_path, 'rb').read()
        self.entries.append((gene_id, pdb_content))


    def store(self):
        self.db.insert_many(t.structure_table_name, entries=self.entries)
        self.entries = []

    def close(self):
        self.db.disconnect()


class Structure(object):

    def __init__(self, args, run=terminal.Run(), progress=terminal.Progress()):
        self.args = args
        self.run = run
        self.progress = progress

        # initialize self.arg parameters
        A = lambda x, t: t(args.__dict__[x]) if x in self.args.__dict__ else None
        null = lambda x: x
        self.contigs_db_path = A('contigs_db', null)
        self.genes_of_interest_path = A('genes_of_interest', null)
        self.splits_of_interest_path = A('splits_of_interest', null)
        self.bin_id = A('bin_id', null)
        self.collection_name = A('collection_name', null)
        self.gene_caller_ids = A('gene_caller_ids', null)
        self.output_dir = A('output_dir', null)
        self.full_output = A('black_no_sugar', bool)

        # MODELLER params
        self.modeller_database = A('database_name', null)
        self.best = A('best', null)
        self.max_matches = A('max_number_templates', null)
        self.min_proper_pident = A('percent_identical_cutoff', null)
        self.num_models = A('num_models', null)
        self.deviation = A('deviation', null)
        self.very_fast = A('very_fast', bool)

        # check output and define absolute path
        self.output_dir = filesnpaths.check_output_directory(self.output_dir, ok_if_exists=False)

        # identify which genes user wants to model structures for
        self.get_genes_of_interest()

        self.sanity_check()
    

    def sanity_check(self):

        # check for genes that do not appear in the contigs database
        bad_gene_caller_ids = [g for g in self.genes_of_interest if g not in self.genes_in_database]
        if bad_gene_caller_ids:
            raise ConfigError(("This gene caller id you provided is" if len(bad_gene_caller_ids) == 1 else \
                               "These gene caller ids you provided are") + " not known to this contigs database: {}.\
                               You have only 2 lives left. 2 more mistakes, and anvi'o will automatically uninstall \
                               itself. Yes, seriously :(".format(", ".join([str(x) for x in bad_gene_caller_ids])))

        # Finally, raise warning if number of genes is greater than 20
        if len(self.genes_of_interest) > 20:
            self.run.warning("Modelling protein structures is no joke. The number of genes you want protein structures for is \
                              {}, which is a lot (of time!). If its taking too long, consider using the --very-fast flag. \
                              CTRL + C to cancel.".format(len(self.genes_of_interest)))

        # if self.percent_identical_cutoff is < 25, you should be careful about accuracy of models
        if self.min_proper_pident < 25:
            self.run.warning("You selected a percent identical cutoff of {}%. Below 25%, you should pay close attention \
                              to the quality of the proteins...".format(self.min_proper_pident))


    def get_genes_of_interest(self):
        """
        nabs the genes of interest based on user arguments (self.args)
        """
        self.genes_of_interest = None

        # identify the gene caller ids of all genes available
        self.genes_in_database = set(dbops.ContigsSuperclass(self.args).genes_in_splits.keys())

        if not self.genes_in_database:
            raise ConfigError("This contigs database does not contain any identified genes...")

        # settling genes of interest
        if self.genes_of_interest_path and self.gene_caller_ids:
            raise ConfigError("You can't provide a gene caller id from the command line, and a list of gene caller ids\
                               as a file at the same time, obviously.")

        if self.gene_caller_ids:
            self.gene_caller_ids = set([x.strip() for x in self.gene_caller_ids.split(',')])

            self.genes_of_interest = []
            for gene in self.gene_caller_ids:
                try:
                    self.genes_of_interest.append(int(gene))
                except:
                    raise ConfigError("Anvi'o does not like your gene caller id '%s'..." % str(gene))

            self.genes_of_interest = set(self.genes_of_interest)

        elif self.genes_of_interest_path:
            filesnpaths.is_file_tab_delimited(self.genes_of_interest_path, expected_number_of_fields=1)

            try:
                self.genes_of_interest = set([int(s.strip()) for s in open(self.genes_of_interest_path).readlines()])
            except ValueError:
                raise ConfigError("Well. Anvi'o was working on your genes of interest ... and ... those gene IDs did not\
                                   look like anvi'o gene caller ids :/ Anvi'o is now sad.")

        if not self.genes_of_interest:
            # no genes of interest are specified. Assuming all, which could be innumerable--raise warning
            self.genes_of_interest = self.genes_in_database
            self.run.warning("You did not specify any genes of interest, so anvi'o will assume all of them are of interest.")


    def process(self):
        """
        This is the workflow for a standard protein search for homologous templates and then
        modelling the structure of the target protein using the homologous protein structures.
        """

        for gene in self.genes_of_interest:

            # MODELLER outputs a lot of stuff into its working directory. A temporary directory is made
            # for each instance of MODELLER (i.e. each protein), and files are moved into
            # self.output_dir afterwards. If --black-no-sugar is provided, everything is moved.
            # Otherwise, only pertinent files are moved. See move_results_to_output_dir()
            self.args.directory = filesnpaths.get_temp_directory_path()
            self.args.target_fasta_path = filesnpaths.get_temp_file_path()

            self.run.warning("Working directory: {}".format(self.args.directory),
                              header='MODELLING STRUCTURE FOR GENE ID {}'.format(gene),
                              lc="cyan")

            dbops.export_aa_sequences_from_contigs_db(self.contigs_db_path, self.args.target_fasta_path, set([gene]), quiet=True)

            self.run_modeller()
            #self.move_results_to_output_dir()

    # This is now deprecated since we are going to use a structured database. But it may be useful for
    # an Nvio program that will be for exporting the structure db into a folder structure
    def move_results_to_output_dir(self):
        """
        if --black-no-sugar, all files from MODELLERs directory are recursively moved into
        output_gene_dir.  Otherwise, the list of files we care about are defined in this function
        and moved into output_gene_dir.
        """
        output_gene_dir = os.path.join(self.output_dir, self.modeller.gene_id)
        filesnpaths.check_output_directory(output_gene_dir)

        if self.full_output:
            shutil.move(self.modeller.directory, output_gene_dir)

        else:
            filesnpaths.gen_output_directory(output_gene_dir)
            list_to_keep = [self.best_structure_filepath,
                            self.modeller.alignment_pap_path,
                            self.modeller.search_results_path,
                            self.modeller.model_info_path]
            for filepath in list_to_keep:
                shutil.move(filepath, output_gene_dir)

        self.run.warning("results folder: {}".format(output_gene_dir),
                     header='FINISHED STRUCTURE FOR GENE ID {}'.format(self.modeller.gene_id),
                     lc="green")


    def run_modeller(self):

        self.modeller = MODELLER.MODELLER(self.args, run=self.run, progress=self.progress)
        self.modeller.process()

#        self.run.warning("Working directory: {}".format(self.modeller.directory),
#                     header='MODELLING STRUCTURE FOR GENE ID {}'.format(self.modeller.gene_id),
#                     lc="cyan")
#
#        try:
#            self.modeller.run_fasta_to_pir()
#
#            self.modeller.check_database()
#
#            self.modeller.run_search()
#
#            self.parse_search_results()
#
#            self.download_structures()
#
#            self.modeller.run_align_to_templates(self.top_seq_matches)
#
#            self.modeller.run_get_model(self.num_models, self.deviation, self.very_fast)
#
#            self.modeller.tidyup()
#
#            if not self.full_output:
#                self.best_structure_filepath = self.modeller.pick_best_model()
#
#            self.modeller.rewrite_model_info()
#
#        except self.modeller.EndModeller as e:
#            print(e)
#            self.modeller.abort()

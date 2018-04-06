# -*- coding: utf-8
# pylint: disable=line-too-long

"""Classes to make sense of genes and variability within the context of protein structure"""

import os
import time
import shutil

import anvio
import anvio.db as db
import anvio.tables as t
import anvio.fastalib as u
import anvio.utils as utils
import anvio.dbops as dbops
import anvio.terminal as terminal
import anvio.constants as constants
import anvio.filesnpaths as filesnpaths
import anvio.drivers.MODELLER as MODELLER

import numpy as np
import pandas as pd

from anvio.errors import ConfigError
from Bio.PDB import PDBParser
from Bio.PDB import DSSP


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
        self.output_db_path = A('output_db_path', null)
        self.full_output = A('dump_dir', null)
        self.skip_DSSP = A('skip_DSSP', bool)
        self.DSSP_executable = None

        # MODELLER params
        self.modeller_database = A('database_name', null)
        self.best = A('best', null)
        self.max_matches = A('max_number_templates', null)
        self.min_proper_pident = A('percent_identical_cutoff', null)
        self.num_models = A('num_models', null)
        self.deviation = A('deviation', null)
        self.very_fast = A('very_fast', bool)

        # check outputs are writable
        filesnpaths.is_output_file_writable(self.output_db_path)
        if self.full_output:
            self.full_output = filesnpaths.check_output_directory(self.full_output, ok_if_exists=False)

        # initialize residue annotation dataframe
        self.residue_annotation_df = pd.DataFrame({})

        # list of methods to be used for residue annotation table
        self.list_of_residue_annotation_methods = [self.run_DSSP]
        if self.skip_DSSP:
            del self.list_of_residue_annotation_methods[self.list_of_residue_annotation_methods.index(self.run_DSSP)]

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

        # check that DSSP exists
        if self.skip_DSSP:
            self.run.warning("You requested to skip amino acid residue annotation with DSSP. A bold move only an expert could justify... \
                              Anvi'o's respect for you increases slightly.")

        else:
            if utils.is_program_exists("mkdssp", dont_raise=True): # mkdssp is newer and preferred
                self.DSSP_executable = "mkdssp"

            if not self.DSSP_executable:
                if utils.is_program_exists("dssp", dont_raise=True):
                    self.DSSP_executable = "dssp"
                else:
                    raise ConfigError("An anvi'o function needs 'mkdssp' or 'dssp' to be installed on your system, but\
                                       neither seem to appear in your path :/ If you are certain you have either on your\
                                       system (for instance you can run either by typing 'mkdssp' or 'dssp' in your terminal\
                                       window), you may want to send a detailed bug report. If you want to skip secondary\
                                       structure and solvent accessibility annotation, provide the flag --skip-DSSP.")

            self.run.info_single("Anvi'o found the DSSP executable `%s`, and will use it."\
                                  % self.DSSP_executable, nl_before=1, nl_after=1)


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
        """

        self.structure_db = StructureDatabase(self.output_db_path, "DUMMYHASH", create_new=True) # FIXME

        for gene_id in self.genes_of_interest:

            # MODELLER outputs a lot of stuff into its working directory. A temporary directory is
            # made for each instance of MODELLER (i.e. each protein), And bits and pieces of this
            # directory are used in the creation of the structure database. If self.full_output is
            # provided, these directories and their contents are moved into self.full_output.
            self.args.directory = filesnpaths.get_temp_directory_path()
            self.args.target_fasta_path = filesnpaths.get_temp_file_path()

            dbops.export_aa_sequences_from_contigs_db(self.contigs_db_path, self.args.target_fasta_path, set([gene_id]), quiet=True)

            # pdb_filepath is the file path for the best structure model
            pdb_filepath = self.run_modeller()

            # annotation_subset_for_gene_id is a dataframe that stores annotations made by all
            # annotation methods (e.g.  DSSP) for the current gene_id. Each time an annotation
            # source is ran, its results are appended as columns to annotation_subset_for_gene_id.
            # All annotation sources must have the index called "residue_index" whose values are
            # anvi'o-indexed, i.e. the methionine has index 0. Each annotation source does NOT have
            # to annotate each residue in the gene. Each annotation source should have distinct
            # column names from all other annotation sources, however if they do not they will be
            # generically uniquified with the internal pandas method `_maybe_dep_names`, e.g. 2
            # columns named "A" will be named "A" and "A.1".
            if self.list_of_residue_annotation_methods:
                annotation_subset_for_gene_id = pd.DataFrame({})

                for annotation_method in self.list_of_residue_annotation_methods:
                    annotation_subset_for_gene_id = pd.concat([annotation_subset_for_gene_id, annotation_method(gene_id, pdb_filepath)], axis=1)
                annotation_subset_for_gene_id.columns = pd.io.parsers.ParserBase({'names':annotation_subset_for_gene_id.columns})._maybe_dedup_names(annotation_subset_for_gene_id.columns)
                annotation_subset_for_gene_id.insert(0, "gene_callers_id", gene_id)

                self.residue_annotation_df = self.residue_annotation_df.append(annotation_subset_for_gene_id)

            self.append_results_to_structure_db(gene_id, pdb_filepath)

            if self.full_output:
                self.dump_results_to_full_output()

        self.structure_db.store()
        self.structure_db.close()


    def append_results_to_structure_db(self, gene_id, pdb_filepath):
        self.structure_db.append(gene_id, pdb_filepath)


    def dump_results_to_full_output(self):
        """
        if self.full_output, all files from MODELLERs temp directory are recursively moved into
        output_gene_dir. Otherwise, the list of files we care about are defined in this function
        and moved into output_gene_dir.
        """
        output_gene_dir = os.path.join(self.full_output, self.modeller.gene_id)
        filesnpaths.check_output_directory(output_gene_dir)
        shutil.move(self.modeller.directory, output_gene_dir)

    def run_DSSP(self, gene_id, pdb_filepath):
        # Determine the model name by loading the structure file
        p = PDBParser()
        structure = p.get_structure(gene_id, pdb_filepath)
        model = structure[0] # pdb files can have multiple models. DSSP assumes the first.

        # run DSSP
        residue_annotation = DSSP(model, pdb_filepath, dssp = self.DSSP_executable, acc_array = "Wilke")

        # convert to a digestible format
        return self.convert_DSSP_output_from_biopython_to_dataframe(residue_annotation)


    def convert_DSSP_output_from_biopython_to_dataframe(self, dssp_biopython_object):
        """
        From the DSSP module in Biopython:
            ============ ===
            Tuple Index  Value
            ============ ===
            0            DSSP index
            1            Amino acid
            2            Secondary structure
            3            Relative ASA
            4            Phi
            5            Psi
            6            NH-->O_1_relidx
            7            NH-->O_1_energy
            8            O-->NH_1_relidx
            9            O-->NH_1_energy
            10           NH-->O_2_relidx
            11           NH-->O_2_energy
            12           O-->NH_2_relidx
            13           O-->NH_2_energy
            ============ ===

        Changes from Biopython format to anvi'o format:
            - residue index converted from 1Met to 0Met
            - aa converted to 3-letter code
            - ss type "-" is converted to coil (C)
            - relative indicies for h-bonds replaced with absolute residue indices
              (e.g. if relative index = -1 for residue 4, the absolute residue index is 3)
        """

        one_to_three = {v: k for k, v in constants.AA_to_single_letter_code.items()}
        columns = ("residue_index",
                   "aa",
                   "sec_struct",
                   "rel_solvent_acc",
                   "phi",
                   "psi",
                   "NH-O_1_index",
                   "NH-O_1_energy",
                   "O-NH_1_index",
                   "O-NH_1_energy",
                   "NH-O_2_index",
                   "NH-O_2_energy",
                   "O-NH_2_index",
                   "O-NH_2_energy")

        # convert biopython object to dictionary d
        d = {}
        for key in dssp_biopython_object.keys():
            d[key] = list(dssp_biopython_object[key])
            d[key][columns.index("residue_index")] = utils.convert_sequence_indexing(d[key][columns.index("residue_index")], source="not anvio", destination="anvio")
            d[key][columns.index("aa")] = one_to_three[d[key][columns.index("aa")]]
            if d[key][columns.index("sec_struct")] == "-":
                d[key][columns.index("sec_struct")] = "C"

            for hbond in ["NH-O_1", "O-NH_1", "NH-O_2", "O-NH_2"]:
                res_index = d[key][columns.index("residue_index")]
                rel_index = d[key][columns.index(hbond+"_index")]
                if rel_index == 0:
                    d[key][columns.index(hbond+"_index")] = np.nan
                    d[key][columns.index(hbond+"_energy")] = np.nan
                else:
                    d[key][columns.index(hbond+"_index")] = res_index + rel_index

        # convert dictionary d to dataframe df
        return pd.DataFrame(d, index=columns).T.set_index("residue_index")


    def run_modeller(self):
        self.modeller = MODELLER.MODELLER(self.args, run=self.run, progress=self.progress)
        return self.modeller.get_best_model()



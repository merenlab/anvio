# -*- coding: utf-8
# pylint: disable=line-too-long

"""Classes to make sense of genes and variability within the context of protein structure"""

import os
import time
import numpy as np
import shutil
import pandas as pd
import sqlite3
import warnings
import datetime

# multiprocess is a fork of multiprocessing that uses the dill serializer instead of pickle
# using the multiprocessing module directly results in a pickling error in Python 3.10 which
# goes like this:
#
#   >>> AttributeError: Can't pickle local object 'SOMEFUNCTION.<locals>.<lambda>' multiprocessing
#
import multiprocess as multiprocessing

from Bio.PDB import DSSP, PDBParser

import anvio
import anvio.db as db
import anvio.utils as utils
import anvio.dbops as dbops
import anvio.tables as t
import anvio.fastalib as u
import anvio.terminal as terminal
import anvio.constants as constants
import anvio.filesnpaths as filesnpaths
import anvio.drivers.MODELLER as MODELLER

from anvio.errors import ConfigError, FilesNPathsError
from anvio.dbops import ContigsSuperclass

J = lambda x, y: os.path.join(x, y)

from Bio.PDB.PDBExceptions import PDBConstructionWarning
warnings.simplefilter(action='ignore', category=PDBConstructionWarning)


class StructureDatabase(object):
    """Structure database operations"""

    def __init__(self, file_path, db_hash=None, create_new=False, ignore_hash=False, run=terminal.Run(), progress=terminal.Progress(), quiet=False):
        self.db_type = 'structure'
        self.db_hash = str(db_hash)
        self.version = anvio.__structure__version__
        self.file_path = file_path
        self.quiet = quiet
        self.run = run
        self.progress = progress
        self.create_new = create_new

        if not db_hash and create_new:
            raise ConfigError("You cannot create a Structure DB without supplying a DB hash.")

        self.db = db.DB(self.file_path, self.version, new_database = create_new)

        if create_new:
            self.create_tables()
            self.genes_with_structure = set([])
            self.genes_without_structure = set([])
            self.genes_queried = set([])
        else:
            self.db_hash = str(self.db.get_meta_value('contigs_db_hash'))
            self.genes_with_structure = self.get_genes_with_structure()
            self.genes_without_structure = self.get_genes_without_structure()
            self.genes_queried = self.get_genes_queried()

            if not len(self.genes_queried):
                raise ConfigError("Interesting... This structure database has no gene caller ids. Anvi'o is"
                                  "not sure how you managed that. Please send a report to the "
                                  "developers. Thank you.")

        if not ignore_hash:
            self.check_hash()

        # entries initialized as empty list are added with insert_many()
        # entries initialized as empty DataFrame are added with insert_rows_from_dataframe()
        self.entries = {
            t.pdb_data_table_name : [],
            t.residue_info_table_name : pd.DataFrame({}),
            t.templates_table_name : pd.DataFrame({}),
            t.models_table_name : pd.DataFrame({}),
        }


    def get_run_params_dict(self):
        """Return dictionary containing all pertinent run parameters that were used during the creation of this DB"""

        run_params_dict = {}

        run_params_dict['modeller_database'] = self.db.get_meta_value('modeller_database', try_as_type_int=False)
        run_params_dict['scoring_method'] = self.db.get_meta_value('scoring_method', try_as_type_int=False)
        run_params_dict['percent_cutoff'] = float(self.db.get_meta_value('percent_cutoff', try_as_type_int=False))
        run_params_dict['alignment_fraction_cutoff'] = float(self.db.get_meta_value('alignment_fraction_cutoff', try_as_type_int=False))
        run_params_dict['very_fast'] = bool(self.db.get_meta_value('very_fast', try_as_type_int=True))
        run_params_dict['deviation'] = float(self.db.get_meta_value('deviation', try_as_type_int=False))
        run_params_dict['max_number_templates'] = self.db.get_meta_value('max_number_templates', try_as_type_int=True)
        run_params_dict['num_models'] = self.db.get_meta_value('num_models', try_as_type_int=True)
        run_params_dict['skip_DSSP'] = bool(self.db.get_meta_value('skip_DSSP', try_as_type_int=True))

        return run_params_dict


    def get_genes_with_structure(self):
        """Returns set of gene caller ids that have a structure in the DB, queried from self table"""

        return set([int(x) for x in self.db.get_meta_value('genes_with_structure', try_as_type_int=False).split(',') if not x == ''])


    def get_genes_without_structure(self):
        """Returns set of gene caller ids that failed to generate a structure, queried from self table"""

        return set([int(x) for x in self.db.get_meta_value('genes_without_structure', try_as_type_int=False).split(',') if not x == ''])


    def get_genes_queried(self):
        """Returns set of all gene caller ids that structures were attempted for, regardless of success or failure

        Queried from self table

        Notes
        =====
        - FIXME This inefficiency could pose problems in 2030 when we have structure dbs with
          millions of structures
        """

        return self.get_genes_with_structure() | self.get_genes_without_structure()


    def update_genes_with_and_without_structure(self):
        """Writes genes_queried, genes_with_structure, and genes_without_structure entries in self table"""

        self.db.set_meta_value('genes_queried', ",".join(str(g) for g in self.genes_queried))
        self.db.set_meta_value('genes_with_structure', ",".join(str(g) for g in self.genes_with_structure))
        self.db.set_meta_value('genes_without_structure', ",".join(str(g) for g in self.genes_without_structure))


    def create_tables(self):
        self.db.set_meta_value('db_type', self.db_type)
        self.db.set_meta_value('contigs_db_hash', self.db_hash)
        self.db.set_meta_value('creation_date', time.time())

        self.db.create_table(t.pdb_data_table_name, t.pdb_data_table_structure, t.pdb_data_table_types)
        self.db.create_table(t.templates_table_name, t.templates_table_structure, t.templates_table_types)
        self.db.create_table(t.models_table_name, t.models_table_structure, t.models_table_types)
        self.db.create_table(t.residue_info_table_name, t.residue_info_table_structure, t.residue_info_table_types)
        self.db.create_table(t.states_table_name, t.states_table_structure, t.states_table_types)


    def store_modeller_params(self, modeller_params):
        """Store all run parameters in the self table

        Parameters
        ==========
        modeller_params : dict
            For e.g. {
                'modeller_database': 'pdb_95',
                'scoring_method': 'DOPE_score',
                'max_number_templates': 5,
                'percent_cutoff': 30,
                'alignment_fraction_cutoff': 80,
                'num_models': 1,
                'deviation': 4.0,
                'very_fast': True
            }
        """

        for param, value in modeller_params.items():
            self.db.set_meta_value(param, value)


    def check_hash(self):
        actual_db_hash = str(self.db.get_meta_value('contigs_db_hash'))
        if self.db_hash != actual_db_hash:
            raise ConfigError('The hash value inside Structure Database "%s" does not match with Contigs Database hash "%s",\
                               these files probably belong to different projects.' % (actual_db_hash, self.db_hash))


    def store(self, table_name):
        """Stores entries placed in self.entries[table_name], then empties self.entries[table_name]"""

        rows_data = self.entries[table_name]

        if type(rows_data) == list:
            self.db.insert_many(table_name, entries=rows_data)
            self.entries[table_name] = []

        elif type(rows_data) == pd.core.frame.DataFrame:
            self.db.insert_rows_from_dataframe(table_name, rows_data, raise_if_no_columns=False)
            self.entries[table_name] = pd.DataFrame({})

        else:
            raise ConfigError("store :: rows_data must be either a list of tuples or a pandas dataframe.")


    def remove_gene(self, corresponding_gene_call, remove_from_self=True):
        """Remove a gene from the structure database"""

        # Remove from tables
        self.db.remove_some_rows_from_table(
            t.residue_info_table_name,
            where_clause="corresponding_gene_call = %d" % corresponding_gene_call,
        )
        self.db.remove_some_rows_from_table(
            t.templates_table_name,
            where_clause="corresponding_gene_call = %d" % corresponding_gene_call,
        )
        self.db.remove_some_rows_from_table(
            t.models_table_name,
            where_clause="corresponding_gene_call = %d" % corresponding_gene_call,
        )

        if remove_from_self:
            # Remove from self entries
            self.genes_queried.remove(corresponding_gene_call)
            self.genes_with_structure.remove(corresponding_gene_call)
            self.update_genes_with_and_without_structure()


    def get_pdb_content(self, corresponding_gene_call):
        """Returns the file content (as a string) of a pdb for a given gene"""

        if not corresponding_gene_call in self.genes_with_structure:
            raise ConfigError('The gene caller id {} was not found in the structure database :('.format(corresponding_gene_call))

        return self.db.get_single_column_from_table(
            t.pdb_data_table_name,
            'pdb_content',
            where_clause="corresponding_gene_call = %d" % corresponding_gene_call,
        )[0].decode('utf-8')


    def export_pdb_content(self, corresponding_gene_call, filepath, ok_if_exists=False):
        """Export the pdb of a gene to a filepath"""

        if filesnpaths.is_output_file_writable(filepath, ok_if_exists=ok_if_exists):
            pdb_content = self.get_pdb_content(corresponding_gene_call)
            with open(filepath, 'w') as f:
                f.write(pdb_content)

        return filepath


    def export_pdbs(self, genes_of_interest, output_dir, ok_if_exists=False):
        """Exports the pdbs of a collection of genes to an output dir (calls self.export_pdb_content)"""

        filesnpaths.check_output_directory(output_dir, ok_if_exists=ok_if_exists)
        filesnpaths.gen_output_directory(output_dir)

        for i, gene in enumerate(genes_of_interest):
            file_path = os.path.join(output_dir, 'gene_{}.pdb'.format(gene))
            self.export_pdb_content(gene, file_path, ok_if_exists=ok_if_exists)

        self.run.info('PDB file output', output_dir)


    def get_structure(self, corresponding_gene_call):
        """Return a anvio.structureops.Structure object for a given gene"""

        return Structure(self.export_pdb_content(corresponding_gene_call, filesnpaths.get_temp_file_path()))


    def get_residue_info_for_gene(self, corresponding_gene_call, drop_null=True):
        """Get residue info for gene as a dataframe

        Parameters
        ==========
        drop_null : bool, True
            Drop columns that have all null values.
        """

        return self.db.get_table_as_dataframe(
            table_name=t.residue_info_table_name,
            where_clause="corresponding_gene_call = %d" % corresponding_gene_call,
            drop_if_null=drop_null,
        )


    def get_residue_info_for_all(self, drop_null=True):
        """Get the full residue info table as a dataframe

        Parameters
        ==========
        drop_null : bool, True
            Drop columns that have all null values.
        """

        return self.db.get_table_as_dataframe(
            table_name=t.residue_info_table_name,
            drop_if_null=drop_null,
        )


    def get_template_info_for_gene(self, corresponding_gene_call):
        """Get template info for gene as a dataframe"""

        return self.db.get_table_as_dataframe(
            table_name=t.templates_table_name,
            where_clause="corresponding_gene_call = %d" % corresponding_gene_call,
        )


    def get_template_info_for_all(self):
        """Get the full template info as a dataframe"""

        return self.db.get_table_as_dataframe(
            table_name=t.templates_table_name,
        )


    def disconnect(self):
        self.db.disconnect()


class StructureSuperclass(object):
    """Structure operations

    Parameters
    ==========
    args : argparse.Namespace

    create : bool, False
        Whether or not the structure DB is going to be made or not. If False, it should already
        exist
    """

    def __init__(self, args, create=False, run=terminal.Run(), progress=terminal.Progress()):
        self.args = args
        self.create = create
        self.run = run
        self.progress = progress

        A = lambda x, t: t(args.__dict__[x]) if x in self.args.__dict__ else None
        null = lambda x: x

        self.contigs_db_path = A('contigs_db', null)
        self.structure_db_path = A('structure_db', null)
        self.external_structures_path = A('external_structures', null)
        self.modeller_executable = A('modeller_executable', null)
        self.list_modeller_params = A('list_modeller_params', null)
        self.full_modeller_output = A('dump_dir', null)

        self.run_mode = 'modeller' if not self.external_structures_path else 'external'

        self.num_threads = A('num_threads', int)
        self.queue_size = self.num_threads * 2
        self.write_buffer_size = self.num_threads * A('write_buffer_size_per_thread', int)

        self.genes_of_interest_path = A('genes_of_interest', null)
        self.gene_caller_ids = A('gene_caller_ids', null)
        self.rerun_genes = A('rerun_genes', null)

        utils.is_contigs_db(self.contigs_db_path)
        contigs_db = dbops.ContigsDatabase(self.contigs_db_path)
        self.contigs_db_hash = contigs_db.meta['contigs_db_hash']
        contigs_db.disconnect()

        # init ContigsSuperClass
        self.contigs_super = ContigsSuperclass(self.args, r=terminal.Run(verbose=False), p=terminal.Progress(verbose=False))

        if self.create:
            # check database output
            if not self.structure_db_path:
                self.structure_db_path = "STRUCTURE.db"
            if not self.structure_db_path.endswith('.db'):
                raise ConfigError("The structure database output file (`-o / --output`) must end with '.db'")
            if filesnpaths.is_file_exists(self.structure_db_path, dont_raise=True):
                raise ConfigError("This structure DB already exists. Anvi'o will not overwrite")

            filesnpaths.is_output_file_writable(self.structure_db_path)


        # init StructureDatabase
        self.structure_db = StructureDatabase(self.structure_db_path, self.contigs_db_hash, create_new=create)

        if self.list_modeller_params:
            params_dict = self.structure_db.get_run_params_dict()
            for param, value in params_dict.items():
                self.run.info(param, value)
            import sys; sys.exit()

        # Determine the modeller parameters and store in db
        # NOTE self.skip_DSSP is down here because get_modeller_params has the potential to
        # overwrite self.args.skip_DSSP if create=False
        self.modeller_params = self.get_modeller_params()
        self.skip_DSSP = A('skip_DSSP', bool)

        if self.create:
            self.structure_db.store_modeller_params(self.modeller_params)
            self.structure_db.db.set_meta_value('skip_DSSP', self.skip_DSSP)

        # init annotation sources
        self.dssp = DSSPClass(skip_sanity_check=self.skip_DSSP) # If we skip DSSP, skip sanity_check

        self.sanity_check()


    def get_modeller_params(self):
        """Parses self.args to return dictionary of modeller parameters

        Returns
        =======
        output : dict

        Notes
        =====
        - If self.create=False, parameters in self.args accessed by this function are first set via
          self.set_prior_modeller_params
        - If self.run_mode != 'modeller', all values are set to null values
        """

        if not self.create:
            # Populates attributes of self.args based on metavalues in database
            self.set_prior_modeller_params()

        A = lambda x, t: t(self.args.__dict__[x]) if x in self.args.__dict__ else None
        null = lambda x: x

        if self.run_mode == 'modeller':
            return {
                'modeller_database': A('modeller_database', null),
                'scoring_method': A('scoring_method', null),
                'max_number_templates': A('max_number_templates', null),
                'percent_cutoff': A('percent_cutoff', null),
                'alignment_fraction_cutoff': A('alignment_fraction_cutoff', null),
                'num_models': A('num_models', null),
                'deviation': A('deviation', null),
                'very_fast': A('very_fast', bool),
            }
        else:
            return {
                'modeller_database': None,
                'scoring_method': None,
                'max_number_templates': 0,
                'percent_cutoff': 0,
                'alignment_fraction_cutoff': 0,
                'num_models': 0,
                'deviation': 0,
                'very_fast': False,
            }


    def set_prior_modeller_params(self):
        """Add the previous run parameters used during database creation into self.args

        If the class is initiated with create=False, it means that the user is going to modify the
        currently existing database, probably by running modeller. To ensure consistency between
        parameters used during creation and during modification, the previously used run parameters
        are determined from the structure database here and then set as arguments to self.args.
        """

        run_params_dict = self.structure_db.get_run_params_dict()
        for param, value in run_params_dict.items():
            setattr(self.args, param, value)


    def sanity_check(self):
        if self.run_mode == 'modeller':
            self.run.warning("Anvi'o will use 'MODELLER' by Webb and Sali (DOI: 10.1002/cpbi.3) to model "
                             "protein structures. If you publish your findings, please do not forget to "
                             "properly credit their work.", lc='green', header="CITATION")

        # if self.percent_cutoff is < 25, you should be careful about accuracy of models
        if self.run_mode == 'modeller' and self.modeller_params['percent_cutoff'] < 25:
            self.run.warning("You selected a percent identical cutoff of {}%. Below 25%, you should pay close attention "
                             "to the quality of the proteins... Keep in mind random sequence are expected to share "
                             "around 10% identity.".format(self.modeller_params['percent_cutoff']))

        if self.skip_DSSP:
            self.run.warning("It was requested that amino acid residue annotation with DSSP be skipped. A bold move only "
                             "an expert could justify... Anvi'o's respect for you increases slightly.")

        if self.run_mode == 'modeller':
            # Perform a rather extensive check on whether the MODELLER executable is going to work. We
            # do this here so we can initiate MODELLER.MODELLER with lazy_init so it does not do this
            # check every time
            self.args.modeller_executable = MODELLER.check_MODELLER(self.modeller_executable)
            self.modeller_executable = self.args.modeller_executable
            self.run.info_single("Anvi'o found the MODELLER executable %s, so will use it" % self.modeller_executable, nl_after=1, nl_before=1, mc='green')

            # Check and populate modeller databases if required
            MODELLER.MODELLER(self.args, filesnpaths.get_temp_file_path(), check_db_only=True)
        elif self.run_mode == 'external':
            self.run.info_single("Anvi'o will attempt to generate a database using external structures", nl_after=1, nl_before=1, mc='green')

            if self.full_modeller_output:
                raise ConfigError("No sense providing a --dump-dir when --external-structures are provided.")

            self.external_structures = ExternalStructuresFile(path=self.external_structures_path, contigs_db_path=self.contigs_db_path)


    def get_genes_of_interest(self, genes_of_interest_path=None, gene_caller_ids=None, raise_if_none=False):
        """Nabs the genes of interest based on genes_of_interest_path, gene_caller_ids, and self.external_structures

        If no genes of interest are provided through either genes_of_interest_path,
        gene_caller_ids, or self.external_structures, all will be assumed

        Parameters
        ==========
        raise_if_none : bool, False
            If True, an error will be raised if genes_of_interest_path and gene_caller_ids are both None
        """

        genes_of_interest = None

        if self.run_mode == 'external':
            if genes_of_interest_path or gene_caller_ids:
                raise ConfigError("You can't provide a --gene-caller-ids or --genes-of-interest concurrently with --external-structures. "
                                  "If you are trying to create a database from a subset of structures in your external structures file, "
                                  "please instead create an external structures file containing only those structures and use that instead."
                                  "Sorry for the inconvenience.")

            genes_of_interest = self.external_structures.content['gene_callers_id']
            return genes_of_interest

        # identify the gene caller ids of all genes available
        genes_in_contigs_database = set(self.contigs_super.genes_in_contigs_dict.keys())

        if not genes_in_contigs_database:
            raise ConfigError("This contigs database does not contain any identified genes...")

        # settling genes of interest
        if genes_of_interest_path and gene_caller_ids:
            raise ConfigError("You can't provide a gene caller id from the command line, and a list of gene caller ids "
                              "as a file at the same time, obviously.")

        if gene_caller_ids:
            gene_caller_ids = set([x.strip() for x in gene_caller_ids.split(',')])

            genes_of_interest = []
            for gene in gene_caller_ids:
                try:
                    genes_of_interest.append(int(gene))
                except:
                    raise ConfigError("Anvi'o does not like your gene caller id '%s'..." % str(gene))

            genes_of_interest = set(genes_of_interest)

        elif genes_of_interest_path:
            filesnpaths.is_file_tab_delimited(genes_of_interest_path, expected_number_of_fields=1)

            try:
                genes_of_interest = set([int(s.strip()) for s in open(genes_of_interest_path).readlines()])
            except ValueError:
                raise ConfigError("Well. Anvi'o was working on your genes of interest ... and ... those gene IDs did not "
                                  "look like anvi'o gene caller ids :/ Anvi'o is now sad.")

        if not genes_of_interest:
            if raise_if_none:
                raise ConfigError("You gotta supply some genes of interest.")

            # no genes of interest are specified. Assuming all, which could be innumerable--raise warning
            genes_of_interest = genes_in_contigs_database
            self.run.warning("You did not specify any genes of interest, so anvi'o will assume all of them are of interest.")

        # Check for genes that do not appear in the contigs database
        bad_gene_caller_ids = [g for g in genes_of_interest if g not in genes_in_contigs_database]
        if bad_gene_caller_ids:
            raise ConfigError(("This gene caller id you provided is" if len(bad_gene_caller_ids) == 1 else \
                               "These gene caller ids you provided are") + " not known to this contigs database: {}.\
                               You have only 2 lives left. 2 more mistakes, and anvi'o will automatically uninstall \
                               itself. Yes, seriously :(".format(", ".join([str(x) for x in bad_gene_caller_ids])))

        # Finally, raise warning if number of genes is greater than 20 FIXME determine average time
        # per gene and describe here
        if len(genes_of_interest) > 20:
            self.run.warning("Modelling protein structures is no joke. The number of genes you want protein structures for is "
                             "{}, which is a lot (of time!). If its taking too long, consider using the --very-fast flag. "
                             "CTRL + C to cancel.".format(len(genes_of_interest)))

        return genes_of_interest


    @staticmethod
    def worker(self, available_index_queue, output_queue):
        while True:
            try:
                corresponding_gene_call = available_index_queue.get(True)
                structure_info = self.process_gene(corresponding_gene_call)
                output_queue.put(structure_info)
            except Exception as e:
                # This thread encountered an error. We send the error back to the main thread which
                # will terminate the job.
                output_queue.put(e)

        # Code never reaches here because worker is terminated by main thread
        return


    def _run(self):
        """Calls either run_multi_thread or run_single_thread"""

        self.run.warning('', header='General info', nl_after=0, lc='green')
        self.run.info('contigs_db', self.contigs_db_path)
        self.run.info('contigs_db_hash', self.contigs_db_hash)
        self.run.info('structure_db_path', self.structure_db_path)
        self.run.info('num_threads', self.num_threads)
        self.run.info('write_buffer_size', self.write_buffer_size)
        self.run.info('genes_of_interest', self.genes_of_interest_path)
        self.run.info('gene_caller_ids', self.gene_caller_ids)
        self.run.info('skip_DSSP', self.skip_DSSP)
        self.run.info('run_mode', self.run_mode)

        if self.run_mode == 'modeller':
            self.run.warning('', header='Modeller parameters', nl_after=0, lc='green')
            self.run.info('modeller_executable', self.modeller_executable)
            for param, value in self.modeller_params.items():
                self.run.info(param, value)
            self.run.info('dump_dir', self.full_modeller_output, nl_after=1)

        if not anvio.DEBUG:
            self.run.warning("Do you want live info about how the modelling procedure is going for "
                             "each gene? Then restart this process with the --debug flag", lc='yellow')

        if self.run_mode == 'modeller' and not self.full_modeller_output:
            self.run.warning("When this finishes, do you want a potentially massive folder that "
                             "contains a murder of unorganized data in volumes that far exceed what "
                             "you could possibly want? Perfect, then restart this process and "
                             "provide a --dump-dir", lc='yellow')

        # Are we creating a new database, or updating an old one?
        update_db = True if not self.create else False

        if not update_db:
            # If creating new db, its simple. Just define the genes of interest
            genes_of_interest = self.get_genes_of_interest(self.genes_of_interest_path, self.gene_caller_ids)
        else:
            # If updating a db, we a bit more to do
            genes_of_interest = self.get_genes_of_interest(self.genes_of_interest_path, self.gene_caller_ids, raise_if_none=True)

            if self.rerun_genes:
                self.run.warning("You've requested to re-run structural modelling if your genes of interest are already "
                                 "present in the database. Just so you know.")
            else:
                genes_already_in_db = []
                for gene in genes_of_interest:
                    if gene in self.structure_db.genes_with_structure:
                        genes_already_in_db.append(gene)

                if len(genes_already_in_db):
                    raise ConfigError("Of your %d genes of interest, %d already exist in the DB. If you want to rerun "
                                      "modelling for these genes, use the flag --rerun. Here are the first 5 such gene "
                                      "IDs that are in your DB already: %s." % (len(genes_of_interest),
                                                                                len(genes_already_in_db),
                                                                                genes_already_in_db[:5]))

        self.run_multi_thread(genes_of_interest) if self.num_threads > 1 else self.run_single_thread(genes_of_interest)


    def run_multi_thread(self, genes_of_interest):
        structure_infos = []

        manager = multiprocessing.Manager()
        available_index_queue = manager.Queue()
        output_queue = manager.Queue(self.queue_size)

        # Get the genes of interest
        num_genes = len(genes_of_interest)

        # put contig indices into the queue to be read from within the worker
        for g in genes_of_interest:
            available_index_queue.put(g)

        processes = []
        for _ in range(self.num_threads):
            processes.append(multiprocessing.Process(target=StructureSuperclass.worker, args=(self, available_index_queue, output_queue)))

        for proc in processes:
            proc.start()

        mem_tracker = terminal.TrackMemory(at_most_every=0.1)
        mem_usage, mem_diff = mem_tracker.start()

        num_with_structure = 0
        num_without_structure = 0
        num_tried = 0

        self.progress.new('Using %d threads' % self.num_threads, progress_total_items=num_genes)
        msg = self.msg() % (num_tried, num_genes, num_with_structure, mem_usage, mem_diff)
        self.progress.update(msg)

        while num_tried < num_genes:
            try:
                structure_info = output_queue.get()

                if isinstance(structure_info, Exception):
                    # If thread returns an exception, we raise it and kill the main thread.
                    raise structure_info

                corresponding_gene_call = structure_info['corresponding_gene_call']

                # Add it to the storage buffer
                structure_infos.append(structure_info)

                num_tried += 1
                self.structure_db.genes_queried.add(corresponding_gene_call)

                if structure_info['has_structure']:
                    num_with_structure += 1
                    self.structure_db.genes_with_structure.add(corresponding_gene_call)
                else:
                    num_without_structure += 1
                    self.structure_db.genes_without_structure.add(corresponding_gene_call)

                if mem_tracker.measure():
                    mem_usage = mem_tracker.get_last()
                    mem_diff = mem_tracker.get_last_diff()

                self.progress.increment(num_tried)
                msg = self.msg() % (num_tried, num_genes, num_with_structure, mem_usage, mem_diff)
                self.progress.update(msg)

                if self.write_buffer_size > 0 and num_with_structure % self.write_buffer_size == 0:
                    self.progress.update('%d/%d GENES PROCESSED | WRITING TO DB 💾 ...' % (num_tried, num_genes))

                    # Store tables
                    self.store_genes(structure_infos)

                    # Empty storage buffer
                    structure_infos = []

                    # Update self table
                    self.structure_db.update_genes_with_and_without_structure()

                    self.progress.update(msg)

            except KeyboardInterrupt:
                self.run.info_single("Anvi'o received SIGINT, terminating all processes...", nl_before=2)
                break

            except Exception as worker_error:
                # An exception was thrown in one of the profile workers. We kill all processes in this case
                self.progress.end()
                for proc in processes:
                    proc.terminate()
                raise worker_error

        for proc in processes:
            proc.terminate()

        self.progress.update('%d/%d GENES PROCESSED | WRITING TO DB 💾 ...' % (num_tried, num_genes))

        # Store tables
        self.store_genes(structure_infos)

        # Empty storage buffer
        structure_infos = []

        # Update self table
        self.structure_db.update_genes_with_and_without_structure()

        self.progress.end(timing_filepath='anvio.debug.timing.txt' if anvio.DEBUG else None)

        if not num_with_structure:
            raise ConfigError("Well this is really sad. No structures were modelled, so there is nothing to do. Bye :'(")

        self.structure_db.disconnect()

        self.run.info("Structure database", self.structure_db_path)


    def run_single_thread(self, genes_of_interest):
        num_genes = len(genes_of_interest)

        mem_tracker = terminal.TrackMemory(at_most_every=0.1)
        mem_usage, mem_diff = mem_tracker.start()

        num_with_structure = 0
        num_without_structure = 0
        num_tried = 0

        self.progress.new('Using 1 thread', progress_total_items=num_genes)
        msg = self.msg() % (num_tried, num_genes, num_with_structure, mem_usage, mem_diff)
        self.progress.update(msg)

        try:
            for corresponding_gene_call in genes_of_interest:
                structure_info = self.process_gene(corresponding_gene_call)
                num_tried += 1

                self.progress.increment(num_tried)
                msg = self.msg() % (num_tried, num_genes, num_with_structure, mem_usage, mem_diff)
                self.progress.update(msg)

                if mem_tracker.measure():
                    mem_usage = mem_tracker.get_last()
                    mem_diff = mem_tracker.get_last_diff()

                self.store_gene(structure_info)

                self.dump_raw_results(structure_info)

                self.structure_db.genes_queried.add(corresponding_gene_call)

                if structure_info['has_structure']:
                    num_with_structure += 1
                    self.structure_db.genes_with_structure.add(corresponding_gene_call)
                else:
                    num_without_structure += 1
                    self.structure_db.genes_without_structure.add(corresponding_gene_call)

                # We update self.table every gene because there is no GIL cost with single thread and it
                # allows user to CTRL+C and still have a valid DB
                self.structure_db.update_genes_with_and_without_structure()
        except KeyboardInterrupt:
            self.run.info_single("Anvi'o received SIGINT, terminating all processes...", nl_before=2)

        self.progress.end(timing_filepath='anvio.debug.timing.txt' if anvio.DEBUG else None)

        if not num_with_structure:
            raise ConfigError("Well this is really sad. No structures were modelled, so there is nothing to do. Bye :'(")

        self.structure_db.disconnect()

        self.run.info("Structure database", self.structure_db_path)


    def msg(self):
        return '│ PROCESSED: %d/%d │ STRUCTURES: %d │ MEMORY: 🧠  %s (%s) │'


    def process_gene(self, corresponding_gene_call):
        structure_info = {
            'corresponding_gene_call': corresponding_gene_call,
            'has_structure': False,
        }

        if self.run_mode == 'modeller':
            directory = filesnpaths.get_temp_directory_path()

            # Export sequence
            target_fasta_path = filesnpaths.get_temp_file_path()
            self.contigs_super.get_sequences_for_gene_callers_ids([corresponding_gene_call],
                                                                  output_file_path=target_fasta_path,
                                                                  report_aa_sequences=True,
                                                                  simple_headers=True)

            try:
                filesnpaths.is_file_fasta_formatted(target_fasta_path)
            except FilesNPathsError:
                self.run.warning("You wanted to model a structure for gene ID %d, but the exported FASTA file "
                                 "is not what anvi'o considers a FASTA formatted file. The reason why this "
                                 "occassionally happens has not been investigated, but if it is any consolation, "
                                 "it is not your fault. You may want to try again, and maybe it will work. Or "
                                 "maybe it will not. Regardless, at this time anvi'o cannot model the gene. "
                                 "Here is the temporary fasta file path: %s " % (corresponding_gene_call, target_fasta_path))
                return structure_info

            if self.skip_gene_if_not_clean(corresponding_gene_call, fasta_path=target_fasta_path):
                return structure_info

            # Model structure
            structure_info['results'] = self.run_modeller(target_fasta_path, directory)

        elif self.run_mode == 'external':

            structure = self.external_structures.get_structure(corresponding_gene_call)
            if self.skip_gene_if_not_clean(corresponding_gene_call, sequence=structure.get_sequence()):
                return structure_info

            structure_info['results'] = self.create_results_dict_for_external_structure(corresponding_gene_call)

        # Annotate residues
        if structure_info['results']['structure_exists']:
            structure_info['has_structure'] = True

            structure_info['residue_info'] = self.get_gene_contribution_to_residue_info_table(
                corresponding_gene_call=corresponding_gene_call,
                pdb_filepath=structure_info['results']['best_model_path'],
            )

        return structure_info


    def create_results_dict_for_external_structure(self, corresponding_gene_call):
        return {
            'templates': {'pdb_id': ['none'], 'chain_id': ['none'], 'proper_percent_similarity': [0], 'percent_similarity': [0], 'align_fraction': [0]},
            'models': {'molpdf': [0], 'GA341_score': [0], 'DOPE_score': [0], 'picked_as_best': [True]},
            'corresponding_gene_call': corresponding_gene_call,
            'structure_exists': True,
            'best_model_path': self.external_structures.get_path(corresponding_gene_call),
            'best_score': None,
            'scoring_method': self.modeller_params['scoring_method'],
            'percent_cutoff': self.modeller_params['percent_cutoff'],
            'alignment_fraction_cutoff': self.modeller_params['alignment_fraction_cutoff'],
            'very_fast': self.modeller_params['very_fast'],
            'deviation': self.modeller_params['deviation'],
        }


    def skip_gene_if_not_clean(self, corresponding_gene_call, fasta_path=None, sequence=None):
        """Do not try modelling gene if it is not clean

        Parameters
        ==========
        corresponding_gene_call : int
            What is the gene callers id?
        fasta_path : str, None
            Provide either the path to the amino acid fasta
        sequence : str, None
            Or the amino acid sequence itself. Don't provide both
        """

        if not (fasta_path or sequence):
            raise ConfigError("skip_gene_if_not_clean :: provide a fasta_path or sequence")
        if fasta_path and sequence:
            raise ConfigError("skip_gene_if_not_clean :: don't provide both a fasta_path and a sequence")

        if fasta_path:
            fasta = u.SequenceSource(fasta_path); next(fasta)
            sequence = fasta.seq

        try:
            utils.is_gene_sequence_clean(sequence, amino_acid=True, can_end_with_stop=False, must_start_with_met=False)
            return False
        except ConfigError as error:
            self.run.warning("You wanted to model a structure for gene ID %d, but it is not what anvi'o "
                             "considers a 'clean gene'. Anvi'o will move onto the next gene. Here is the "
                             "error that was raised: \"%s\"" % (corresponding_gene_call, error.e))
            return True


    def get_gene_contribution_to_residue_info_table(self, corresponding_gene_call, pdb_filepath):
        results = [
            self.dssp.run(pdb_filepath, dont_run=self.skip_DSSP, drop=['aa']),
            self.run_contact_map_annotation(pdb_filepath),
            self.run_residue_identity_annotation(corresponding_gene_call, pdb_filepath),
        ]
        residue_annotation_for_gene = pd.concat(results, axis=1, sort=True)

        # add corresponding_gene_call and codon_order_in_gene as 0th and 1st columns
        residue_annotation_for_gene.insert(0, "corresponding_gene_call", corresponding_gene_call)
        residue_annotation_for_gene.insert(1, "codon_order_in_gene", residue_annotation_for_gene.index)

        return residue_annotation_for_gene


    def run_residue_identity_annotation(self, corresponding_gene_call, pdb_filepath):
        """A small routine to return a data frame containing codon numbers, codons, and amino acids"""

        nt_sequence = self.contigs_super.get_sequences_for_gene_callers_ids([corresponding_gene_call])
        nt_sequence = nt_sequence[1][corresponding_gene_call]['sequence']

        seq_dict = {"codon_order_in_gene": [],
                    "codon_number":        [],
                    "codon":               [],
                    "amino_acid":          []}

        last_codon = nt_sequence[:-3]
        if last_codon in ['TAA', 'TAG', 'TGA']:
            gene_length_in_codons = len(nt_sequence)//3 - 1 # subtract 1 because it's the stop codon
        else:
            gene_length_in_codons = len(nt_sequence)//3

        for codon_order_in_gene in range(gene_length_in_codons):
            seq_dict["codon_order_in_gene"].append(codon_order_in_gene)
            seq_dict["codon_number"].append(codon_order_in_gene+1)
            seq_dict["codon"].append(nt_sequence[3*codon_order_in_gene:3*(codon_order_in_gene + 1)])
            seq_dict["amino_acid"].append(constants.codon_to_AA[nt_sequence[3*codon_order_in_gene:3*(codon_order_in_gene + 1)]])

        return pd.DataFrame(seq_dict).set_index("codon_order_in_gene")


    def run_contact_map_annotation(self, pdb_path):
        """Returns the contact map in compressed form with index 'codon_order_in_gene' and column 'contact_numbers'"""

        # Run ContactMap class
        contact_map_long_form = Structure(pdb_path).get_boolean_contact_map(compressed=True, c='number')

        # Customize for this class and return
        column_rename = {
            'codon_number': 'codon_order_in_gene',
            'contacts': 'contact_numbers',
        }

        contact_map_long_form['codon_number'] -= 1

        return contact_map_long_form.rename(columns=column_rename).set_index('codon_order_in_gene')


    def run_modeller(self, target_fasta_path, directory):
        """Calls and returns results of MODELLER.MODELLER driver"""

        return MODELLER.MODELLER(self.args, target_fasta_path, directory=directory, lazy_init=True, skip_warnings=True).process()


    def dump_raw_results(self, structure_info):
        """Dump all raw modeller output into output_gene_dir if self.full_modeller_output"""

        if not self.full_modeller_output:
            return

        if 'results' not in structure_info:
            return

        output_gene_dir = os.path.join(self.full_modeller_output, structure_info['results']['corresponding_gene_call'])
        shutil.move(structure_info['results']['directory'], output_gene_dir)


    def store_gene(self, structure_info):
        """Store a gene's info into the structure database"""

        if 'results' not in structure_info:
            # There is nothing to store
            return

        corresponding_gene_call = structure_info["corresponding_gene_call"]
        results = structure_info['results']

        # If the gene is present in the database, remove it first
        if corresponding_gene_call in self.structure_db.genes_with_structure:
            # We do not remove the gene from self because that is handled in _run
            self.structure_db.remove_gene(corresponding_gene_call, remove_from_self=False)

        # templates is always added, even when structure was not modelled
        templates = pd.DataFrame(results['templates'])
        templates.insert(0, 'corresponding_gene_call', corresponding_gene_call)
        self.structure_db.entries[t.templates_table_name] = \
            self.structure_db.entries[t.templates_table_name].append(templates)
        self.structure_db.store(t.templates_table_name)

        # entries that are only added if a structure was modelled
        if results['structure_exists']:

            # models
            models = pd.DataFrame(results['models'])
            models.insert(0, 'corresponding_gene_call', corresponding_gene_call)
            self.structure_db.entries[t.models_table_name] = \
                self.structure_db.entries[t.models_table_name].append(models)
            self.structure_db.store(t.models_table_name)

            # pdb file data
            pdb_file = open(results['best_model_path'], 'rb')
            pdb_contents = pdb_file.read()
            pdb_file.close()
            pdb_table_entry = (corresponding_gene_call, pdb_contents)
            self.structure_db.entries[t.pdb_data_table_name].append(pdb_table_entry)
            self.structure_db.store(t.pdb_data_table_name)

            # residue_info
            self.structure_db.entries[t.residue_info_table_name] = \
                self.structure_db.entries[t.residue_info_table_name].append(structure_info['residue_info'])
            self.structure_db.store(t.residue_info_table_name)



    def store_genes(self, structure_infos):
        """Store info for a list of structure_infos. Calls self.store_gene"""

        for structure_info in structure_infos:
            self.store_gene(structure_info)
            self.dump_raw_results(structure_info)


class DSSPClass(object):
    def __init__(self, executable=None, skip_sanity_check=False):
        self.executable = executable
        self.skip_sanity_check = skip_sanity_check

        self.fields = [
            'codon_order_in_gene',
            'aa',
            'sec_struct',
            'rel_solvent_acc',
            'phi',
            'psi',
            'NH_O_1_index',
            'NH_O_1_energy',
            'O_NH_1_index',
            'O_NH_1_energy',
            'NH_O_2_index',
            'NH_O_2_energy',
            'O_NH_2_index',
            'O_NH_2_energy',
        ]

        if not self.skip_sanity_check:
            self.set_executable()
            utils.is_program_exists(self.executable)
            self.is_executable_a_working_DSSP_program()


    def run(self, pdb_filepath, drop=[], dont_run=False):
        """Run DSSP through Biopython. Return a dataframe

        Parameters
        ==========
        pdb_filepath : str
            Path to the pdb file you wish to run on. Assumes first model in file

        drop : list, []
            Which of self.fields do you want to removed from the dataframe?

        dont_run : bool, False
            If True, DSSP is not actually run. Instead, an empty dataframe is populated
        """

        if dont_run:
            # We don't actually run DSSP, we just make a dataframe populated with null values
            return self.get_null_dataframe(pdb_filepath, drop=drop)

        # Determine the model name by loading the structure file
        p = PDBParser()
        structure = p.get_structure(None, pdb_filepath)
        model = structure[0] # pdb files can have multiple models. DSSP assumes the first.

        # run DSSP
        residue_annotation = DSSP(model, pdb_filepath, dssp=self.executable, acc_array='Wilke')

        # convert to a digestible format
        return self.convert_DSSP_output_from_biopython_to_dataframe(residue_annotation, drop=drop)


    def get_null_dataframe(self, pdb_filepath, drop=[]):
        """Return a correctly sized dataframe with all null values, and codon_order_in_gene as the index"""

        # Get number of residues from PDB file
        num_residues = len(PDBParser().get_structure(None, pdb_filepath)[0].child_list[0])

        # Fields that take on null values
        fields = [field for field in self.fields if field not in drop and field != 'codon_order_in_gene']

        # Dict to define which column gets which null value
        fields_to_null_value_dict = {
            field: {'integer': np.nan, 'text': '', 'real': np.nan}[typ]
            for field, typ in zip(t.residue_info_table_structure, t.residue_info_table_types)
            if field in fields
        }

        d = {field: [fields_to_null_value_dict[field]] * num_residues for field in fields}

        return pd.DataFrame(d)


    def set_executable(self):
        """Determine which DSSP executables exist and should be used. Set to self.executable"""

        if self.executable:
            utils.is_program_exists(self.executable, dont_raise=True)
            return

        # Determine what DSSP program should be used. Tries mkdssp and then dssp, and raises
        # error if neither are found. mkdssp is newer and preferred
        if utils.is_program_exists("mkdssp", dont_raise=True):
            self.executable = "mkdssp"
        elif utils.is_program_exists("dssp", dont_raise=True):
            self.executable = "dssp"
        else:
            raise ConfigError("'mkdssp' or 'dssp' must be installed on your system, but "
                              "neither seem to appear in your path :/ If you are certain you have either on your "
                              "system (for instance you can run either by typing 'mkdssp' or 'dssp' in your terminal "
                              "window), you may want to send a detailed bug report. If you want to install DSSP, "
                              "check out http://merenlab.org/2016/06/18/installing-third-party-software/#dssp. "
                              "If you want to skip secondary structure and solvent accessibility annotation, "
                              "provide the flag --skip-DSSP.")


    def is_executable_a_working_DSSP_program(self):
        test_input = os.path.join(os.path.dirname(anvio.__file__), 'tests/sandbox/mock_data_for_structure/STRUCTURES/0.pdb')

        p = PDBParser()
        test_structure = p.get_structure('test', test_input)
        test_model = test_structure[0] # pdb files can have multiple models. DSSP assumes the first.

        # run DSSP
        try:
            test_residue_annotation = DSSP(test_model, test_input, dssp = self.executable, acc_array = "Wilke")
        except Exception:
            raise ConfigError('Your executable of DSSP, `{}`, doesn\'t appear to be working. For information on how to test '
                              'that your version is working correctly, please visit '
                              'http://merenlab.org/2016/06/18/installing-third-party-software/#dssp'\
                               .format(self.executable))

        if not len(test_residue_annotation.keys()):
            raise ConfigError("Your executable of DSSP, `{}`, exists but didn't return any meaningful output. This "
                              "is a known issue with certain distributions of DSSP. For information on how to test "
                              "that your version is working correctly, please visit "
                              "http://merenlab.org/2016/06/18/installing-third-party-software/#dssp"\
                               .format(self.executable))

        try:
            self.convert_DSSP_output_from_biopython_to_dataframe(test_residue_annotation)
        except:
            import pickle
            with open('troubleshoot_DDSP_output.pickle', 'wb') as output:
                pickle.dump(test_residue_annotation, output, pickle.HIGHEST_PROTOCOL)
            raise ConfigError('Your executable of DSSP ran and produced an output, but anvi\'o wasn\'t able to correctly parse it. '
                              'This is probably our fault. In your working directory should exist a file named "troubleshoot_DDSP_output.pickle". '
                              'Please send this to an anvio developer so that we can help identify what went wrong.')


    def convert_DSSP_output_from_biopython_to_dataframe(self, dssp_biopython_object, drop=[]):
        """Convert the output of the DSSP module in Biopython to a dataframe

        From the DSSP module in Biopython:

            ============ ==================== ================
            Tuple Index  Biopython            Anvi'o
            ============ ==================== ================
            0            DSSP index           codon_order_in_gene
            1            Amino acid           aa
            2            Secondary structure  sec_struct
            3            Relative ASA         rel_solvent_acc
            4            Phi                  phi
            5            Psi                  psi
            6            NH__>O_1_relidx      NH_O_1_index
            7            NH__>O_1_energy      NH_O_1_energy
            8            O__>NH_1_relidx      O_NH_1_index
            9            O__>NH_1_energy      O_NH_1_energy
            10           NH__>O_2_relidx      NH_O_2_index
            11           NH__>O_2_energy      NH_O_2_energy
            12           O__>NH_2_relidx      O_NH_2_index
            13           O__>NH_2_energy      O_NH_2_energy
            ============ ==================== ================

        Changes from Biopython format to anvi'o format:
            - residue index converted from 1Met to 0Met
            - aa converted to 3-letter code
            - ss type "-" is converted to coil (C)
            - relative indicies for h-bonds replaced with absolute residue indices
              (e.g. if relative index = -1 for residue 4, the absolute residue index is 3)

        Parameters
        ==========
        dssp_biopython_object : output of Bio.PDB.DSSP

        drop : list, []
            drop these columns (i.e. ['sec_struct', 'psi'])

        Returns
        =======
        output : pandas DataFrame
        """

        one_to_three = {v: k for k, v in constants.AA_to_single_letter_code.items()}

        # convert biopython object to dictionary d
        d = {}
        for key in dssp_biopython_object.keys():
            d[key] = list(dssp_biopython_object[key])
            d[key][self.fields.index('codon_order_in_gene')] = utils.convert_sequence_indexing(d[key][self.fields.index('codon_order_in_gene')], source='M1', destination='M0')
            d[key][self.fields.index('aa')] = one_to_three[d[key][self.fields.index('aa')]]

            if d[key][self.fields.index('sec_struct')] == '-':
                d[key][self.fields.index('sec_struct')] = 'C'

            for hbond in ['NH_O_1', 'O_NH_1', 'NH_O_2', 'O_NH_2']:
                res_index = d[key][self.fields.index('codon_order_in_gene')]
                rel_index = d[key][self.fields.index(hbond + '_index')]

                if rel_index == 0:
                    d[key][self.fields.index(hbond + '_index')] = np.nan
                    d[key][self.fields.index(hbond + '_energy')] = np.nan

                else:
                    d[key][self.fields.index(hbond + '_index')] = res_index + rel_index

        # convert dictionary d to dataframe df
        return pd.DataFrame(d, index=self.fields).T.set_index('codon_order_in_gene').drop(drop, axis=1)


class PDBDatabase(object):
    """This class manages the PDB database generated and managed by anvi-setup-pdb-database

    This class does not model structures, and does not create 'structure databases'. It generates
    and manages an exhaustive repository of available PDB structures. The database therefore allows
    offline access to the entire 'globe' of available structures.

    The collection of structures this database is a moving target as more structures become solved,
    but our friends at the Sali lab continually update this collection. For more information, see
    https://salilab.org/modeller/supplemental.html and for the list of PDB codes in the current
    version as well as when it was last updated please see
    https://salilab.org/modeller/downloads/pdb_95.cod.gz

    Parameters
    ==========
    pdb_database_path : str, None
        If None, the default data/misc/PDB.db filepath is used

    reset : bool, False
        If True, the DB is totally wiped and recreated from scratch

    update : bool, True
        Check to see if the DB is up to date. If it is not, store the latest additions. This can
        also be used if the previous run was cancelled prematurely.
    """

    def __init__(self, args, run=terminal.Run(), progress=terminal.Progress()):
        self.args = args
        self.run = run
        self.progress = progress

        A = lambda x, t: t(args.__dict__[x]) if x in self.args.__dict__ else None
        null = lambda x: x

        self.reset = A('reset', null)
        self.update = A('update', null)
        self.db_path = A('pdb_database_path', null) or constants.default_pdb_database_path
        self.skip_modeller_update = A('skip_modeller_update', null)

        self.num_threads = A('num_threads', int) or 1
        self.queue_size = self.num_threads * 5
        self.buffer_size = self.num_threads * 10

        self.exists = os.path.exists(self.db_path)
        self.clusters = pd.DataFrame({})


    def process(self):
        """Main method"""

        if self.reset:
            self.reset_db()

        if (self.exists and self.update) or (not self.exists and not self.update):
            self.check_or_create_db()
            self.database_summary()
            self._run()
        elif self.exists and not self.update:
            self.run.warning("Your database already exists, so anvi'o doesn't know what you "
                             "want to do. If you want to update your already existing database "
                             "use --update. If you want to delete it and start fresh use --reset. "
                             "Anvi'o is now in the process of displaying some stats for your "
                             "currently existing database. Afterwards, it will exit.",
                             header = "NOTHING TO DO", lc='yellow')
            self.check_or_create_db()
            self.database_summary()
        elif not self.exists and self.update:
            raise ConfigError("You asked to update the database, but the database was not found.")
        else:
            raise ConfigError("HOW?!")


    def get_clusters_dataframe(self):
        self.progress.new('Downloading')
        self.progress.update('95% clusters')

        clusters_path = self.download_clusters()
        clusters = self.parse_clusters_file(clusters_path)

        self.progress.end()

        self.run.warning("", header="DOWNLOADED CLUSTERS INFO", lc='green')
        self.run.info('Total structures', clusters.shape[0])
        self.run.info('Number of clusters', clusters['cluster'].nunique(), nl_after=1)

        return clusters


    def update_modeller_database(self):
        """Update the search database used by MODELLER so it shares the same proteins as this database"""

        modeller_database_dir = J(os.path.dirname(anvio.__file__), 'data/misc/MODELLER/db')
        pir_db = J(modeller_database_dir, 'pdb_95.pir')
        bin_db = J(modeller_database_dir, 'pdb_95.bin')
        dmnd_db = J(modeller_database_dir, 'pdb_95.dmnd')

        if self.skip_modeller_update or not filesnpaths.is_output_file_writable(pir_db):
            return

        self.progress.new('Updating')
        self.progress.update('MODELLER\'s search DB')

        if filesnpaths.is_file_exists(pir_db, dont_raise=True):
            os.remove(pir_db)

        if filesnpaths.is_file_exists(bin_db, dont_raise=True):
            os.remove(bin_db)

        if filesnpaths.is_file_exists(dmnd_db, dont_raise=True):
            os.remove(dmnd_db)

        db_download_path = os.path.join(modeller_database_dir, "pdb_95.pir.gz")
        utils.download_file("https://salilab.org/modeller/downloads/pdb_95.pir.gz", db_download_path)
        utils.run_command(['gzip', '-d', db_download_path], log_file_path=filesnpaths.get_temp_file_path())

        self.progress.end()

        self.run.info('MODELLER search DB updated', db_download_path.rstrip('.gz'), nl_after=1)


    def _run(self):
        self.clusters = self.get_clusters_dataframe()
        self.update_modeller_database()

        self.run.info_single("Anvi'o will now download structures from the RSCB PDB server. Press CTRL+C once to cancel.", nl_after=1)

        self.progress.new('Setup')
        self.progress.update('setting up threads')

        failed = []
        structures = []

        manager = multiprocessing.Manager()
        available_index_queue = manager.Queue()
        output_queue = manager.Queue(self.queue_size)

        self.progress.update('Determining the set of pdb ids')
        # Consider only PDB ids that aren't already stored
        # NOTE We store as list so order is retained when debugging
        pdb_ids = self.get_representative_ids(self.clusters)
        already_stored = self.get_stored_structure_ids()
        pdb_ids = [pdb_id for pdb_id in pdb_ids if pdb_id not in already_stored]

        num_structures = len(pdb_ids)

        if num_structures == 0:
            raise self.run.warning("The database is up to date :)")

        # put contig indices into the queue to be read from within the worker
        for p in pdb_ids:
            available_index_queue.put(p)

        processes = []
        for _ in range(self.num_threads):
            processes.append(multiprocessing.Process(target=PDBDatabase.get_structure, args=(self, available_index_queue, output_queue)))

        for proc in processes:
            proc.start()

        done = 0
        db_size = self.size_of_database()

        self.progress.end()
        self.progress.new('Using %d threads' % self.num_threads, progress_total_items=num_structures)
        self.progress.update('%d / %d processed | total DB size %s' % (done, num_structures, db_size))

        while done < num_structures:
            try:
                structure = output_queue.get()
                structures.append(structure)

                if structure['failed']:
                    failed.append(structure['id'])

                done += 1

                self.progress.increment(done)
                self.progress.update('%d / %d processed | total DB size %s' % (done, num_structures, db_size))

                if len(structures) > self.buffer_size:
                    self.store_buffer(structures)
                    db_size = self.size_of_database()
                    structures = []

            except KeyboardInterrupt:
                self.run.warning("", header="!!! DON'T TOUCH ANYTHING !!!", nl_before=2, nl_after=0, lc='yellow')
                self.run.info_single("Ending upon user request. Please wait patiently for anvi'o to "
                                     "die gracefully, or else risk a corrupted DB", nl_before=0)
                break

        for proc in processes:
            proc.terminate()

        self.store_buffer(structures)

        if len(failed):
            self.run.warning("The following structures (%d in total) couldn't be downloaded for whatever reason. Anvi'o "
                             "suggests you run this again with the --update flag to see if you can download "
                             "them. If you can't don't worry. %s" % (len(failed), failed))

        self.progress.end()


    @staticmethod
    def get_structure(self, available_index_queue, output_queue):
        while True:
            pdb_id = available_index_queue.get(True)

            structure = {}
            structure['id'] = pdb_id
            structure['pdb_content'] = self.download_and_return_pdb_content(pdb_id)
            structure['failed'] = False if structure['pdb_content'] else True
            structure['cluster_info'] = self.get_cluster_info(pdb_id)

            output_queue.put(structure)

        # Never reaches here because process is killed by main thread
        return


    def get_cluster_info(self, representative_id):
        cluster_id = self.clusters.loc[self.clusters['id'] == representative_id, 'cluster'].iloc[0]

        return self.clusters.loc[self.clusters['cluster'] == cluster_id, :]


    def download_and_return_pdb_content(self, pdb_id):
        """Download and return binary pdb content. pdb_id is a 5-letter code (PDB code + chain ID)"""

        temp_path = filesnpaths.get_temp_file_path()

        four_letter_code = pdb_id[:4]
        chain_id = pdb_id[-1]

        path = utils.download_protein_structure(four_letter_code, output_path=temp_path, chain=chain_id, raise_if_fail=False)

        if path:
            with open(path, 'rb') as f:
                pdb_content = f.read()
        else:
            # The structure could not be downloaded. Interesting :\
            self.run.warning("%s, chain %s was not downloaded :\\" % (pdb_id[:4], pdb_id[-1]))
            pdb_content = None

        if filesnpaths.is_file_exists(temp_path, dont_raise=True):
            os.remove(temp_path)

        return pdb_content


    def store_buffer(self, structures):
        pdb_db = self.load_db(write=True)

        # Remove any already existing entries
        new_clusters = set([])
        new_structure_ids = set([])

        for structure in structures:
            new_clusters.add(structure['cluster_info']['cluster'].iloc[0])
            new_structure_ids.add(structure['id'])

        if len(new_clusters):
            pdb_db.remove_some_rows_from_table('clusters', 'cluster IN (%s)' % ','.join(['"' + x + '"' for x in new_clusters]))

        if len(new_structure_ids):
            pdb_db.remove_some_rows_from_table('structures', 'representative_id IN (%s)' % ','.join(['"' + x + '"' for x in new_structure_ids]))

        # Add the new entries
        for structure in structures:
            # store the pdb content
            pdb_db.insert('structures', (structure['id'], structure['pdb_content']))

            # store the cluster
            pdb_db.insert_rows_from_dataframe('clusters', structure['cluster_info'])

        pdb_db.set_meta_value('last_update', str(datetime.datetime.now()))
        pdb_db.disconnect()


    def parse_clusters_file(self, path):
        """Parse a clusters file and return a clusters dataframe

        Returns
        =======
        output : pandas DataFrame
            output has columns Index(['cluster', 'id', 'representative'])
        """

        d = {'cluster': [], 'id': [], 'representative': []}

        with open(path, 'r') as f:
            for cluster, line in enumerate(f.readlines()):
                line = line.strip()
                rep_id, ids = line.split(':')
                rep_id = rep_id.strip()
                ids = ids.strip().split(' ')

                for pdb_id in ids:
                    d['cluster'].append('cluster_%06d' % cluster)
                    d['id'].append(pdb_id)
                    d['representative'].append(1 if pdb_id == rep_id else 0)

        return pd.DataFrame(d)


    def download_clusters(self, directory=None):
        """Downloads the current 95% cluster groups

        Downloads https://salilab.org/modeller/downloads/pdb_95.grp.gz and then decompresses the
        file with `gzip -d`

        Parameters
        ==========
        directory : str
            Directory that the file will be placed in. Temporary directory assumed if None

        Returns
        =======
        output : str
            The path of the downloaded and decompressed file.
        """

        if not directory:
            directory = filesnpaths.get_temp_directory_path()

        path = os.path.join(directory, "pdb_95.grp.gz")
        utils.download_file("https://salilab.org/modeller/downloads/pdb_95.grp.gz", path)
        utils.run_command(['gzip', '-d', path], log_file_path=filesnpaths.get_temp_file_path())

        return path.rstrip('.gz')


    def load_db(self, write=False):
        """Load the PDB database

        Parameters
        ==========
        write : bool, False
            Should the database be opened with write permissions enabled?
        """

        return db.DB(self.db_path, '0.1', new_database=False, read_only=(not write))


    def check_or_create_db(self):
        """Check the DB. Complain if it sucks. Create it if it doesn't exist"""

        if self.exists:
            try:
                pdb_db = self.load_db()
            except:
                raise ConfigError("Anvi'o is not convinced %s was generated with anvi-setup-pdb-database" % self.db_path)

            table_names = pdb_db.get_table_names()

            if 'structures' not in table_names:
                pdb_db.disconnect()
                raise ConfigError("Anvi'o was expecting to find the table with the name 'structures' in %s" % self.db_path)

            if 'clusters' not in table_names:
                pdb_db.disconnect()
                raise ConfigError("Anvi'o was expecting to find the table with the name 'clusters' in %s" % self.db_path)


        if not self.exists:
            pdb_db = db.DB(self.db_path, '0.1', new_database=True)
            pdb_db.set_meta_value('last_update', str(datetime.datetime.now()))
            pdb_db.set_meta_value('clustered_at', '95%')
            pdb_db.create_table('structures', ['representative_id', 'pdb_content'], ['text', 'blob'])
            pdb_db.create_table('clusters', ['cluster', 'id', 'representative'], ['text', 'text', 'integer'])

        try:
            # Create an index for `representative_id` for fast lookup
            pdb_db._exec("CREATE INDEX chain_index ON structures (representative_id);")
        except sqlite3.OperationalError:
            # The index has already been set
            pass

        pdb_db.disconnect()


    def reset_db(self):
        if not self.exists:
            raise ConfigError("You asked to --reset the database, but the database doesn't exist")

        self.run.warning("You want to --reset the database (size %s), which seems potentially rash. "
                         "Anvi'o is putting you in time out for 20 seconds before deleting. Press "
                         "CTRL+C at any time during your time out to cancel the deletion." \
                            % self.size_of_database())
        time.sleep(20)
        os.remove(self.db_path)
        self.exists = False


    def database_summary(self):
        self.run.warning("", header="DATABASE INFO", lc='green')
        self.run.info('Previous database found', self.exists)
        self.run.info('Database path', self.db_path)
        self.run.info('DB size', self.size_of_database())

        pdb_db = self.load_db()
        self.run.info('Last updated', pdb_db.get_meta_value('last_update'))
        self.run.info('Number of structures', pdb_db.get_row_counts_from_table('structures'), nl_after=1)
        pdb_db.disconnect()


    def get_representative_ids(self, clusters):
        """Get representative IDs from a clusters dataframe"""

        return clusters.loc[clusters['representative'] == 1, 'id'].tolist()


    def get_stored_structure_ids(self):
        """Get structure IDs of those stored in DB"""

        self.stored_structure_ids = self.get_representative_ids(self.get_clusters())

        return self.stored_structure_ids


    def get_clusters(self):
        """Get 'clusters' table as a dataframe"""

        pdb_db = self.load_db()
        clusters = pdb_db.get_table_as_dataframe('clusters', error_if_no_data=False)
        pdb_db.disconnect()

        return clusters


    def export_pdb(self, pdb_id, output_path=None):
        if not output_path:
            output_path = filesnpaths.get_temp_file_path()

        pdb_db = self.load_db()

        as_bytes = pdb_db.get_some_rows_from_table_as_dict('structures', 'representative_id = "%s"' % pdb_id)[pdb_id]['pdb_content']

        if as_bytes is None:
            raise ConfigError("PDB ID %s was found in the database, but its content is empty. That's weird--sorry." % pdb_db)

        with open(output_path, 'wb') as f:
            f.write(as_bytes)

        pdb_db.disconnect()

        return output_path


    def size_of_database(self):
        return utils.human_readable_file_size(os.path.getsize(self.db_path))


class Structure(object):
    def __init__(self, pdb_path):
        """Object to handle the analysis of PDB files"""

        self.distances_methods_dict = {
            'CA': self.calc_CA_dist, # alpha-carbon distances
            'COM': self.calc_COM_dist, # side chain centroid (center-of-mass) distances
        }

        self.path = pdb_path
        self._load_pdb_file(self.path)


    def _load_pdb_file(self, pdb_path, name_id='structure', model_index=0, chain_index=0):
        """Loads exactly one chain from one model of a PDB"""

        p = PDBParser()
        structure_obj = p.get_structure(name_id, pdb_path)

        if not len(structure_obj.get_list()):
            raise FilesNPathsError(f"Biopython was unable to parse {pdb_path} as a PDB file...")

        model = structure_obj.get_list()[model_index]
        self.structure = model.get_list()[chain_index]


    def get_sequence(self):
        return ''.join([constants.AA_to_single_letter_code[res.resname.capitalize()] for res in self.structure.get_list()])


    def get_contact_map(self, distance_method='CA', compressed=False, c='order'):
        """Returns a contact map (pairwise distances in Angstroms)

        Parameters
        ==========
        distance_method : str, 'CA'
            Which distance method? 'CA' is for alpha carbon distances. See
            self.distances_methods_dict for options

        compressed : bool, False
            Converts to long-form format (3 columns instead of square matrix)

        c : str, 'order'
            Determines whether contacts are defined according to codon_order_in_gene (i.e. 1st Met
            is 0) or codon_number (i.e. 1st Met is 1). Choose 'order' for codon_order_in_gene and
            'number' for codon_number. Valid only if compressed == True.

        Returns
        =======
        output: NxN 2-dim array
            Where N is the # of AAs

        Notes
        =====
        - See also `get_boolean_contact_map`
        """

        contact_map = np.zeros((len(self.structure), len(self.structure)))
        for i, residue1 in enumerate(self.structure):
            for j, residue2 in enumerate(self.structure):
                if i > j:
                    contact_map[i, j] = contact_map[j, i]
                else:
                    contact_map[i, j] = self.get_residue_to_residue_distance(residue1, residue2, distance_method)

        if compressed:
            return self.get_compressed_representation(contact_map, c=c)
        else:
            return contact_map


    def get_boolean_contact_map(self, distance_method='CA', threshold=6, compressed=False, c='order'):
        """Returns a boolean contact map (1 for touching, 0 for not)

        Parameters
        ==========
        distance_method : str, 'CA'
            Which distance method? 'CA' is for alpha carbon distances. See
            self.distances_methods_dict for options

        threshold : int or float, 6
            Residues are considered touching if their distances are less than or equal to this
            amount. Units are in Angstroms

        compressed : bool, False
            Converts to long-form format (3 columns instead of square matrix)

        c : str, 'order'
            Determines whether contacts are defined according to codon_order_in_gene (i.e. 1st Met
            is 0) or codon_number (i.e. 1st Met is 1). Choose 'order' for codon_order_in_gene and
            'number' for codon_number. Valid only if compressed == True.

        Returns
        =======
        output: NxN 2-dim array
            Where N is the # of AAs

        Notes
        =====
        - See also `get_contact_map`
        """

        contact_map = self.get_contact_map(distance_method=distance_method, compressed=False)

        contact_map[contact_map <= threshold] = 1
        contact_map[contact_map >  threshold] = 0

        if compressed:
            return self.get_compressed_representation(contact_map, c=c)
        else:
            return contact_map


    def get_compressed_representation(self, contact_map, c='order'):
        """Converts contact map into condensed representation

        Parameters
        ==========
        c : str, 'order'
            Determines whether contacts are defined according to codon_order_in_gene (i.e. 1st Met
            is 0) or codon_number (i.e. 1st Met is 1). Choose 'order' for codon_order_in_gene and
            'number' for codon_number

        Returns
        =======
        output : pandas DataFrame
            A dataframe with 2 columns, 'codon_order_in_gene'/'codon_number' (see Parameters for
            which it will be), and 'contacts'
        """

        col_name = 'codon_order_in_gene' if c == 'order' else 'codon_number'

        contacts_dict = {
            col_name: [],
            'contacts': [],
        }

        for codon_order_in_gene in range(contact_map.shape[0]):
            contacts = np.add(np.where(contact_map[codon_order_in_gene, :] == 1)[0], 1)

            # logic that handles if user wants in terms of codon_order_in_gene or codon_number
            if c == 'order':
                contacts = utils.convert_sequence_indexing(contacts, source='M1', destination='M0')
                index = codon_order_in_gene
            else:
                index = utils.convert_sequence_indexing(codon_order_in_gene, source='M0', destination='M1')

            # residues are not contacts with themselves
            contacts = contacts[contacts != index]

            contacts_dict[col_name].append(index)
            contacts_dict['contacts'].append(",".join([str(x) for x in contacts]))

        return pd.DataFrame(contacts_dict)


    def get_residue_to_residue_distance(self, residue1, residue2, method='CA'):
        return self.distances_methods_dict[method](residue1, residue2)


    def get_residue(self, codon_order_in_gene):
        """Return the residue object for a given codon_order_in_gene aka the 0-indexed residue position"""

        return self.structure.child_list[codon_order_in_gene]


    def get_residue_center_of_mass(self, residue):
        """Return the CoM of the sidechain of a Residue object (alpha-carbon included)"""

        M = 0
        weighted_coords = np.zeros(3)

        backbone = {'C', 'N', 'O'}

        for name, atom in residue.child_dict.items():
            if name in backbone:
                continue

            M += atom.mass
            weighted_coords += atom.coord * atom.mass

        return weighted_coords/M


    def calc_CA_dist(self, residue1, residue2):
        """Returns the C-alpha distance between two residues"""

        return residue1["CA"] - residue2["CA"]


    def calc_COM_dist(self, residue1, residue2):
        """Returns the distance between the sidechain center of masses between two residues"""

        COM1 = self.get_residue_center_of_mass(residue1)
        COM2 = self.get_residue_center_of_mass(residue2)

        return np.sqrt(np.sum((COM1 - COM2)**2))


class ExternalStructuresFile(object):
    def __init__(self, path, contigs_db_path, lazy=False, p=terminal.Progress(), r=terminal.Run()):
        """Check the integrity of an external structures file and provide contents as the attribute self.content

        Parameters
        ==========
        contigs_db_path : str, None
            The path to the corresponding contigs database.
        lazy : bool, False
            If false, each structure file will be opened and the sequence therein will be explicitly compared to
            the amino acid of the gene callers id found in the contigs database. If True, only superficial checks
            will be carried out, like making sure the file is tab-delimited and that all files pointed to actually
            exist.
        """

        self.run, self.progress = r, p

        self.path = path
        self.contigs_db_path = contigs_db_path

        utils.is_contigs_db(self.contigs_db_path)
        filesnpaths.is_file_tab_delimited(self.path)

        self.content = pd.read_csv(self.path, sep='\t')

        self.is_header_ok()
        self.is_duplicates()
        self.is_gene_caller_ids_ok()
        self.is_files_exist()
        if not lazy:
            self.test_integrity()


    def get_structure(self, gene_callers_id):
        """Return Structure object for given gene callers id"""

        path = self.get_path(gene_callers_id)
        return Structure(path)


    def get_path(self, gene_callers_id):
        """Return Structure object for given gene callers id"""

        result = self.content.loc[self.content['gene_callers_id'] == gene_callers_id, 'path']
        if result.empty:
            raise ConfigError(f"Structure.get_path :: Can't find gene callers id '{gene_callers_id}'.")
        return result.iloc[0]


    def is_header_ok(self):
        headers_proper = ['gene_callers_id', 'path']
        with open(self.path, 'r') as input_file:
            headers = input_file.readline().strip().split('\t')
            missing_headers = [h for h in headers_proper if h not in headers]

            if len(headers) != 2:
                raise FilesNPathsError("Your external structures file does not contain the right number of columns :/ Here are "
                                       "what the header columns should be called, in this order: '%s'." % ', '.join(headers_proper))

            if len(missing_headers):
                raise FilesNPathsError("Your external structures file has the wrong headers. They should be: '%s', not '%s'." % (', '.join(headers_proper), ', '.join(headers)))

        return True


    def is_duplicates(self):
        counts = self.content['gene_callers_id'].value_counts()
        multiple = counts[counts > 1].index.tolist()
        if len(multiple):
            raise FilesNPathsError(f"Only one structure can be assigned to each gene callers id. But the following genes are present "
                                   f"multiple times in your external structures file: {multiple}")


    def is_files_exist(self):
        """Check that all files pointed to in the file actually exist"""
        for _, row in self.content.iterrows():
            gene_callers_id, path = row['gene_callers_id'], row['path']
            if not filesnpaths.is_file_exists(path, dont_raise=True):
                raise FilesNPathsError(f"This is kind of an issue. Your external structures file points to the following path: {path}. "
                                       f"for gene callers id {gene_callers_id}. Well that path is not a file :\\")

        return True


    def is_gene_caller_ids_ok(self):
        """Returns True if all gene_callers_ids in external structures file are in the contigs database and are coding"""

        contigs_db = db.DB(self.contigs_db_path, client_version=None, ignore_version=True)
        table = contigs_db.get_table_as_dataframe('gene_amino_acid_sequences')
        genes_in_contigs_db_with_aa_seqs = table.loc[table['sequence'] != '', 'gene_callers_id'].tolist()
        missing_in_contigs = [x for x in self.content['gene_callers_id'] if x not in genes_in_contigs_db_with_aa_seqs]

        if len(missing_in_contigs):
            raise ConfigError(f"Some gene caller ids in your external structures file are either missing from your contigs database "
                              f"or are non-coding (they have no corresponding amino acid sequence). This is a show stopper. "
                              f"Here are the gene caller ids: {missing_in_contigs}")


    def test_integrity(self):
        """Parse the sequence contents of each PDB and ensure it matches the sequences in the contigs database"""

        # Fetch the amino acid sequences found in contigs database
        contigs_db = db.DB(self.contigs_db_path, client_version=None, ignore_version=True)
        table = contigs_db.get_table_as_dataframe('gene_amino_acid_sequences')
        amino_acid_sequences = dict(zip(table['gene_callers_id'], table['sequence']))

        self.progress.new('External structures', progress_total_items=self.content.shape[0])
        self.progress.update('Testing personal integrity')

        for _, row in self.content.iterrows():
            gene_callers_id, path = row['gene_callers_id'], row['path']
            s = Structure(path)
            aa_seq_structure = s.get_sequence()
            aa_seq_contigs = amino_acid_sequences[gene_callers_id]

            if aa_seq_structure != aa_seq_contigs:
                self.progress.end()
                raise ConfigError(f"The sequence in the structure for gene callers id {gene_callers_id} ({path}) does not match the sequence "
                                  f"found for this gene in the contigs database. Here is the sequence found in the structure: {aa_seq_structure}. "
                                  f"And here is the sequence in the contigs database that anvi'o was expecting: {aa_seq_contigs}.")

            self.progress.increment()
            self.progress.update(self.progress.msg)

        self.progress.end()
        return True



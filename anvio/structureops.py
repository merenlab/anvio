# -*- coding: utf-8
# pylint: disable=line-too-long

"""Classes to make sense of genes and variability within the context of protein structure"""

import os
import time
import numpy as np
import shutil
import pandas as pd

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

from anvio.errors import ConfigError
from anvio.dbops import ContigsSuperclass


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
        else:
            self.db_hash = str(self.db.get_meta_value('contigs_db_hash'))
            self.genes_with_structure = self.get_genes_with_structure()
            self.genes_without_structure = self.get_genes_without_structure()
            self.genes_queried = self.get_genes_queried()

            if not len(self.genes_queried):
                raise ConfigError("Interesting... This structure database has no gene caller ids. Anvi'o is"
                                  "not sure how you managed that. please send a report to the "
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
        run_params_dict['percent_identical_cutoff'] = float(self.db.get_meta_value('percent_identical_cutoff', try_as_type_int=False))
        run_params_dict['very_fast'] = bool(self.db.get_meta_value('very_fast', try_as_type_int=True))
        run_params_dict['deviation'] = float(self.db.get_meta_value('deviation', try_as_type_int=False))
        run_params_dict['max_number_templates'] = self.db.get_meta_value('max_number_templates', try_as_type_int=True)
        run_params_dict['num_models'] = self.db.get_meta_value('num_models', try_as_type_int=True)
        run_params_dict['skip_DSSP'] = bool(self.db.get_meta_value('skip_DSSP', try_as_type_int=True))

        return run_params_dict


    def get_genes_with_structure(self):
        """Returns list of gene caller ids that have a structure in the DB"""

        return [int(x) for x in self.db.get_meta_value('genes_with_structure', try_as_type_int=False).split(',') if not x == '']


    def get_genes_without_structure(self):
        """Returns list of gene caller ids that failed to generate a structure"""

        return [int(x) for x in self.db.get_meta_value('genes_without_structure', try_as_type_int=False).split(',') if not x == '']


    def get_genes_queried(self):
        """Returns list of all gene caller ids that structures were attempted for, regardless of success or failure

        Notes
        =====
        - FIXME This inefficiency could pose problems in 2030 when we have structure dbs with
          millions of structures
        """

        return self.get_genes_with_structure() + self.get_genes_without_structure()


    def create_tables(self):
        self.db.set_meta_value('db_type', self.db_type)
        self.db.set_meta_value('contigs_db_hash', self.db_hash)
        self.db.set_meta_value('creation_date', time.time())

        self.db.create_table(t.pdb_data_table_name, t.pdb_data_table_structure, t.pdb_data_table_types)
        self.db.create_table(t.templates_table_name, t.templates_table_structure, t.templates_table_types)
        self.db.create_table(t.models_table_name, t.models_table_structure, t.models_table_types)
        self.db.create_table(t.residue_info_table_name, t.residue_info_table_structure, t.residue_info_table_types)
        self.db.create_table(t.states_table_name, t.states_table_structure, t.states_table_types)


    def check_hash(self):
        actual_db_hash = str(self.db.get_meta_value('contigs_db_hash'))
        if self.db_hash != actual_db_hash:
            raise ConfigError('The hash value inside Structure Database "%s" does not match with Contigs Database hash "%s",\
                               these files probably belong to different projects.' % (actual_db_hash, self.db_hash))


    def store(self, table_name, key=None):
        rows_data = self.entries[table_name]

        if type(rows_data) == list:
            self.db.insert_many(table_name, entries=rows_data)
            self.entries[table_name] = []

        elif type(rows_data) == pd.core.frame.DataFrame:
            self.db.insert_rows_from_dataframe(table_name, rows_data, raise_if_no_columns=False, key=key)
            self.entries[table_name] = pd.DataFrame({})

        else:
            raise ConfigError("store :: rows_data must be either a list of tuples or a pandas dataframe.")


    def get_pdb_content(self, corresponding_gene_call):
        """Returns the file content (as a string) of a pdb for a given gene"""

        if not corresponding_gene_call in self.genes_with_structure:
            raise ConfigError('The gene caller id {} was not found in the structure database :('.format(corresponding_gene_call))

        return self.db.get_single_column_from_table(t.pdb_data_table_name,
            'pdb_content', where_clause="corresponding_gene_call = %d" % corresponding_gene_call)[0].decode('utf-8')


    def export_pdb_content(self, corresponding_gene_call, filepath, ok_if_exists=False):
        """Export the pdb of a gene to a filepath"""

        if filesnpaths.is_output_file_writable(filepath, ok_if_exists=ok_if_exists):
            pdb_content = self.get_pdb_content(corresponding_gene_call)
            with open(filepath, 'w') as f:
                f.write(pdb_content)


    def export_pdbs(self, genes_of_interest, output_dir, ok_if_exists=False):
        """Exports the pdbs of a collection of genes to an output dir (calls self.export_pdb_content)"""

        filesnpaths.check_output_directory(output_dir, ok_if_exists=ok_if_exists)
        filesnpaths.gen_output_directory(output_dir)

        for i, gene in enumerate(genes_of_interest):
            file_path = os.path.join(output_dir, 'gene_{}.pdb'.format(gene))
            self.export_pdb_content(gene, file_path, ok_if_exists=ok_if_exists)

        self.run.info('PDB file output', output_dir)


    def get_summary_for_interactive(self, corresponding_gene_call):
        """A very specific use case function. FIXME Should be moved to interactive.StructureInteractive"""

        summary = {}
        summary['pdb_content'] = self.get_pdb_content(corresponding_gene_call)
        summary['residue_info'] = self.db.get_table_as_dataframe(t.residue_info_table_name,
            where_clause = "corresponding_gene_call = %d" % corresponding_gene_call).to_json(orient='index')

        return summary


    def disconnect(self):
        self.db.disconnect()


class Structure(object):
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
        self.structure_db_path = A('structure_db_path', null)
        self.modeller_executable = A('modeller_executable', null)
        self.full_modeller_output = A('dump_dir', null)

        self.genes_of_interest_path = A('genes_of_interest', null)
        self.gene_caller_ids = A('gene_caller_ids', null)
        self.splits_of_interest_path = A('splits_of_interest', null)
        self.bin_id = A('bin_id', null)
        self.collection_name = A('collection_name', null)

        utils.is_contigs_db(self.contigs_db_path)
        self.contigs_db = dbops.ContigsDatabase(self.contigs_db_path)
        self.contigs_db_hash = self.contigs_db.meta['contigs_db_hash']

        if self.create:
            # MODELLER params FIXME may not be needed
            self.modeller_database = A('modeller_databaseth ', null)
            self.scoring_method = A('scoring_method', null)
            self.max_number_templates = A('max_number_templates', null)
            self.percent_identical_cutoff = A('percent_identical_cutoff', null)
            self.num_models = A('num_models', null)
            self.deviation = A('deviation', null)
            self.very_fast = A('very_fast', bool)

            self.skip_DSSP = A('skip_DSSP', bool)

            # check database output
            if not self.structure_db_path:
                self.structure_db_path = "STRUCTURE.db"
            if not self.structure_db_path.endswith('.db'):
                raise ConfigError("The structure database output file (`-o / --output`) must end with '.db'")
            filesnpaths.is_output_file_writable(self.structure_db_path)

            # identify which genes user wants to model structures for
            self.genes_of_interest = self.get_genes_of_interest(self.genes_of_interest_path, self.gene_caller_ids)

        else:
            self.genes_of_interest = None
            self.set_modeller_params()

        self.sanity_check()

        # init StructureDatabase
        self.structure_db = StructureDatabase(
            self.structure_db_path,
            self.contigs_db_hash,
            create_new=create,
        )

        # init ContigsSuperClass
        self.contigs_super = ContigsSuperclass(self.args)

        # init annotation sources
        self.contactmap = ContactMap()

        if not self.skip_DSSP:
            self.dssp = DSSPClass()


    def sanity_check(self):
        if self.genes_of_interest:
            # Check for genes that do not appear in the contigs database
            bad_gene_caller_ids = [g for g in self.genes_of_interest if g not in self.genes_in_contigs_database]
            if bad_gene_caller_ids:
                raise ConfigError(("This gene caller id you provided is" if len(bad_gene_caller_ids) == 1 else \
                                   "These gene caller ids you provided are") + " not known to this contigs database: {}.\
                                   You have only 2 lives left. 2 more mistakes, and anvi'o will automatically uninstall \
                                   itself. Yes, seriously :(".format(", ".join([str(x) for x in bad_gene_caller_ids])))

            # Finally, raise warning if number of genes is greater than 20 FIXME determine average time
            # per gene and describe here
            if len(self.genes_of_interest) > 20:
                self.run.warning("Modelling protein structures is no joke. The number of genes you want protein structures for is "
                                 "{}, which is a lot (of time!). If its taking too long, consider using the --very-fast flag. "
                                 "CTRL + C to cancel.".format(len(self.genes_of_interest)))

        # if self.percent_identical_cutoff is < 25, you should be careful about accuracy of models
        if self.percent_identical_cutoff < 25:
            self.run.warning("You selected a percent identical cutoff of {}%. Below 25%, you should pay close attention "
                             "to the quality of the proteins...".format(self.percent_identical_cutoff))

        if self.skip_DSSP:
            self.run.warning("You requested to skip amino acid residue annotation with DSSP. A bold move only an expert could justify... "
                             "Anvi'o's respect for you increases slightly.")

        # Perform a rather extensive check on whether the MODELLER executable is going to work. We
        # do this here so we can initiate MODELLER.MODELLER with lazy_init so it does not do this
        # check every time
        self.args.modeller_executable = MODELLER.check_MODELLER(self.modeller_executable)
        self.modeller_executable = self.args.modeller_executable
        self.run.info_single("Anvi'o found the MODELLER executable %s, so will use it" % self.modeller_executable, nl_after=1, mc='green')


    def get_genes_of_interest(self, genes_of_interest_path=None, gene_caller_ids=None):
        """Nabs the genes of interest based on user arguments (self.args)"""

        genes_of_interest = None

        # identify the gene caller ids of all genes available
        self.genes_in_contigs_database = set(dbops.ContigsSuperclass(self.args).genes_in_splits.keys())

        if not self.genes_in_contigs_database:
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
            # no genes of interest are specified. Assuming all, which could be innumerable--raise warning
            genes_of_interest = self.genes_in_contigs_database
            self.run.warning("You did not specify any genes of interest, so anvi'o will assume all of them are of interest.")

        return genes_of_interest


    def skip_gene_if_not_clean(self, corresponding_gene_call, fasta_path):
        """Do not try modelling gene if it is not clean"""

        fasta = u.SequenceSource(fasta_path); next(fasta)
        try:
            utils.is_gene_sequence_clean(fasta.seq, amino_acid=True, can_end_with_stop=False, must_start_with_met=False)
            return False
        except ConfigError as error:
            self.run.warning("You wanted to model a structure for gene ID %d, but it is not what anvi'o "
                             "considers a 'clean gene'. Anvi'o will move onto the next gene. Here is the "
                             "error that was raised: \"%s\"" % (corresponding_gene_call, error.e))
            return True


    def process(self):
        # which genes had structures and which did not. this information is added to the structure database self table
        has_structure = {True: [], False: []}

        num_genes_tried = 0
        num_genes_to_try = len(self.genes_of_interest)

        for corresponding_gene_call in self.genes_of_interest:
            # MODELLER outputs a lot of stuff into its working directory. A temporary directory is
            # made for each instance of MODELLER (i.e. each protein), And bits and pieces of this
            # directory are used in the creation of the structure database. If self.full_modeller_output is
            # provided, these directories and their contents are moved into self.full_modeller_output.
            self.args.directory = filesnpaths.get_temp_directory_path()
            self.args.target_fasta_path = filesnpaths.get_temp_file_path()

            # Export sequence
            dbops.export_aa_sequences_from_contigs_db(self.contigs_db_path,
                                                      self.args.target_fasta_path,
                                                      set([corresponding_gene_call]),
                                                      quiet = True)

            if self.skip_gene_if_not_clean(corresponding_gene_call, self.args.target_fasta_path):
                num_genes_tried += 1
                has_structure[False].append(corresponding_gene_call)
                continue

            # Model structure
            modeller_out = self.run_modeller()
            if modeller_out["structure_exists"]:
                self.run.info_single("Gene successfully modelled!", nl_after=1, mc="green")

            has_structure[modeller_out["structure_exists"]].append(str(corresponding_gene_call))

            # Annotate residues
            residue_info_dataframe = None
            if modeller_out["structure_exists"]:
                residue_info_dataframe = self.get_residue_info_table_for_gene(corresponding_gene_call,
                                                                              modeller_out["best_model_path"])
            self.store_gene(modeller_out, residue_info_dataframe)

            self.dump_raw_results()

            num_genes_tried += 1

        if not has_structure[True]:
            raise ConfigError("Well this is really sad. No structures were modelled, so there is nothing to do. Bye :'(")

        self.structure_db.db.set_meta_value('genes_queried', ",".join([str(g) for g in self.genes_of_interest]))
        self.structure_db.db.set_meta_value('genes_with_structure', ",".join(has_structure[True]))
        self.structure_db.db.set_meta_value('genes_without_structure', ",".join(has_structure[False]))
        self.structure_db.db.set_meta_value('modeller_database', self.modeller.modeller_database)
        self.structure_db.db.set_meta_value('scoring_method', self.scoring_method)
        self.structure_db.db.set_meta_value('percent_identical_cutoff', str(self.percent_identical_cutoff))
        self.structure_db.db.set_meta_value('very_fast', str(int(self.very_fast)))
        self.structure_db.db.set_meta_value('deviation', self.deviation)
        self.structure_db.db.set_meta_value('max_number_templates', self.max_number_templates)
        self.structure_db.db.set_meta_value('num_models', self.num_models)
        self.structure_db.db.set_meta_value('skip_DSSP', self.skip_DSSP)

        self.structure_db.disconnect()
        self.run.info("Structure database", self.structure_db_path)


    def get_residue_info_table_for_gene(self, corresponding_gene_call, pdb_filepath):
        results = [
            self.dssp.run(pdb_filepath, id_for_protein=corresponding_gene_call),
            self.run_contact_map_annotation(pdb_filepath),
            self.run_residue_identity_annotation(corresponding_gene_call, pdb_filepath),
        ]
        residue_annotation_for_gene = pd.concat(results, axis=1, sort=True)

        # add corresponding_gene_call and codon_order_in_gene as 0th and 1st columns
        residue_annotation_for_gene.insert(0, "entry_id", list(range(residue_annotation_for_gene.shape[0])))
        residue_annotation_for_gene.insert(1, "corresponding_gene_call", corresponding_gene_call)
        residue_annotation_for_gene.insert(2, "codon_order_in_gene", residue_annotation_for_gene.index)

        return residue_annotation_for_gene


    def run_residue_identity_annotation(self, corresponding_gene_call, pdb_filepath):
        nt_sequence = self.contigs_super.get_sequences_for_gene_callers_ids([corresponding_gene_call])
        nt_sequence = nt_sequence[1][corresponding_gene_call]['sequence']

        seq_dict = {"codon_order_in_gene": [],
                    "codon_number":        [],
                    "codon":               [],
                    "amino_acid":          []}

        gene_length_in_codons = len(nt_sequence)//3 - 1 # subtract 1 because it's the stop codon
        for codon_order_in_gene in range(gene_length_in_codons):
            seq_dict["codon_order_in_gene"].append(codon_order_in_gene)
            seq_dict["codon_number"].append(codon_order_in_gene+1)
            seq_dict["codon"].append(nt_sequence[3*codon_order_in_gene:3*(codon_order_in_gene + 1)])
            seq_dict["amino_acid"].append(constants.codon_to_AA[nt_sequence[3*codon_order_in_gene:3*(codon_order_in_gene + 1)]])

        return pd.DataFrame(seq_dict).set_index("codon_order_in_gene")


    def run_contact_map_annotation(self, pdb_path):
        """Returns the contact map in compressed form with index 'codon_order_in_gene' and column 'contact_numbers'"""

        # Run ContactMap class
        contact_map = self.contactmap.get_boolean_contact_map(pdb_path)
        compressed_rep = self.contactmap.get_compressed_representation(contact_map, c='number')

        # Customize for this class and return
        column_rename = {
            'codon_number': 'codon_order_in_gene',
            'contacts': 'contact_numbers',
        }

        compressed_rep['codon_number'] -= 1

        return compressed_rep.rename(columns=column_rename).set_index('codon_order_in_gene')


    def run_modeller(self):
        self.modeller = MODELLER.MODELLER(self.args, lazy_init=True)
        modeller_out = self.modeller.process()

        return modeller_out


    def dump_raw_results(self):
        """Dump all raw modeller output into self.output_gene_dir if self.full_modeller_output"""

        if not self.full_modeller_output:
            return

        output_gene_dir = os.path.join(self.full_modeller_output, self.modeller.corresponding_gene_call)
        shutil.move(self.modeller.directory, output_gene_dir)


    def store_gene(self, modeller_out, residue_info_dataframe):
        """Store a gene's info into the structure database"""

        corresponding_gene_call = modeller_out["corresponding_gene_call"]

        # templates is always added, even when structure was not modelled
        templates = pd.DataFrame(modeller_out["templates"])
        templates.insert(0, "corresponding_gene_call", corresponding_gene_call)
        templates = templates.reset_index().rename(columns={"index": "entry_id"})
        self.structure_db.entries[t.templates_table_name] = \
            self.structure_db.entries[t.templates_table_name].append(templates)
        self.structure_db.store(t.templates_table_name, key="entry_id")

        # entries that are only added if a structure was modelled
        if modeller_out["structure_exists"]:

            # models
            models = pd.DataFrame(modeller_out["models"])
            models.insert(0, "corresponding_gene_call", corresponding_gene_call)
            models = models.reset_index().rename(columns={"index": "entry_id"})
            self.structure_db.entries[t.models_table_name] = \
                self.structure_db.entries[t.models_table_name].append(models)
            self.structure_db.store(t.models_table_name, key="entry_id")

            # pdb file data
            pdb_file = open(modeller_out["best_model_path"], 'rb')
            pdb_contents = pdb_file.read()
            pdb_file.close()
            pdb_table_entry = (corresponding_gene_call, pdb_contents)
            self.structure_db.entries[t.pdb_data_table_name].append(pdb_table_entry)
            self.structure_db.store(t.pdb_data_table_name)

            # residue_info
            self.structure_db.entries[t.residue_info_table_name] = \
                self.structure_db.entries[t.residue_info_table_name].append(residue_info_dataframe)
            self.structure_db.store(t.residue_info_table_name, key="entry_id")


class StructureUpdate(Structure):
    def __init__(self, args, run=terminal.Run(), progress=terminal.Progress()):
        self.args = args
        self.run = run
        self.progress = progress

        # initialize self.arg parameters
        A = lambda x, t: t(args.__dict__[x]) if x in self.args.__dict__ else None
        null = lambda x: x
        self.contigs_db_path = A('contigs_db', null)
        self.structure_db_path = A('structure_db', null)
        self.genes_to_remove = A('genes_to_remove', null)
        self.genes_to_remove_path = A('genes_to_remove_file', null)
        self.genes_to_add = A('genes_to_add', null)
        self.genes_to_add_path = A('genes_to_add_file', null)
        self.full_modeller_output = A('dump_dir', null)
        self.modeller_executable = A('modeller_executable', null)
        self.skip_genes_if_already_present = A('skip_genes_if_already_present', bool)
        self.DSSP_executable = None

        utils.is_contigs_db(self.contigs_db_path)
        self.contigs_db = dbops.ContigsDatabase(self.contigs_db_path)
        self.contigs_db_hash = self.contigs_db.meta['contigs_db_hash']
        self.structure_db_path = self.structure_db_path

        # init ContigsSuperClass
        self.contigs_super = ContigsSuperclass(self.args)

        if not any([self.genes_to_remove, self.genes_to_remove_path, self.genes_to_add, self.genes_to_add_path]):
            raise ConfigError("Please specify some genes to add or remove to your database.")

        if self.genes_to_remove and self.genes_to_remove_path:
            raise ConfigError("Provide either --genes-to-remove or --genes-to-remove-path. You provided both.")

        if self.genes_to_add and self.genes_to_add_path:
            raise ConfigError("Provide either --genes-to-add or --genes-to-add-path. You provided both.")

        if self.genes_to_remove or self.genes_to_remove_path:
            self.run.warning("Removing genes...", header="Updating %s" % self.structure_db_path, lc='green')
            self.load_structure_db()
            remove = self.parse_genes(self.genes_to_remove, self.genes_to_remove_path)
            self.remove_genes(remove)
            self.structure_db.disconnect()

        if self.genes_to_add or self.genes_to_add_path:
            self.run.warning("Adding genes...", header="Updating %s" % self.structure_db_path, lc='green')
            self.load_structure_db()
            self.add_genes()


    def load_structure_db(self):
        utils.is_structure_db(self.structure_db_path)
        self.structure_db = StructureDatabase(self.structure_db_path,
                                              self.contigs_db_hash,
                                              create_new=False)

    def add_genes(self):
        # identify which genes user wants to model structures for
        self.genes_of_interest = self.get_genes_of_interest(self.genes_to_add_path, self.genes_to_add)

        if self.skip_genes_if_already_present:
            redundant_gene_caller_ids = [g for g in self.genes_of_interest if g in self.structure_db.genes_queried]
            if redundant_gene_caller_ids:
                self.run.info("Redundant gene caller ids that will be skipped", ",".join([str(x) for x in redundant_gene_caller_ids]))
                self.genes_of_interest = [g for g in self.genes_of_interest if g not in redundant_gene_caller_ids]
                if not self.genes_of_interest:
                    raise ConfigError("Every gene you wanted to add is already in the database. Since you provided "
                                      "the --skip-genes-if-already-present flag, there is nothing to do :)")

        self.run.info("Gene caller ids to be added", ",".join([str(x) for x in self.genes_of_interest]))

        self.get_MODELLER_params_used_when_db_was_created()

        self.sanity_check_for_adding_genes()

        # residue annotation
        self.residue_annotation_sources_info = self.get_residue_annotation_sources_info()
        self.residue_annotation_df = pd.DataFrame({})

        if self.full_modeller_output:
            self.full_modeller_output = filesnpaths.check_output_directory(self.full_modeller_output, ok_if_exists=True)

        self.process()
        self.run.info_single("Anvi'o attempted to add the requested genes. The above log can inform you which were successful.", nl_after=1, nl_before=1)


    def remove_genes(self, remove):
        self.progress.new("Removing genes from structure database")

        bad_ids = [x for x in remove if x not in self.structure_db.genes_queried]
        if len(bad_ids):
            if len(bad_ids) == len(remove):
                self.run.warning("All of the gene caller IDs you asked to remove are missing from "
                                 "the structure database, so there's no genes to remove. Here they "
                                 "are: [{}]. Anvi'o's trust in you decreases significantly."\
                                      .format(",".join([str(x) for x in bad_ids])))
                self.progress.end()
                return

            self.run.warning("Some of the gene caller ids you asked to remove aren't in the "
                             "structure database. Here they are: [{}].".format(",".join([str(x) for x in bad_ids])))

        remove = set([x for x in remove if x not in bad_ids])

        if remove == set(self.structure_db.genes_queried):
            raise ConfigError("You want to remove every gene in your structure database. No.")

        self.run.info("Gene caller ids to be removed", ",".join([str(x) for x in remove]))

        # remove ids from the three meta-keys in which they can appear
        new_genes_queried = [x for x in self.structure_db.genes_queried if x not in remove]
        self.structure_db.db.update_meta_value('genes_queried', ",".join([str(x) for x in new_genes_queried]))

        new_genes_with_structure = [x for x in self.structure_db.genes_with_structure if x not in remove]
        self.structure_db.db.update_meta_value('genes_with_structure', ",".join([str(x) for x in new_genes_with_structure]))

        new_genes_without_structure = [x for x in self.structure_db.genes_without_structure if x not in remove]
        self.structure_db.db.update_meta_value('genes_without_structure', ",".join([str(x) for x in new_genes_without_structure]))

        # remove all rows of tables in which corrresponding_gene_call matches the ids to remove
        where_clause = 'corresponding_gene_call IN (%s)' % ','.join(['{}'.format(x) for x in remove])
        for table_name in self.structure_db.db.get_table_names():
            if 'corresponding_gene_call' in self.structure_db.db.get_table_structure(table_name):
                self.structure_db.db.remove_some_rows_from_table(table_name, where_clause)

        self.run.info_single("The requested genes have been successfully removed.", nl_after=1)
        self.progress.end()


    def parse_genes(self, comma_delimited_genes=None, genes_filepath=None):
        if comma_delimited_genes:
            gene_caller_ids = set([x.strip() for x in comma_delimited_genes.split(',')])
            genes = []
            for gene in gene_caller_ids:
                try:
                    genes.append(int(gene))
                except:
                    raise ConfigError("Anvi'o does not like your gene caller id '%s'..." % str(gene))

        elif genes_filepath:
            filesnpaths.is_file_tab_delimited(genes_filepath, expected_number_of_fields=1)

            try:
                genes = set([int(s.strip()) for s in open(genes_filepath).readlines()])
            except ValueError:
                raise ConfigError("Well. Anvi'o was working on your genes in `%s` ... and ... those gene IDs did not "
                                  "look like anvi'o gene caller ids :/ Anvi'o is now sad." % genes_filepath)

        return set(genes)


    def get_MODELLER_params_used_when_db_was_created(self):
        self.progress.new("Determining parameters used during structure database creation")

        meta_table_dict = self.structure_db.db.get_table_as_dict("self")
        modeller_params = [
            ('modeller_database', str),
            ('scoring_method', str),
            ('percent_identical_cutoff', float),
            ('very_fast', lambda x: bool(int(x))),
            ('deviation', float),
            ('max_number_templates', int),
            ('num_models', int),
            ('skip_DSSP', lambda x: bool(int(x))),
        ]

        self.run.info_single("Previous parameters used for creating the structure database", nl_before=1)
        for param, cast_type in modeller_params:
            setattr(self, param, cast_type(meta_table_dict[param]["value"])) # set to self
            setattr(self.args, param, cast_type(meta_table_dict[param]["value"])) # set to self.args (passed to MODELLER)
            self.run.info(param, getattr(self, param))

        self.progress.end()


    def sanity_check_for_adding_genes(self):

        # check for genes that do not appear in the contigs database
        bad_gene_caller_ids = [g for g in self.genes_of_interest if g not in self.genes_in_contigs_database]
        if bad_gene_caller_ids:
            raise ConfigError(("This gene caller id you" if len(bad_gene_caller_ids) == 1 else \
                               "These gene caller ids you") + " want to add to the structure database\
                               are not known to the contigs database: {}. You have only 2 lives\
                               left. 2 more mistakes, and anvi'o will automatically uninstall\
                               itself. Yes, seriously :(".format(",".join([str(x) for x in bad_gene_caller_ids])))

        # check for genes that do already appear in the structure database
        redundant_gene_caller_ids = [g for g in self.genes_of_interest if g in self.structure_db.genes_queried]
        if redundant_gene_caller_ids and not self.skip_genes_if_already_present:
            raise ConfigError(("This gene caller id you" if len(redundant_gene_caller_ids) == 1 else \
                               "These gene caller ids you") + " want to add to the structure database\
                               is already in the structure database: {}. If you want to re-do the\
                               modelling, then first remove it with --genes-to-remove or\
                               --genes-to-remove-file (you can do it in the same\
                               anvi-update-structure-database command).".\
                                   format(",".join([str(x) for x in redundant_gene_caller_ids])))

        # raise warning if number of genes is greater than 20
        if len(self.genes_of_interest) > 20:
            self.run.warning("Modelling protein structures is no joke. The number of genes you want "
                             "to append to the structure database is {}, which is a lot (of time!). "
                             "CTRL + C to cancel.".format(len(self.genes_of_interest)))

        if not self.skip_DSSP:
            if utils.is_program_exists("mkdssp", dont_raise=True): # mkdssp is newer and preferred
                self.DSSP_executable = "mkdssp"

            if not self.DSSP_executable:
                if utils.is_program_exists("dssp", dont_raise=True):
                    self.DSSP_executable = "dssp"
                else:
                    raise ConfigError("An anvi'o function needs 'mkdssp' or 'dssp' to be installed on your system, but "
                                      "neither seem to appear in your path :/ If you are certain you have either on your "
                                      "system (for instance you can run either by typing 'mkdssp' or 'dssp' in your terminal "
                                      "window), you may want to send a detailed bug report. If you want to install DSSP, "
                                      "check out http://merenlab.org/2016/06/18/installing-third-party-software/#dssp. "
                                      "If you want to skip secondary structure and solvent accessibility annotation, "
                                      "provide the flag --skip-DSSP.")

            self.run.info_single("Anvi'o found the DSSP executable `%s`, and will use it."\
                                  % self.DSSP_executable, nl_before=1, nl_after=1)

        genes_missing_from_structure_db = [gene for gene in self.genes_of_interest if gene not in self.structure_db.genes_with_structure]
        if genes_missing_from_structure_db:
            show_a_few = genes_missing_from_structure_db if len(genes_missing_from_structure_db) <= 10 else genes_missing_from_structure_db[:10]
            raise ConfigError("{} gene(s) were specified by you but don't exist in the structure database. Here are some of their IDs: {}".
                format(len(genes_missing_from_structure_db), ', '.join([str(x) for x in show_a_few])))

        return genes_of_interest


    def load_structure_db(self):
        utils.is_structure_db(self.structure_db_path)
        self.structure_db = StructureDatabase(self.structure_db_path,
                                              ignore_hash=True,
                                              create_new=False)


class DSSPClass(object):
    def __init__(self, executable=None, check_executable=True):
        self.executable = executable
        self.check_executable = check_executable

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

        self.set_executable()
        utils.is_program_exists(self.executable)
        if self.check_executable:
            self.is_executable_a_working_DSSP_program()


    def run(self, pdb_filepath, id_for_protein=None):
        # Determine the model name by loading the structure file
        p = PDBParser()
        structure = p.get_structure(id_for_protein, pdb_filepath)
        model = structure[0] # pdb files can have multiple models. DSSP assumes the first.

        # run DSSP
        residue_annotation = DSSP(model, pdb_filepath, dssp = self.executable, acc_array = "Wilke")

        # convert to a digestible format
        return self.convert_DSSP_output_from_biopython_to_dataframe(residue_annotation, drop=['aa'])


    def set_executable(self):
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
        except Exception as e:
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


class ContactMap(object):
    def __init__(self, p=terminal.Progress(), r=terminal.Run()):
        self.distances_methods_dict = {
            'CA': self.calc_CA_dist,
        }


    def load_pdb_file(self, pdb_path, name_id = 'structure'):
        p = PDBParser()
        model = p.get_structure(name_id, pdb_path)[0] # [0] = first model
        structure = model['A'] # ['A'] = get first chain in model

        return structure


    def get_contact_map(self, pdb_path, distance_method='CA'):
        """Returns a contact map (pairwise distances in Angstroms)

        Parameters
        ==========
        pdb_path : str
            The PDB file (takes first chain)

        distance_method : str, 'CA'
            Which distance method? 'CA' is for alpha carbon distances. See
            self.distances_methods_dict for options

        Returns
        =======
        output: NxN 2-dim array
            Where N is the # of AAs

        Notes
        =====
        - See also `get_boolean_contact_map`
        """

        structure = self.load_pdb_file(pdb_path)

        contact_map = np.zeros((len(structure), len(structure)))
        for i, residue1 in enumerate(structure):
            for j, residue2 in enumerate(structure):
                contact_map[i, j] = self.distances_methods_dict[distance_method](residue1, residue2)

        return contact_map


    def get_boolean_contact_map(self, pdb_path, distance_method='CA', threshold=6):
        """Returns a boolean contact map (1 for touching, 0 for not)

        Parameters
        ==========
        pdb_path : str
            The PDB file (takes first chain)

        distance_method : str, 'CA'
            Which distance method? 'CA' is for alpha carbon distances. See
            self.distances_methods_dict for options

        threshold : int or float, 6
            Residues are considered touching if their distances are less than or equal to this
            amount. Units are in Angstroms

        Returns
        =======
        output: NxN 2-dim array
            Where N is the # of AAs

        Notes
        =====
        - See also `get_contact_map`
        """

        contact_map = self.get_contact_map(pdb_path, distance_method=distance_method)

        contact_map[contact_map <= threshold] = 1
        contact_map[contact_map >  threshold] = 0

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


    def calc_CA_dist(self, residue1, residue2):
        """Returns the C-alpha distance between two residues"""

        diff_vector = residue1["CA"].coord - residue2["CA"].coord
        return np.sqrt(np.sum(diff_vector**2))



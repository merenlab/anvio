# coding: utf-8
"""
Interface to MODELLER (https://salilab.org/modeller/).
"""

import os
import sys
import tempfile

import anvio
import shutil
import pandas as pd
import anvio.utils as utils
import anvio.fastalib as u
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths
import Bio.PDB as PDB

from anvio.errors import ConfigError, FilesNPathsError


__author__ = "A. Murat Eren"
__copyright__ = "Copyright 2016, The anvio Project"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "A. Murat Eren"
__email__ = "a.murat.eren@gmail.com"


run = terminal.Run()
progress = terminal.Progress()
pp = terminal.pretty_print

class MODELLER:

    def __init__(self, target_fasta_path, database="pdb_95", directory=None, run=run, progress=progress):
        self.run = run
        self.progress = progress

        # check that MODELLER exists
        utils.is_program_exists('mod9.19') # FIXME
        self.mod = "mod9.19"

        # All MODELLER scripts are housed in self.script_folder
        self.scripts_folder = "/Users/evan/Software/anvio/anvio/data/misc/MODELLER/scripts" # FIXME
        if utils.filesnpaths.is_dir_empty(self.scripts_folder):
            raise ConfigError("Anvi'o houses all its MODELLER scripts in {}, but your directory \
                               contains no scripts. You should change that.")

        # All MODELLER databases are housed in self.database_folder
        self.database = database
        self.database_folder = "/Users/evan/Software/anvio/anvio/data/misc/MODELLER/db" # FIXME

        # does target_fasta_path point to a fasta file?
        self.target_fasta_path = target_fasta_path
        utils.filesnpaths.is_file_fasta_formatted(self.target_fasta_path)
        run.info("target fasta", self.target_fasta_path)

        # make sure target_fasta is valid
        target_fasta = u.SequenceSource(target_fasta_path, lazy_init=False)
        if target_fasta.total_seq != 1:
            raise ConfigError("MODELLER::The input FASTA file must have exactly one sequence.\
                               You provided one with {}.".format(target_fasta.total_seq))

        # get gene_id while target_fasta is opened
        while next(target_fasta):
            self.gene_id = target_fasta.id
        target_fasta.close()

        # self.directory is defined because there are a lot of extraneous output and log
        # files output by MODELLER so we cd into directory, do our business, and then cd back
        # into the starting directory.
        self.directory = directory
        if not self.directory:
            self.directory = os.getcwd()
        self.start_dir = os.getcwd()
        os.chdir(self.directory)


    def process(self):
        """
        This is the workflow for a standard protein search and model.
        """
        try:
            self.run_fasta_to_pir()

            self.check_database()

            self.run_search()

            self.parse_search_results()

            self.download_structures()

            self.run_align_to_templates()

            self.run_get_model()

            self.close()

        except self.EndProcess:
            self.close()

    def run_get_model(self):
        """
        This is the magic of MODELLER. Based on the template alignment file, the structures of the
        templates, and satisfaction of physical constraints, the target protein structure is
        approximated.
        """
        script_name = "get-model.py"

        # check script exists, then copy the script into the working directory
        self.copy_script_to_directory(script_name)

        command = [self.mod,
                   script_name,
                   self.alignment_to_templates_pir,
                   self.gene_id,
                   self.best_template_ids_fname,
                   str(10),
                   str(3.0)]

        command = " ".join(command)

        self.run_command(command)
        run.info("get_model", "done")


    def parse_search_results(self, max_matches=5, min_proper_pident=30.0):
        """
        Parses search results and finds best matches to align against.
        """
        # put names to the columns
        column_names = (      "idx"      ,  "code_and_chain"  ,       "type"     ,  "iteration_num"  ,
                            "seq_len"    , "start_pos_target" , "end_pos_target" , "start_pos_dbseq" ,
                        "send_pos_dbseq" ,   "align_length"   ,   "seq_identity" ,      "evalue"      )

        # load the table as a dataframe
        search_df = pd.read_csv(self.search_results, sep="\s+", comment="#", names=column_names, index_col=False)

        # add some useful columns
        search_df["proper_pident"] = search_df["seq_identity"] * search_df["align_length"] / \
                                     search_df.iloc[0, search_df.columns.get_loc("seq_len")]
        search_df["code"] = search_df["code_and_chain"].str[:-1]
        search_df["chain"] = search_df["code_and_chain"].str[-1]

        # filter results by min_proper_pident FIXME I assume number of hits returned is not 0. Need
        # to learn about handling exceptions so process ends if length is zero. Order them and take
        # the first max_matches.
        search_df = search_df[search_df["proper_pident"] >= min_proper_pident]

        if not len(search_df):
            run.warning("Gene {} did not have a search result with percent identicalness above {}. No \
                         structure will be modelled".format(self.gene_id, min_proper_pident))
            raise self.EndProcess

        # of those filtered, get up to <max_matches> of those with the highest proper_ident scores.
        search_df = search_df.sort_values("proper_pident", ascending=False)
        search_df = search_df.iloc[:min([len(search_df), max_matches])]

        # Get their chain and 4-letter ids
        self.top_seq_matches = list(zip(search_df["code"], search_df["chain"]))


    def download_structures(self):
        """
        Downloads structure files for self.top_seq_seq_matches using Bioppython
        If the 4-letter code is `wxyz`, the downloaded file is `pdbwxyz.ent`.
        """
        # get complete list of PDB codes
        pdb_list = PDB.PDBList()

        # make a folder to store the template PDBs
        self.template_pdbs = os.path.join(self.directory, "TEMPLATE_PDBS")
        os.mkdir(self.template_pdbs)

        # function to get the filepath of downloaded pdb, given its 4-letter code
        func_fname = lambda x: os.path.join(self.template_pdbs, "pdb" + x + ".ent")

        for match in self.top_seq_matches:

            # download structure into self.directory; suppress Biopython output
            with HiddenPrints(): # FIXME the anvi'o SuppressAllOutput does not work
                pdb_list.retrieve_pdb_file(match[0], file_format="pdb", pdir=self.template_pdbs, overwrite=True)

            # if structure was not downloaded, that sucks. We remove it from self.top_seq_matches
            if utils.filesnpaths.is_file_exists(func_fname(match[0]), dont_raise=True):
                run.info(match[0], "downloaded")
            else:
                run.warning("{} could not be downloaded. Moving onto the next match.".format(match[0]))
                self.top_seq_matches.remove(match)

        if not len(self.top_seq_matches):
            run.warning("No structures of the candidate proteins were downloadable. Probably something is wrong. Stopping here.")
            raise self.EndProcess


    def run_align_to_templates(self):
        """
        After identifying best candidate proteins based on sequence data, this function aligns the
        candidate proteins based on structural and sequence information as well as to the target
        protein. This alignment file is the make input (besides structures) for the homology
        modelling aspect of MODELLER.
        """
        script_name = "align_to_templates.py"

        # check script exists, then copy the script into the working directory
        self.copy_script_to_directory(script_name)

        # First, write ids and chains to file read by align_to_templates.py MODELLER script
        self.best_template_ids_fname = "{}_best_template_ids.txt".format(self.gene_id)
        f = open(self.best_template_ids_fname, "w")
        for match in self.top_seq_matches:
            f.write("{}\t{}\n".format(match[0], match[1]))
        f.close()

        # name of the output. .pir is the standard format for MODELLER, .pap is human readable
        # protein_family computes a matrix comparing the different templates agianst one another
        self.alignment_to_templates_pir = "{}_align_to_templates.ali".format(self.gene_id)
        self.alignment_to_templates_pap = "{}_align_to_templates.pap".format(self.gene_id)
        self.template_family_matrix     = "{}_protein_family.mat".format(self.gene_id)

        command = [self.mod,
                   script_name,
                   self.target_pir,
                   self.best_template_ids_fname,
                   self.alignment_to_templates_pir,
                   self.alignment_to_templates_pap,
                   self.template_family_matrix]
        command = " ".join(command)

        self.run_command(command, check_output=[self.alignment_to_templates_pir, 
                                                self.alignment_to_templates_pap,
                                                self.template_family_matrix])
        run.info("run_align_to_templates", "done")


    def run_search(self):
        """
        Align the protein against the database based on sequence alone
        """
        script_name = "search.py"

        # check script exists, then copy the script into the working directory
        self.copy_script_to_directory(script_name)

        # name pir file by the gene_id (i.e. defline of the fasta)
        self.search_results = "{}_search_results.prf".format(self.gene_id)

        command = [self.mod,
                   script_name,
                   self.target_pir,
                   self.database_path,
                   self.search_results]
        command = " ".join(command)

        self.run_command(command, check_output=[self.search_results])
        run.info("search_results", self.search_results)


    def check_database(self):
        """
        Checks for the .bin version of database. If it only finds the .pir version, it binarizes it.
        Sets the db filepath.
        """
        extensionless, extension = os.path.splitext(self.database)
        if extension not in [".bin",".pir",""]:
            raise ConfigError("MODELLER :: The only possible database extensions are .bin and .pir")

        bin_db_path = os.path.join(self.database_folder, extensionless+".bin")
        pir_db_path = os.path.join(self.database_folder, extensionless+".pir")
        bin_exists = utils.filesnpaths.is_file_exists(bin_db_path, dont_raise=True)
        pir_exists = utils.filesnpaths.is_file_exists(pir_db_path, dont_raise=True)

        if pir_exists and bin_exists:
            self.database_path = bin_db_path
            return

        if not pir_exists and bin_exists:
            self.database_path = bin_db_path
            return

        if pir_exists and not bin_exists:
            self.run_binarize_database(pir_db_path, bin_db_path)
            self.database_path = bin_db_path
            return

        if not pir_exists and not bin_exists:
            raise ConfigError("The database {} does not exist in {}".format(self.database, self.database_folder))


    def run_binarize_database(self, pir_db_path, bin_db_path):
        """
        Databases can be read in .pir format, but can be more quickly read in binarized format. This
        does that.
        """
        script_name = "binarize-database.py"

        # check script exists, then copy the script into the working directory
        self.copy_script_to_directory(script_name)

        # name pir file by the gene_id (i.e. defline of the fasta)
        self.target_pir = "{}.pir".format(self.gene_id)

        command = [self.mod,
                   script_name,
                   pir_db_path,
                   bin_db_path]
        command = " ".join(command)

        self.run_command(command, check_output=[bin_db_path])
        run.info("binarized database", bin_db_path)


    def copy_script_to_directory(self, script_name):
        """
        All MODELLER scripts are housed in anvio/data/misc/MODELLER/scripts/. This function checks
        that script_name is in anvio/data/misc/MODELLER/scripts/ and then copies the script into
        self.directory. Why copy into self.directory? Whenever a script is ran by MODELLER, a log
        file is output in the directory of the script. By copying the script into self.directory,
        the log is written there instead of anvio/data/misc/MODELLER/scripts/. 
        """
        script_path = os.path.join(self.scripts_folder, script_name)
        try:
            utils.filesnpaths.is_file_exists(script_path)
        except:
            raise ConfigError("MODELLER :: The script {} is not in {}".format(script_name, self.scripts_folder))

        # If all is well, copy script to self.directory
        shutil.copy2(script_path, self.directory)


    def run_fasta_to_pir(self):
        """
        MODELLER uses their own .pir format for search and alignment instead of .fasta. This script
        does the conversion. An example is found at https://salilab.org/modeller/tutorial/basic.html
        """
        script_name = "fasta_to_pir.py"

        # check script exists, then copy the script into the working directory
        self.copy_script_to_directory(script_name)

        # name pir file by the gene_id (i.e. defline of the fasta)
        self.target_pir = "{}.pir".format(self.gene_id)

        command = [self.mod,
                   script_name,
                   self.target_fasta_path,
                   self.target_pir]
        command = " ".join(command)

        self.run_command(command, check_output=[self.target_pir])
        run.info("target pir", self.target_pir)


    def run_command(self, command, check_output=None):
        # FIXME jazz this up with some self.progress updates
        os.system(command)

        if check_output:
            for output in check_output:
                utils.filesnpaths.is_file_exists(output)


    def close(self):
        """
        Change directories back to self.start_dir. Could also have options selectively bring back
        log files and PDB files to self.start_dir if self.directory != self.start_dir.
        """
        # make back to starting directory
        os.chdir(self.start_dir)


    class EndProcess(Exception): 
        pass


class HiddenPrints:
    """
    Has your pesky external library hardcoded print statements that you want to supress? Are you
    sick and tired of seeing messages you aren't in control of? Do your friends laugh at you behind
    your back? Then you need the HiddenPrints class! Billy Mays here, introducing a class that
    executes chunks of code while suppressing all print statements. Simply call your code like:

    with HiddenPrints():
        <your code here!>
    """
    def __enter__(self):
        self._original_stdout = sys.stdout
        sys.stdout = open(os.devnull, 'w')

    def __exit__(self, exc_type, exc_val, exc_tb):
        sys.stdout = self._original_stdout















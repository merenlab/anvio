# coding: utf-8
"""
Interface to MODELLER (https://salilab.org/modeller/).
"""

import os
import sys
import anvio
import shutil
import subprocess

import pandas as pd
import anvio.utils as utils
import anvio.fastalib as u
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths

from anvio.errors import ConfigError, ModellerError


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
    """
    This class is a driver to run MODELLER scripts. MODELLER scripts are written
    in python 2.3 which is the language MODELLER uses to, you can make as many
    MODELLER scripts as you want, and they are all stored in
    anvio/data/misc/MODELLER/scripts. each script should have its own function
    in this class called run_<script_name>(self, <params>) which initializes the
    parameters required to run the script.
    """

    def __init__(self, target_fasta_path, database_name="pdb_95", directory=None, run=run, progress=progress):
        """ 
        PARAMS:
        =======

        target_fasta_path : str
            The fasta file for the protein you are trying to predict the
            structure of.  This is the only parameter needed for intialization.
            There should be exactly 1 sequence in the fasta, and it should be an
            amino acid sequence.

        database_name : str (default = "pdb_95")
            This is the name of the database used to finding homologous proteins
            from.  If you want to use your own database (it should be in .pir or
            .bin format), just put it in anvio/data/misc/MODELLER/db. There is
            no need to specify the extension when providing this parameter;
            anvi'o automatically uses a binary version of your database, and if
            only .pir format is available, anvi'o creates a binary version
            automatically. The default value is "pdb_95", which can be
            downloaded (and should updated periodically) from
            https://salilab.org/modeller/supplemental.html.  It is made by
            taking all structures from the PDB, clustering them into clusters of
            95% sequence identity, and then taking a single representative from
            each cluster.

        directory : str (default = current working directory)
            MODELLER outputs a lot of files into the current working directory.
            To manage this, you can supply a directory name where all MODELLER
            scripts from this instance will run in. Whenever a script method is
            called, anvi'o cds into directory, runs the MODELLER script, and
            then cds back into the original directory. By default everything is
            output into the current working directory, so it is highly
            recommended you specify a directory to be created (it will be
            created if it doesn't already exist)
        """

        self.run = run
        self.progress = progress

        self.alignment_pap_path = None
        self.alignment_pir_path = None
        self.get_template_path = None
        self.search_results_path = None
        self.target_pir_path = None
        self.target_fasta_path = None
        self.template_family_matrix_path = None
        self.template_info_path = None
        self.template_pdbs = None
        self.model_info = None

        self.logs = {}
        self.scripts = {}

        # check that MODELLER exists
        utils.is_program_exists('mod9.19') # FIXME
        self.mod = "mod9.19"

        # All MODELLER scripts are housed in self.script_folder
        self.scripts_folder = os.path.abspath("/Users/evan/Software/anvio/anvio/data/misc/MODELLER/scripts") # FIXME
        if utils.filesnpaths.is_dir_empty(self.scripts_folder):
            raise ConfigError("Anvi'o houses all its MODELLER scripts in {}, but your directory \
                               contains no scripts. You should change that.")

        # All MODELLER databases are housed in self.database_dir
        self.database_name = database_name
        self.database_dir = os.path.abspath("/Users/evan/Software/anvio/anvio/data/misc/MODELLER/db") # FIXME


        # does target_fasta_path point to a fasta file?
        self.target_fasta_path = target_fasta_path
        utils.filesnpaths.is_file_fasta_formatted(self.target_fasta_path)

        # make sure target_fasta is valid
        target_fasta = u.SequenceSource(target_fasta_path, lazy_init=False)
        if target_fasta.total_seq != 1:
            raise ConfigError("MODELLER::The input FASTA file must have exactly one sequence.\
                               You provided one with {}.".format(target_fasta.total_seq))

        # get gene_id while target_fasta is opened
        while next(target_fasta):
            self.gene_id = target_fasta.id
        target_fasta.close()

        # the directory files will be dumped into
        self.directory = directory
        if not self.directory:
            self.directory = os.getcwd()
        self.directory = filesnpaths.check_output_directory(self.directory, ok_if_exists=True)
        filesnpaths.gen_output_directory(self.directory)

        # copy fasta into the working directory
        try:
            shutil.copy2(self.target_fasta_path, self.directory)
            self.target_fasta_path = os.path.join(self.directory, self.target_fasta_path)
        except shutil.SameFileError:
            pass

        # store the original directory so we can cd back and forth between
        # self.directory and self.start_dir
        self.start_dir = os.getcwd()


    def pick_best_model(self, criteria):
        """
        The user is not interested in all of this output. They just want a pdb file for a given
        self.gene_id. renames as <gene_id>.pdb.

        criteria must be one of the ["molpdf", "GA341_score", "DOPE_score", "average"]
        """

        if criteria == "average":
            best_basename = "average.pdb"

        if criteria in ["molpself.model_info", "DOPE_score"]:
            best_basename = self.model_info.loc[self.model_info[criteria].idxmin(axis=0), "name"]

        if criteria == "GA341_score":
            best_basename = self.model_info.loc[self.model_info[criteria].idxmax(axis=0), "name"]

        new_best_file_path = os.path.join(self.directory, "{}.pdb".format(self.gene_id))
        os.rename(os.path.join(self.directory, best_basename), new_best_file_path)
        return new_best_file_path


    def abort(self):
        """
        For whatever reason, this gene was not modelled. Make sure we are in the directory we
        started in and remove any folders and files that just aren't needed anymore.
        """
        os.chdir(self.start_dir)

        stuff_to_remove = [self.template_pdbs,
                           self.scripts.get("align_to_templates.py"),
                           self.scripts.get("binarize_database.py"),
                           self.scripts.get("search.py"),
                           self.scripts.get("fasta_to_pir.py"),
                           self.scripts.get("get_model.py")]

        for thing in stuff_to_remove:
            try:
                if os.path.isfile(thing):
                    os.remove(thing)
                if os.path.isdir(thing):
                    shutil.rmtree(thing)
            except TypeError:
                continue


    def rewrite_model_info(self):
        """
        If changes rename changes or additional columns have been added to self.info_model, those
        changes can be reflected in model_info.txt if this function is ran.
        """
        self.model_info.to_csv(self.model_info_path, sep="\t", index=False)


    def tidyup(self): 
        """
        get_model.py has been ran. Some of the files in here are unnecessary, some of the names are
        disgusting. rename from "2.B99990001.pdb" to "model_1.pdb" if normal model. Rename from
        "cluster.ini" to "average_raw.pdb" and "cluster.opt" to "average.pdb"
        """
        if not "get_model.py" in self.scripts.keys():
            raise ConfigError("You are out of line calling tidyup without running get_model.py")

        # remove all copies of all scrips that were ran
        for script_name, file_path in self.scripts.items():
            os.remove(file_path)

        # what are the pdb file names? load model_info.txt to find out
        self.model_info = pd.read_csv(self.model_info_path, sep="\t", index_col=False)

        for model in self.model_info.index:
            basename = self.model_info.loc[model, "name"]

            # determine the new basename
            if basename == "cluster.opt":
                new_basename = "average.pdb"
            elif basename == "cluster.ini":
                new_basename = "average_raw.pdb"
            else:
                model_num = os.path.splitext(basename)[0][-1]
                new_basename = "model_{}.pdb".format(model_num)

            # rename the files (an reflect changes in self.model_info)
            file_path = os.path.join(self.directory, basename)
            new_file_path = os.path.join(self.directory, new_basename)
            os.rename(file_path, new_file_path)
            self.model_info.loc[model, "name"] = new_basename


    def run_get_model(self, num_models, deviation, very_fast):
        """
        This is the magic of MODELLER. Based on the template alignment file, the structures of the
        templates, and satisfaction of physical constraints, the target protein structure is
        modelled.
        """
        script_name = "get_model.py"

        # check script exists, then copy the script into the working directory
        self.copy_script_to_directory(script_name)

        # model info
        self.model_info_path = os.path.join(self.directory, "{}_model_info.txt".format(self.gene_id))

        self.run.info("Number of models", num_models)
        self.run.info("deviation", str(deviation) + " angstroms")
        self.run.info("fast optimization", str(very_fast))

        if not deviation and num_models > 1:
            raise ConfigError("run_get_model::deviation must be > 0 if num_models > 1.")

        command = [self.mod,
                   script_name,
                   self.alignment_pir_path,
                   self.gene_id,
                   self.template_info_path,
                   str(num_models),
                   str(deviation),
                   str(int(very_fast)),
                   self.model_info_path]

        self.run_command(command, 
                         script_name = script_name,
                         progress_name = "CALCULATING 3D MODEL")

        self.run.info("model info", os.path.basename(self.model_info_path))


    def run_align_to_templates(self, templates_info):
        """
        After identifying best candidate proteins based on sequence data, this function aligns the
        candidate proteins based on structural and sequence information as well as to the target
        protein. This alignment file is the make input (besides structures) for the homology
        modelling aspect of MODELLER.

        templates_info is a list of 2-tuples; the zeroth element is the 4-letter protein code
        and the first element is the chain number. an example is [('4sda', 'A'), ('iq8p', 'E')]
        """
        script_name = "align_to_templates.py"

        # check script exists, then copy the script into the working directory
        self.copy_script_to_directory(script_name)

        # First, write ids and chains to file read by align_to_templates.py MODELLER script
        self.template_info_path = os.path.join(self.directory, "{}_best_template_ids.txt".format(self.gene_id))
        f = open(self.template_info_path, "w")
        for match in templates_info:
            f.write("{}\t{}\n".format(match[0], match[1]))
        f.close()

        # name of the output. .pir is the standard format for MODELLER, .pap is human readable
        # protein_family computes a matrix comparing the different templates agianst one another
        self.alignment_pir_path = os.path.join(self.directory, "{}_alignment.ali".format(self.gene_id))
        self.alignment_pap_path = os.path.join(self.directory, "{}_alignment.pap".format(self.gene_id))
        self.template_family_matrix_path = os.path.join(self.directory, "{}_protein_family.mat".format(self.gene_id))

        command = [self.mod,
                   script_name,
                   self.target_pir_path,
                   self.gene_id,
                   self.template_info_path,
                   self.alignment_pir_path,
                   self.alignment_pap_path,
                   self.template_family_matrix_path]

        self.run_command(command, 
                         script_name = script_name, 
                         progress_name = "CREATING MSA OF HOMOLOGS",
                         check_output = [self.alignment_pir_path, 
                                         self.alignment_pap_path,
                                         self.template_family_matrix_path])

        self.run.info("Similarity matrix of templates", os.path.basename(self.template_family_matrix_path))
        self.run.info("Target alignment to templates", ", ".join([os.path.basename(self.alignment_pir_path),
                                                                  os.path.basename(self.alignment_pap_path)]))


    def run_search(self):
        """
        Align the protein against the database based on sequence alone
        """
        script_name = "search.py"

        # check script exists, then copy the script into the working directory
        self.copy_script_to_directory(script_name)

        # name pir file by the gene_id (i.e. defline of the fasta)
        self.search_results_path = os.path.join(self.directory, "{}_search_results.prf".format(self.gene_id))

        command = [self.mod,
                   script_name,
                   self.target_pir_path,
                   self.database_path,
                   self.search_results_path]

        self.run_command(command, 
                         script_name = script_name, 
                         progress_name = "DB SEARCH FOR HOMOLOGS",
                         check_output = [self.search_results_path])

        run.info("Search results for similar AA sequences", os.path.basename(self.search_results_path))


    def check_database(self):
        """
        Checks for the .bin version of database. If it only finds the .pir version, it binarizes it.
        Sets the db filepath.
        """
        extensionless, extension = os.path.splitext(self.database_name)
        if extension not in [".bin",".pir",""]:
            raise ConfigError("MODELLER :: The only possible database extensions are .bin and .pir")

        bin_db_path = os.path.join(self.database_dir, extensionless+".bin")
        pir_db_path = os.path.join(self.database_dir, extensionless+".pir")
        bin_exists = utils.filesnpaths.is_file_exists(bin_db_path, dont_raise=True)
        pir_exists = utils.filesnpaths.is_file_exists(pir_db_path, dont_raise=True)

        if pir_exists and bin_exists:
            self.database_path = bin_db_path
            return

        if not pir_exists and bin_exists:
            self.database_path = bin_db_path
            return

        if pir_exists and not bin_exists:
            self.run.warning("Your database is not in binary format. That means accessing its contents is slower \
                              than it could be. Anvi'o is going to make a binary format. Just FYI")
            self.run_binarize_database(pir_db_path, bin_db_path)
            self.database_path = bin_db_path
            return

        if not pir_exists and not bin_exists:
            raise ConfigError("Anvi'o looked in {} for a database with the name {} and with an extension \
                               of either .bin or .pir, but didn't find anything matching that criteria. Here \
                               are the files anvi'o did find: {}. You should take a look at 00_README, which gives \
                               some easy instructions for obtaining a database :) Here is its full path: {}".\
                               format(self.database_dir, 
                                      self.database_name, 
                                      ", ".join(os.listdir(self.database_dir)),
                                      os.path.abspath(os.path.join(self.database_dir, "00_README"))))



    def run_binarize_database(self, pir_db_path, bin_db_path):
        """
            Databases can be read in .pir format, but can be more quickly read in binarized format. This
            does that.
        """
        script_name = "binarize_database.py"

        # check script exists, then copy the script into the working directory
        self.copy_script_to_directory(script_name)

        # name pir file by the gene_id (i.e. defline of the fasta)
        self.target_pir_path = "{}.pir".format(self.gene_id)

        command = [self.mod,
                   script_name,
                   pir_db_path,
                   bin_db_path]

        self.run_command(command, 
                         script_name = script_name, 
                         progress_name = "BINARIZING DATABASE",
                         check_output=[bin_db_path])

        self.run.info("new database", bin_db_path)


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

        # add script to scripts dictionary
        self.scripts[script_name] = os.path.join(self.directory, script_name)

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
        self.target_pir_path = os.path.join(self.directory, "{}.pir".format(self.gene_id))

        command = [self.mod,
                   script_name,
                   self.target_fasta_path,
                   self.target_pir_path]

        self.run_command(command, 
                         script_name = script_name, 
                         progress_name = "CONVERT SEQUENCE TO MODELLER FORMAT", 
                         check_output = [self.target_pir_path])

        self.run.info("Target alignment file", os.path.basename(self.target_pir_path))


    def run_command(self, command, script_name, progress_name, check_output=None):
        """
        Runs a script. Must provide script_name (e.g. "script.py") and progress_name (e.g.
        "BINARIZING DATABASE"). Optionally can provide list of output files whose existences are
        checked to make sure command was successfully ran. We ALWAYS cd into the MODELLER's
        directory before running a script, and we ALWAYS cd back into the original directory after
        running a script.
        """
        # first things first, we CD into MODELLER's directory
        os.chdir(self.directory)

        # run the command
        self.progress.new(progress_name)
        self.progress.update("Executing MODELLER script: {}".format(script_name))

        # try and execute the command
        process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        output, error = process.communicate()

        # if MODELLER script gave a traceback, it is caught here and everything is stopped
        if process.returncode: 
            #format the error
            error = str(error).replace("\\n", "\n").replace("\\'","\'")[2:-1].strip()
            self.progress.end()
            self.run.warning(error, header="\nTraceback message for {}".format(script_name), lc='red', raw=True)
            raise ModellerError("The MODELLER script {} did not execute properly. Hopefully it is clear \
                                 from the above error message what went wrong. If you think you may have \
                                 accidentally messed with this script, you can go to \
                                 https://github.com/merenlab/anvio/tree/master/anvio/data/misc/MODELLER/scripts \
                                 to replace it. Otherwise, you have identified a bug and should report \
                                 it to the anvi'o developers. Congratulations *streamers*".format(script_name))

        # If we made it this far, the MODELLER script ran to completion. Now check outputs exist
        if check_output:
            for output in check_output:
                utils.filesnpaths.is_file_exists(output)

        # MODELLER outputs a log that we rename right here, right now
        old_log_name = os.path.splitext(script_name)[0] + ".log"
        new_log_name = "{}_{}".format(self.gene_id, old_log_name)
        os.rename(old_log_name, new_log_name)

        # add to logs
        self.logs[script_name] = new_log_name

        self.progress.end()
        self.run.info("log of {}".format(script_name), new_log_name)

        # last things last, we CD back into the starting directory
        os.chdir(self.start_dir)


    class EndModeller(Exception): 
        pass


class HiddenPrints():
    """
    Has your pesky external library hardcoded print statements that you want to supress? Are you
    sick and tired of seeing messages you aren't in control of? Do your friends laugh at you behind
    your back? Does Meren's terminal.SuppressAllOutput lead to a traceback? Then you need the
    HiddenPrints class! Billy Mays here, introducing a class that executes chunks of code while
    suppressing all print statements. Simply call your code like:

    with HiddenPrints():
        <your code here!>
    """
    def __enter__(self):
        self._original_stdout = sys.stdout
        sys.stdout = open(os.devnull, 'w')

    def __exit__(self, exc_type, exc_val, exc_tb):
        sys.stdout = self._original_stdout




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


pp = terminal.pretty_print

J = lambda x, y: os.path.join(x, y)

class MODELLER:
    """
    This class is a driver to run MODELLER scripts. MODELLER scripts are written
    in python 2.3 which is the language MODELLER uses to, you can make as many
    MODELLER scripts as you want, and they are all stored in
    anvio/data/misc/MODELLER/scripts. each script should have its own function
    in this class called run_<script_name>(self, <params>) which initializes the
    parameters required to run the script.
    """

    def __init__(self, args, run=terminal.Run(), progress=terminal.Progress()):

        self.args = args
        self.run = run
        self.progress = progress

        A = lambda x, t: t(args.__dict__[x]) if x in self.args.__dict__ else None
        null = lambda x: x
        self.best = A('best', str)
        self.deviation = A('deviation', float)
        self.directory = A('directory', str)
        self.very_fast = A('very_fast', bool)
        self.executable = A('executable', str) or "mod9.19"
        self.num_models = A('num_models', int)
        self.target_fasta_path = A('target_fasta_path', str)
        self.database_name = A('database_name', str) or "pdb_95"

        self.alignment_pap_path = None
        self.alignment_pir_path = None
        self.get_template_path = None
        self.search_results_path = None
        self.target_pir_path = None
        self.template_family_matrix_path = None
        self.template_info_path = None
        self.template_pdbs = None
        self.model_info = None

        self.logs = {}
        self.scripts = {}

        self.sanity_check()

        # All MODELLER databases are housed in self.database_dir
        self.database_dir = J(os.path.dirname(anvio.__file__), 'data/misc/MODELLER/db')

        # copy fasta into the working directory
        try:
            shutil.copy2(self.target_fasta_path, self.directory)
            self.target_fasta_path = J(self.directory, self.target_fasta_path)
        except shutil.SameFileError:
            pass

        # store the original directory so we can cd back and forth between
        # self.directory and self.start_dir
        self.start_dir = os.getcwd()


    def sanity_check(self):
        A = lambda x, t: t(args.__dict__[x]) if x in self.args.__dict__ else None
        null = lambda x: x

        # the directory files will be dumped into (can exist but must be empty)
        if filesnpaths.is_file_exists(self.directory, dont_raise=True):
            filesnpaths.is_output_dir_writable(self.directory)
            if not filesnpaths.is_dir_empty(self.directory):
                raise ModellerError("You cannot give MODELLER a non-empty directory to work in.")
        else:
            filesnpaths.gen_output_directory(self.directory)

        # check that MODELLER exists
        utils.is_program_exists(self.executable)

        if A('executable', null):
            self.run.warning("As per your request, anvi'o will use %s to run MODELER." % self.executable)
        else:
            self.run.info_single("Anvi'o found the default executable for MODELLER, %s, and will\
                                  use it." % self.executable, nl_before=1)

        # All MODELLER scripts are housed in self.script_folder
        self.scripts_folder = J(os.path.dirname(anvio.__file__), 'data/misc/MODELLER/scripts')
        if utils.filesnpaths.is_dir_empty(self.scripts_folder):
            raise ConfigError("Anvi'o houses all its MODELLER scripts in {}, but your directory \
                               contains no scripts. Why you do dat?")

        # does target_fasta_path point to a fasta file?
        utils.filesnpaths.is_file_fasta_formatted(self.target_fasta_path)

        # make sure target_fasta is valid
        target_fasta = u.SequenceSource(self.target_fasta_path, lazy_init=False)
        if target_fasta.total_seq != 1:
            raise ConfigError("MODELLER::The input FASTA file must have exactly one sequence.\
                               You provided one with {}.".format(target_fasta.total_seq))

        # (not sanity check but we get self.gene_id since target_fasta is opened)
        while next(target_fasta):
            self.gene_id = target_fasta.id
        target_fasta.close()

        # parameter consistencies
        if self.deviation < 0.5 or self.deviation > 20:
            self.run.warning("You realize that deviation is given in angstroms, right? You chose {}".format(self.deviation))

        if self.very_fast and self.num_models > 1:
            self.run.warning("Since you chose --very-fast, there will be little difference, if at all, between models. You \
                              can potentially save a lot of time by setting --num-models to 1.")


    def pick_best_model(self):
        """
        The user is not interested in all of this output. They just want a pdb file for a given
        self.gene_id. renames as gene_<gene_id>.pdb.

        self.best must be one of the ["molpdf", "GA341_score", "DOPE_score", "average"]
        """

        # initialize new model_info column
        self.model_info["picked_as_best"] = False

        if self.best == "average":
            best_basename = "average.pdb"
            self.model_info.loc[self.model_info["name"] == "average.pdb", "picked_as_best"] = True

        # For these scores, lower is better
        if self.best in ["molpself.model_info", "DOPE_score"]:
            best_basename = self.model_info.loc[self.model_info[self.best].idxmin(axis=0), "name"]
            self.model_info.loc[self.model_info[self.best].idxmin(axis=0), "picked_as_best"] = True

        # For these scores, higher is better
        if self.best == "GA341_score":
            best_basename = self.model_info.loc[self.model_info[self.best].idxmax(axis=0), "name"]
            self.model_info.loc[self.model_info[self.best].idxmax(axis=0), "picked_as_best"] = True

        new_best_file_path = J(self.directory, "gene_{}.pdb".format(self.gene_id))
        os.rename(J(self.directory, best_basename), new_best_file_path)

        # add new column to self.model_info
        self.model_info.loc

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
        self.model_info.to_csv(self.model_info_path, sep="\t", index=False, na_rep="N/A")


    def tidyup(self): 
        """
        get_model.py has been ran. Some of the files in here are unnecessary, some of the names are
        disgusting. rename from "2.B99990001.pdb" to "gene_2_Model001.pdb" if normal model. Rename from
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

            # The default names are bad. This is where they are defined
            if basename == "cluster.opt":
                new_basename = "gene_{}_ModelAvg.pdb".format(self.gene_id)
            else:
                model_num = os.path.splitext(basename)[0][-3:]
                new_basename = "gene_{}_Model{}.pdb".format(self.gene_id, model_num)

            # rename the files (an reflect changes in self.model_info)
            file_path = J(self.directory, basename)
            new_file_path = J(self.directory, new_basename)
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
        self.model_info_path = J(self.directory, "gene_{}_ModelInfo.txt".format(self.gene_id))

        self.run.info("Number of models", num_models)
        self.run.info("deviation", str(deviation) + " angstroms")
        self.run.info("fast optimization", str(very_fast))

        if not deviation and num_models > 1:
            raise ConfigError("run_get_model::deviation must be > 0 if num_models > 1.")

        command = [self.executable,
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
        self.template_info_path = J(self.directory, "gene_{}_BestTemplateIDs.txt".format(self.gene_id))
        f = open(self.template_info_path, "w")
        for match in templates_info:
            f.write("{}\t{}\n".format(match[0], match[1]))
        f.close()

        # name of the output. .pir is the standard format for MODELLER, .pap is human readable
        # protein_family computes a matrix comparing the different templates agianst one another
        self.alignment_pir_path = J(self.directory, "gene_{}_Alignment.ali".format(self.gene_id))
        self.alignment_pap_path = J(self.directory, "gene_{}_Alignment.pap".format(self.gene_id))
        self.template_family_matrix_path = J(self.directory, "gene_{}_ProteinFamily.mat".format(self.gene_id))

        command = [self.executable,
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
        self.search_results_path = J(self.directory, "gene_{}_SearchResults.prf".format(self.gene_id))

        command = [self.executable,
                   script_name,
                   self.target_pir_path,
                   self.database_path,
                   self.search_results_path]

        self.run_command(command, 
                         script_name = script_name, 
                         progress_name = "DB SEARCH FOR HOMOLOGS",
                         check_output = [self.search_results_path])

        self.run.info("Search results for similar AA sequences", os.path.basename(self.search_results_path))


    def check_database(self):
        """
        Checks for the .bin version of database. If it only finds the .pir version, it binarizes it.
        Sets the db filepath.
        """
        extensionless, extension = os.path.splitext(self.database_name)
        if extension not in [".bin",".pir",""]:
            raise ConfigError("MODELLER :: The only possible database extensions are .bin and .pir")

        bin_db_path = J(self.database_dir, extensionless+".bin")
        pir_db_path = J(self.database_dir, extensionless+".pir")
        bin_exists = utils.filesnpaths.is_file_exists(bin_db_path, dont_raise=True)
        pir_exists = utils.filesnpaths.is_file_exists(pir_db_path, dont_raise=True)

        self.database_path = bin_db_path

        if bin_exists:
            return

        if not pir_exists and not bin_exists:
            self.run.warning("Anvi'o looked in {} for a database with the name {} and with an extension \
                              of either .bin or .pir, but didn't find anything matching that \
                              criteria. We'll try and download the best database we know of from \
                              https://salilab.org/modeller/downloads/pdb_95.pir.gz and use that.".\
                              format(self.database_dir, 
                                     self.database_name))

            db_download_path = os.path.join(self.database_dir, "pdb_95.pir.gz")
            utils.download_file("https://salilab.org/modeller/downloads/pdb_95.pir.gz", db_download_path)
            utils.run_command(['gzip', '-d', db_download_path], log_file_path=filesnpaths.get_temp_file_path())

            pir_exists = utils.filesnpaths.is_file_exists(pir_db_path, dont_raise=True)

        if pir_exists and not bin_exists:
            self.run.warning("Your database is not in binary format. That means accessing its contents is slower \
                              than it could be. Anvi'o is going to make a binary format. Just FYI")
            self.run_binarize_database(pir_db_path, bin_db_path)
            return


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

        command = [self.executable,
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
        script_path = J(self.scripts_folder, script_name)
        try:
            utils.filesnpaths.is_file_exists(script_path)
        except:
            raise ConfigError("MODELLER :: The script {} is not in {}".format(script_name, self.scripts_folder))

        # add script to scripts dictionary
        self.scripts[script_name] = J(self.directory, script_name)

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
        self.target_pir_path = J(self.directory, "{}.pir".format(self.gene_id))

        command = [self.executable,
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
            self.progress.end()
            error = error.decode('utf-8').strip()

            is_licence_key_error = True if error.find('Invalid license key') > -1 else False            

            if is_licence_key_error:
                license_target_file = error.split('\n')[-1]
                raise ModellerError("MODELLER could not find your licence key. Please go to https://salilab.org/modeller/ and \
                                    get a new license. After you receive an e-mail with your key, please open '%s' \
                                    and replace 'XXXXX' with your key, save the file and try again. " % license_target_file)
            else:
                error = "\n".join(error.split('\n')[2:-1])
                print(terminal.c(error, color='red'))
                raise ModellerError("The MODELLER script {} did not execute properly. Hopefully it is clear \
                                     from the above error message".format(script_name))

        # If we made it this far, the MODELLER script ran to completion. Now check outputs exist
        if check_output:
            for output in check_output:
                utils.filesnpaths.is_file_exists(output)

        # MODELLER outputs a log that we rename right here, right now
        old_log_name = os.path.splitext(script_name)[0] + ".log"
        new_log_name = "gene_{}_{}".format(self.gene_id, old_log_name)
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




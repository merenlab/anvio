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

from anvio.errors import ConfigError, ModellerError, ModellerScriptError


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

    def __init__(self, args, run=terminal.Run(), progress=terminal.Progress(), progress_title=None):

        self.args = args
        self.run = run
        self.progress = progress

        up_to_date_modeller_exec = "mod9.20" # default exec to use

        A = lambda x, t: t(args.__dict__[x]) if x in self.args.__dict__ else None
        null = lambda x: x
        self.scoring_method = A('scoring_method', str)
        self.deviation = A('deviation', float)
        self.directory = A('directory', str)
        self.very_fast = A('very_fast', bool)
        self.executable = A('modeller_executable', null) or up_to_date_modeller_exec
        self.num_models = A('num_models', int)
        self.target_fasta_path = A('target_fasta_path', str)
        self.modeller_database = A('modeller_database', str) or "pdb_95"
        self.max_number_templates = A('max_number_templates', null)
        self.percent_identical_cutoff = A('percent_identical_cutoff', null)
        self.deviation = A('deviation', null)

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

        # as reward, whoever called this class will receive self.out when they run self.process()
        self.out = {
            "templates"                : {"pdb_id": [],"chain_id": [],"ppi": []},
            "models"                   : {"molpdf": [],"GA341_score": [],"DOPE_score": [],"picked_as_best": []},
            "corresponding_gene_call"  : self.corresponding_gene_call,
            "structure_exists"         : False,
            "best_model_path"          : None,
            "best_score"               : None,
            "scoring_method"           : self.scoring_method,
            "percent_identical_cutoff" : self.percent_identical_cutoff,
            "very_fast"                : self.very_fast,
            "deviation"                : self.deviation,
            }

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

        self.progress_title = progress_title
        if not self.progress_title:
            self.progress_title = "Running MODELLER for gene id {}".format(self.corresponding_gene_call)


    def process(self):

        self.run.warning("Working directory: {}".format(self.directory),
                         header='Modelling structure for gene ID {}'.format(self.corresponding_gene_call),
                         lc="green")

        self.progress.new(self.progress_title)

        try:
            self.run_fasta_to_pir()

            self.check_database()

            self.run_search()

            self.parse_search_results()

            self.download_structures()

            self.run_align_to_templates(self.top_seq_matches)

            self.run_get_model(self.num_models, self.deviation, self.very_fast)

            self.tidyup()

            self.pick_best_model()

            self.out["structure_exists"] = True

        except self.EndModeller as e:
            print(e)

        except ModellerScriptError as e:
            print(e)

        finally:
            self.abort()

        self.progress.end()
        return self.out


    def download_structures(self):
        """
        Downloads structure files for self.top_seq_seq_matches using Biopython
        If the 4-letter code is `wxyz`, the downloaded file is `pdbwxyz.ent`.
        """
        self.progress.update("Downloading homologs from PDB")

        # define directory path name to store the template PDBs (it can already exist)
        self.template_pdbs = os.path.join(self.directory, "{}_TEMPLATE_PDBS".format(self.corresponding_gene_call))

        downloaded = utils.download_protein_structures([code[0] for code in self.top_seq_matches], self.template_pdbs)

        # redefine self.top_seq_matches in case not all were downloaded
        self.top_seq_matches = [(code, chain_code) for code, chain_code in self.top_seq_matches if code in downloaded]

        if not len(self.top_seq_matches):
            self.progress.end()
            self.run.warning("No structures of the homologous proteins (templates) were downloadable. Probably something \
                              is wrong. Maybe you are not connected to the internet. Stopping here.")
            raise self.EndModeller

        self.run.info("Structures downloaded for", ", ".join([code[0] for code in self.top_seq_matches]), progress=self.progress)


    def parse_search_results(self):
        """
        Parses search results and filters for best homologs to use as structure templates.

        parameters used :: self.percent_identical_cutoff, self.max_number_templates
        """
        if not self.percent_identical_cutoff or not self.max_number_templates:
            raise ConfigError("parse_search_results::You initiated this class without providing values for percent_identical_cutoff \
                               and max_number_templates, which is required for this function.")

        self.progress.update("Parsing and filtering homologs")

        # put names to the columns
        column_names = (      "idx"      ,  "code_and_chain"  ,       "type"     ,  "iteration_num"  ,
                            "seq_len"    , "start_pos_target" , "end_pos_target" , "start_pos_dbseq" ,
                        "send_pos_dbseq" ,   "align_length"   ,   "seq_identity" ,      "evalue"      )

        # load the table as a dataframe
        search_df = pd.read_csv(self.search_results_path, sep="\s+", comment="#", names=column_names, index_col=False)

        # matches found (-1 because target is included)
        matches_found = len(search_df) - 1

        if not matches_found:
            self.progress.end()
            self.run.warning("No proteins with homologous sequence were found for {}. No structure will be modelled".format(self.corresponding_gene_call))
            raise self.EndModeller

        # add some useful columns
        search_df["proper_pident"] = search_df["seq_identity"] * search_df["align_length"] / \
                                     search_df.iloc[0, search_df.columns.get_loc("seq_len")]
        search_df["code"] = search_df["code_and_chain"].str[:-1]
        search_df["chain"] = search_df["code_and_chain"].str[-1]

        # filter results by self.percent_identical_cutoff.
        max_pident_found = search_df["proper_pident"].max()
        id_of_max_pident = tuple(search_df.loc[search_df["proper_pident"].idxmax(), ["code", "chain"]].values)
        search_df = search_df[search_df["proper_pident"] >= self.percent_identical_cutoff]

        # Order them and take the first self.modeller.max_number_templates.
        matches_after_filter = len(search_df)
        if not matches_after_filter:
            self.progress.end()
            self.run.warning("Gene {} did not have a search result with proper percent identicalness above or equal \
                              to {}%. The best match was chain {} of https://www.rcsb.org/structure/{}, which had a\
                              proper percent identicalness of {:.2f}%. No structure will be modelled.".\
                              format(self.corresponding_gene_call,
                                     self.percent_identical_cutoff,
                                     id_of_max_pident[1],
                                     id_of_max_pident[0],
                                     max_pident_found))
            raise self.EndModeller

        # of those filtered, get up to self.modeller.max_number_templates of those with the highest proper_ident scores.
        search_df = search_df.sort_values("proper_pident", ascending=False)
        search_df = search_df.iloc[:min([len(search_df), self.max_number_templates])]

        # Get their chain and 4-letter ids
        self.top_seq_matches = list(zip(search_df["code"], search_df["chain"]))

        self.run.info("Max number of templates allowed", self.max_number_templates, progress=self.progress)
        self.run.info("Number of candidate templates", matches_found, progress=self.progress)
        self.run.info("After >{}% identical filter".format(self.percent_identical_cutoff), matches_after_filter, progress=self.progress)
        self.run.info("Number accepted as templates", len(self.top_seq_matches), progress=self.progress)

        # update user on which templates are used, and write the templates to self.out
        for i in range(len(self.top_seq_matches)):
            pdb_id, chain_id = self.top_seq_matches[i]
            ppi = search_df["proper_pident"].iloc[i]

            self.out["templates"]["pdb_id"].append(pdb_id)
            self.out["templates"]["chain_id"].append(chain_id)
            self.out["templates"]["ppi"].append(ppi)

            self.run.info("Template {}".format(i+1),
                          "Protein ID: {}, Chain {} ({:.1f}% identical)".format(pdb_id, chain_id, ppi),
                          progress=self.progress)


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

        # All MODELLER scripts are housed in self.script_folder
        self.scripts_folder = J(os.path.dirname(anvio.__file__), 'data/misc/MODELLER/scripts')
        if utils.filesnpaths.is_dir_empty(self.scripts_folder):
            raise ConfigError("Anvi'o houses all its MODELLER scripts in {}, but your directory \
                               contains no scripts. Why you do dat?")

        # check that MODELLER exists
        if self.args.__dict__['modeller_executable'] if 'modeller_executable' in self.args.__dict__ else None:
            self.run.info_single("As per your request, anvi'o will use `%s` to run MODELLER." % self.executable, nl_before=1)
            utils.is_program_exists(self.executable)
        else:
            try:
                utils.is_program_exists(self.executable)
            except ConfigError as e:
                raise ConfigError("Anvi'o needs a MODELLER program to be installed on your system. You didn't specify one\
                                   (which can be done with `--modeller-executable`), so anvi'o tried the most recent version\
                                   it knows about: '%s'. If you are certain you have it on your system (for instance you can run it\
                                   by typing '%s' in your terminal window), you may want to send a detailed bug report. If you\
                                   don't have it on your system, check out these installation instructions on our website:\
                                   http://merenlab.org/2016/06/18/installing-third-party-software/#modeller" % (self.executable, self.executable))

            self.run.info_single("Anvi'o found the default executable for MODELLER, `%s`, and will\
                                  use it." % self.executable, nl_before=1)
        self.is_executable_a_MODELLER_program()

        # does target_fasta_path point to a fasta file?
        utils.filesnpaths.is_file_fasta_formatted(self.target_fasta_path)

        # make sure target_fasta is valid
        target_fasta = u.SequenceSource(self.target_fasta_path, lazy_init=False)
        if target_fasta.total_seq != 1:
            raise ConfigError("MODELLER::The input FASTA file must have exactly one sequence.\
                               You provided one with {}.".format(target_fasta.total_seq))

        # (not sanity check but we get self.corresponding_gene_call since target_fasta is opened)
        while next(target_fasta):
            self.corresponding_gene_call = target_fasta.id
        target_fasta.close()

        # parameter consistencies
        if self.deviation < 0.5 or self.deviation > 20:
            self.run.warning("You realize that deviation is given in angstroms, right? You chose {}".format(self.deviation))

        if self.very_fast and self.num_models > 1:
            self.run.warning("Since you chose --very-fast, there will be little difference, if at all, between models. You \
                              can potentially save a lot of time by setting --num-models to 1.")

        if self.percent_identical_cutoff <= 20:
            self.run.warning("Two completely unrelated sequences of same length can expect to have around 10% proper \
                              percent identicalness... Having this parameter below 20% is probably a bad idea.")


    def pick_best_model(self):
        """
        The user is not interested in all of this output. They just want a pdb file for a given
        self.corresponding_gene_call. renames as gene_<corresponding_gene_call>.pdb.

        self.scoring_method must be one of the ["molpdf", "GA341_score", "DOPE_score"]
        """

        # initialize new model_info column
        self.model_info["picked_as_best"] = False

        # For these scores, lower is better
        if self.scoring_method in ["molpself.model_info", "DOPE_score"]:
            best_basename = self.model_info.loc[self.model_info[self.scoring_method].idxmin(axis=0), "name"]
            self.model_info.loc[self.model_info[self.scoring_method].idxmin(axis=0), "picked_as_best"] = True

        # For these scores, higher is better
        if self.scoring_method == "GA341_score":
            best_basename = self.model_info.loc[self.model_info[self.scoring_method].idxmax(axis=0), "name"]
            self.model_info.loc[self.model_info[self.scoring_method].idxmax(axis=0), "picked_as_best"] = True

        new_best_file_path = J(self.directory, "gene_{}.pdb".format(self.corresponding_gene_call))
        os.rename(J(self.directory, best_basename), new_best_file_path)

        # append model information to self.out
        for model_index in self.model_info.index:
            self.out["models"]["molpdf"].append(self.model_info.loc[model_index, "molpdf"])
            self.out["models"]["GA341_score"].append(self.model_info.loc[model_index, "GA341_score"])
            self.out["models"]["DOPE_score"].append(self.model_info.loc[model_index, "DOPE_score"])
            self.out["models"]["picked_as_best"].append(self.model_info.loc[model_index, "picked_as_best"])

        # append pdb path to self.out
        self.out["best_model_path"] = new_best_file_path

        # append the best score to self.out
        self.out["best_score"] = self.model_info.loc[self.model_info["picked_as_best"] == True, self.scoring_method]


    def abort(self):
        """
        For whatever reason, this gene was not modelled. Make sure we are in the directory we
        started in.
        """
        os.chdir(self.start_dir)


    def tidyup(self): 
        """
        get_model.py has been ran. Some of the files in here are unnecessary, some of the names are
        disgusting. rename from "2.B99990001.pdb" to "gene_2_Model001.pdb" if normal model. Rename
        from "cluster.opt" to "gene_2_ModelAvg.pdb"
        """
        if not "get_model.py" in self.scripts.keys():
            raise ConfigError("You are out of line calling tidyup without running get_model.py")

        # remove all copies of all scrips that were ran
        for script_name, file_path in self.scripts.items():
            os.remove(file_path)

        for model in self.model_info.index:
            basename = self.model_info.loc[model, "name"]

            # The default names are bad. This is where they are defined
            if basename == "cluster.opt":
                new_basename = "gene_{}_ModelAvg.pdb".format(self.corresponding_gene_call)
            else:
                model_num = os.path.splitext(basename)[0][-3:]
                new_basename = "gene_{}_Model{}.pdb".format(self.corresponding_gene_call, model_num)

            # rename the files (an reflect changes in self.model_info)
            file_path = J(self.directory, basename)
            new_file_path = J(self.directory, new_basename)
            os.rename(file_path, new_file_path)
            self.model_info.loc[model, "name"] = new_basename


    def run_get_model(self, num_models, deviation, very_fast):
        """
        This is the magic of MODELLER. Based on the template alignment file, the structures of the
        templates, and satisfaction of physical constraints, the target protein structure is
        modelled without further user input.
        """
        script_name = "get_model.py"

        # check script exists, then copy the script into the working directory
        self.copy_script_to_directory(script_name)

        # model info
        self.model_info_path = J(self.directory, "gene_{}_ModelInfo.txt".format(self.corresponding_gene_call))

        self.run.info("Number of models", num_models, progress=self.progress)
        self.run.info("Deviation", str(deviation) + " angstroms", progress=self.progress)
        self.run.info("Fast optimization", str(very_fast), progress=self.progress)

        if not deviation and num_models > 1:
            raise ConfigError("run_get_modeli :: deviation must be > 0 if num_models > 1.")

        command = [self.executable,
                   script_name,
                   self.alignment_pir_path,
                   self.corresponding_gene_call,
                   self.template_info_path,
                   str(num_models),
                   str(deviation),
                   str(int(very_fast)),
                   self.model_info_path]

        self.run_command(command, 
                         script_name = script_name,
                         progress_update = "Calculating 3D model(s)",
                         check_output = [self.model_info_path])

        # load the model results information as a dataframe
        self.model_info = pd.read_csv(self.model_info_path, sep="\t", index_col=False)

        self.run.info("Model info", os.path.basename(self.model_info_path), progress=self.progress)


    def run_align_to_templates(self, templates_info):
        """
        After identifying best candidate proteins based on sequence data, this function aligns the
        protein. This alignment file is the main input (besides structures) for the homology
        protein.

        templates_info is a list of 2-tuples; the zeroth element is the 4-letter protein code
        and the first element is the chain number. an example is [('4sda', 'A'), ('iq8p', 'E')]
        """
        script_name = "align_to_templates.py"

        # check script exists, then copy the script into the working directory
        self.copy_script_to_directory(script_name)

        # First, write ids and chains to file read by align_to_templates.py MODELLER script
        self.template_info_path = J(self.directory, "gene_{}_BestTemplateIDs.txt".format(self.corresponding_gene_call))
        f = open(self.template_info_path, "w")
        for match in templates_info:
            f.write("{}\t{}\n".format(match[0], match[1]))
        f.close()

        # name of the output. .pir is the standard format for MODELLER, .pap is human readable
        # protein_family computes a matrix comparing the different templates agianst one another
        self.alignment_pir_path = J(self.directory, "gene_{}_Alignment.ali".format(self.corresponding_gene_call))
        self.alignment_pap_path = J(self.directory, "gene_{}_Alignment.pap".format(self.corresponding_gene_call))
        self.template_family_matrix_path = J(self.directory, "gene_{}_ProteinFamily.mat".format(self.corresponding_gene_call))

        command = [self.executable,
                   script_name,
                   self.target_pir_path,
                   self.corresponding_gene_call,
                   self.template_info_path,
                   self.alignment_pir_path,
                   self.alignment_pap_path,
                   self.template_family_matrix_path]

        self.run_command(command, 
                         script_name = script_name, 
                         progress_update = "Aligning sequence to template structures",
                         check_output = [self.alignment_pir_path, 
                                         self.alignment_pap_path,
                                         self.template_family_matrix_path])

        self.run.info("Similarity matrix of templates", os.path.basename(self.template_family_matrix_path), progress=self.progress)
        self.run.info("Target alignment to templates", ", ".join([os.path.basename(self.alignment_pir_path),
                                                                  os.path.basename(self.alignment_pap_path)]),
                                                                  progress=self.progress)


    def run_search(self):
        """
        Align the protein against the database based on sequence alone
        """
        script_name = "search.py"

        # check script exists, then copy the script into the working directory
        self.copy_script_to_directory(script_name)

        # name pir file by the corresponding_gene_call (i.e. defline of the fasta)
        self.search_results_path = J(self.directory, "gene_{}_SearchResults.prf".format(self.corresponding_gene_call))

        command = [self.executable,
                   script_name,
                   self.target_pir_path,
                   self.database_path,
                   self.search_results_path]

        self.run_command(command, 
                         script_name = script_name, 
                         progress_update = "Searching DB for sequence homologs",
                         check_output = [self.search_results_path])

        self.run.info("Search results for similar AA sequences", os.path.basename(self.search_results_path), progress=self.progress)


    def check_database(self):
        """
        Checks for the .bin version of database. If it only finds the .pir version, it binarizes it.
        Sets the db filepath.
        """
        extensionless, extension = os.path.splitext(self.modeller_database)
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
            self.progress.clear()
            self.run.warning("Anvi'o looked in {} for a database with the name {} and with an extension \
                              of either .bin or .pir, but didn't find anything matching that \
                              criteria. We'll try and download the best database we know of from \
                              https://salilab.org/modeller/downloads/pdb_95.pir.gz and use that. \
                              You can checkout https://salilab.org/modeller/ for more info about the pdb_95 \
                              database".format(self.database_dir, self.modeller_database))

            db_download_path = os.path.join(self.database_dir, "pdb_95.pir.gz")
            utils.download_file("https://salilab.org/modeller/downloads/pdb_95.pir.gz", db_download_path)
            utils.run_command(['gzip', '-d', db_download_path], log_file_path=filesnpaths.get_temp_file_path())

            pir_exists = utils.filesnpaths.is_file_exists(pir_db_path, dont_raise=True)

        if pir_exists and not bin_exists:
            self.progress.clear()
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

        # name pir file by the corresponding_gene_call (i.e. defline of the fasta)
        self.target_pir_path = "{}.pir".format(self.corresponding_gene_call)

        command = [self.executable,
                   script_name,
                   pir_db_path,
                   bin_db_path]

        self.run_command(command, 
                         script_name = script_name, 
                         progress_update = "Binarizing database",
                         check_output=[bin_db_path])

        self.run.info("New database", bin_db_path, progress=self.progress)


    def copy_script_to_directory(self, script_name, add_to_scripts_dict=True, directory=None):
        """
        All MODELLER scripts are housed in anvio/data/misc/MODELLER/scripts/. This function checks
        that script_name is in anvio/data/misc/MODELLER/scripts/ and then copies the script into
        self.directory. Why copy into self.directory? Whenever a script is ran by MODELLER, a log
        file is output in the directory of the script. By copying the script into self.directory,
        the log is written there instead of anvio/data/misc/MODELLER/scripts/. 
        """
        if not directory:
            directory = self.directory

        script_path = J(self.scripts_folder, script_name)
        try:
            utils.filesnpaths.is_file_exists(script_path)
        except:
            raise ConfigError("MODELLER :: The script {} is not in {}".format(script_name, self.scripts_folder))

        # add script to scripts dictionary
        if add_to_scripts_dict:
            self.scripts[script_name] = J(directory, script_name)

        # If all is well, copy script to directory
        shutil.copy2(script_path, directory)


    def is_executable_a_MODELLER_program(self):
        # temp_dir created because log file outputs to wherever fasta_to_pir.py is
        temp_dir = filesnpaths.get_temp_directory_path()
        self.copy_script_to_directory('fasta_to_pir.py', add_to_scripts_dict=False, directory=temp_dir)
        test_script = J(temp_dir, 'fasta_to_pir.py')
        test_input = os.path.abspath(J(os.path.dirname(anvio.__file__), '../tests/sandbox/mock_data_for_structure/proteins.fa'))
        test_output = J(temp_dir, 'test_out')

        command = [self.executable,
                   test_script,
                   test_input,
                   test_output]

        # try and execute the command
        process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        output, error = process.communicate()

        if process.returncode:
            # modeller has failed
            error = error.decode('utf-8').strip()

            is_licence_key_error = True if error.find('Invalid license key') > -1 else False
            if is_licence_key_error:
                # its a valid modeller program with no license key
                license_target_file = error.split('\n')[-1]
                raise ConfigError("You're making progress and anvi'o is proud of you! You just need to validate your MODELLER\
                                   with a license key (it's free). Please go to https://salilab.org/modeller/registration.html\
                                   to register for a new license. After you receive an e-mail with your key, please open '%s'\
                                   and replace the characters XXXXX with your own key. Save the file and try again. " % license_target_file)

            else:
                error = "\n" + "\n".join(error.split('\n'))
                print(terminal.c(error, color='red'))
                raise ConfigError("The executable you requested is called `%s`, but anvi'o doesn't agree with you that\
                                   it is a working MODELLER program. That was determined by running the command `%s`, which raised the\
                                   error seen above. If you want to specify a specific MODELLER program, you can specify it with\
                                   `--modeller-executable`."
                                       % (self.executable, " ".join(command)))

        # no error was raised. now check if output file exists
        try:
            filesnpaths.is_file_exists(test_output)
        except FilesNPathsError:
            raise ConfigError("The executable you requested is called `%s`, but anvi'o doesn't agree with you that\
                               it is a working MODELLER program. That was determined by running the command `%s`, which did not\
                               output the file expected. If you want to specify a specific MODELLER program, you can specify it with\
                               `--modeller-executable`." % (self.executable, " ".join(command)))


    def run_fasta_to_pir(self):
        """
        MODELLER uses their own .pir format for search and alignment instead of .fasta. This script
        does the conversion. An example is found at https://salilab.org/modeller/tutorial/basic.html
        """
        script_name = "fasta_to_pir.py"

        # check script exists, then copy the script into the working directory
        self.copy_script_to_directory(script_name)

        # name pir file by the corresponding_gene_call (i.e. defline of the fasta)
        self.target_pir_path = J(self.directory, "{}.pir".format(self.corresponding_gene_call))

        command = [self.executable,
                   script_name,
                   self.target_fasta_path,
                   self.target_pir_path]

        self.run_command(command, 
                         script_name = script_name, 
                         progress_update = "Convert FASTA to MODELLER format", 
                         check_output = [self.target_pir_path])

        self.run.info("Target alignment file", os.path.basename(self.target_pir_path), progress=self.progress)


    def run_command(self, command, script_name, progress_update, check_output=None):
        """
        Runs a script. Must provide script_name (e.g. "script.py") and progress_update (e.g.
        "BINARIZING DATABASE"). Optionally can provide list of output files whose existences are
        checked to make sure command was successfully ran. We ALWAYS cd into the MODELLER's
        directory before running a script, and we ALWAYS cd back into the original directory after
        running a script.
        """
        # first things first, we CD into MODELLER's directory
        os.chdir(self.directory)

        # run the command
        self.progress.update(progress_update)

        # try and execute the command
        process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        output, error = process.communicate()

        # if MODELLER script gave a traceback, it is caught here and everything is stopped
        if process.returncode: 
            self.progress.end()
            error = error.decode('utf-8').strip()

            error = "\n" + "\n".join(error.split('\n'))
            print(terminal.c(error, color='red'))
            self.out["structure_exists"] = False
            raise ModellerScriptError("The MODELLER script {} did not execute properly. Hopefully it is clear \
                                       from the above error message. No structure is going to be modelled."\
                                       .format(script_name))

        # If we made it this far, the MODELLER script ran to completion. Now check outputs exist
        if check_output:
            for output in check_output:
                utils.filesnpaths.is_file_exists(output)

        # MODELLER outputs a log that we rename right here, right now
        old_log_name = os.path.splitext(script_name)[0] + ".log"
        new_log_name = "gene_{}_{}".format(self.corresponding_gene_call, old_log_name)
        os.rename(old_log_name, new_log_name)

        # add to logs
        self.logs[script_name] = new_log_name

        self.run.info("Log of {}".format(script_name), new_log_name, progress=self.progress)

        # last things last, we CD back into the starting directory
        os.chdir(self.start_dir)


    class EndModeller(Exception):
        pass



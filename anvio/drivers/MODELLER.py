# coding: utf-8
"""
Interface to MODELLER (https://salilab.org/modeller/).
"""

import os
import anvio
import shutil
import argparse
import subprocess

import pandas as pd
import anvio.fastalib as u
import anvio.terminal as terminal
import anvio.constants as constants
import anvio.filesnpaths as filesnpaths

from anvio.drivers import diamond
from anvio.errors import ConfigError, ModellerError, ModellerScriptError, FilesNPathsError
from anvio.utils.commandline import run_command
from anvio.utils.network import download_file, download_protein_structure
from anvio.utils.system import is_program_exists

__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "Evan Kiefl"
__email__ = "kiefl.evan@gmail.com"


up_to_date_modeller_exec = "mod10.1" # default exec to use

J = lambda x, y: os.path.join(x, y)


class MODELLER:
    """Driver class for MODELLER

    This class is a driver to run MODELLER scripts. MODELLER scripts are written
    in python 2.3 which is the language MODELLER used when this driver was written.

    Parameters
    ==========
    args : argparse.Namespace object
        Check __init__ for allowable attributes

    target_fasta_path: str
        Path to amino acid sequence fasta file with 1 sequence, the gene to be modelled. The defline
        should be an integer (This class will assume this integer is the genes gene caller id)

    directory: str, None
        Path to directory that MODELLER will be run in. If None, temp dir will be created

    lazy_init : bool, False
        If True, check_MODELLER will not be called

    skip_warnings : bool, False
        If True, all warnings will be suppressed

    Notes
    =====
    - You can add MODELLER scripts by storing them in anvio/data/misc/MODELLER/scripts. Each script
      should have its own function in this class. For example, align_to_templates.py is a script
      anvi'o has found in that directory and has a corresponding function in this class called
      self.run_align_to_templates. Please see that method if you want to add your own script.
    """

    def __init__(self, args, target_fasta_path, directory=None, run=terminal.Run(),
                 lazy_init=False, skip_warnings=False, check_db_only=False):

        self.args = args
        self.run = run
        if skip_warnings and not anvio.DEBUG:
            self.run.verbose = False
        self.lazy_init = lazy_init

        self.check_db_only = check_db_only
        self.target_fasta_path = target_fasta_path
        self.directory = directory if directory else filesnpaths.get_temp_directory_path()

        A = lambda x, t: t(args.__dict__[x]) if x in self.args.__dict__ else None
        null = lambda x: x
        self.scoring_method = A('scoring_method', str) or 'DOPE_score'
        self.very_fast = A('very_fast', bool) or False
        self.executable = A('modeller_executable', null) or up_to_date_modeller_exec
        self.num_models = A('num_models', int) or 5
        self.modeller_database = A('modeller_database', str) or 'pdb_95'
        self.max_number_templates = A('max_number_templates', null) or 5
        self.percent_cutoff = A('percent_cutoff', null) or 30
        self.alignment_fraction_cutoff = A('alignment_fraction_cutoff', null) or 0.80
        self.deviation = A('deviation', null) or 4
        self.pdb_db_path = A('pdb_db', null)
        self.offline_mode = A('offline_mode', null)

        # All MODELLER scripts are housed in self.script_folder
        self.scripts_folder = constants.default_modeller_scripts_dir

        self.alignment_pap_path = None
        self.alignment_pir_path = None
        self.get_template_path = None
        self.target_pir_path = None
        self.template_family_matrix_path = None
        self.template_info_path = None
        self.template_pdb_dir = None
        self.model_info = None
        self.pdb_db = None
        self.use_pdb_db = False

        self.logs = {}
        self.scripts = {}

        # All MODELLER databases are housed in self.database_dir
        self.database_dir = constants.default_modeller_database_dir

        # store the original directory so we can cd back and forth between
        # self.directory and self.start_dir
        self.start_dir = os.getcwd()

        if self.check_db_only:
            self.check_database()
            return

        self.sanity_check()
        self.corresponding_gene_call = self.get_corresponding_gene_call_from_target_fasta_path()

        # as reward, whoever called this class will receive self.out when they run self.process()
        self.out = {
            "templates"                 : {"pdb_id": [], "chain_id": [], "proper_percent_similarity": [], "percent_similarity": [], "align_fraction":[]},
            "models"                    : {"molpdf": [], "GA341_score": [], "DOPE_score": [], "picked_as_best": []},
            "corresponding_gene_call"   : self.corresponding_gene_call,
            "structure_exists"          : False,
            "best_model_path"           : None,
            "best_score"                : None,
            "scoring_method"            : self.scoring_method,
            "percent_cutoff"            : self.percent_cutoff,
            "alignment_fraction_cutoff" : self.alignment_fraction_cutoff,
            "very_fast"                 : self.very_fast,
            "deviation"                 : self.deviation,
            "directory"                 : self.directory,
        }

        # copy fasta into the working directory
        try:
            shutil.copy2(self.target_fasta_path, self.directory)
            self.target_fasta_path = J(self.directory, self.target_fasta_path)
        except shutil.SameFileError:
            pass


    def get_corresponding_gene_call_from_target_fasta_path(self):
        """corresponding_gene_call is assumed to be the defline of self.args.target_fasta_path"""

        target_fasta = u.SequenceSource(self.target_fasta_path, lazy_init=False)
        while next(target_fasta):
            corresponding_gene_call = target_fasta.id
        target_fasta.close()

        return corresponding_gene_call


    def load_pdb_db(self):
        """Try loading a PDB database with path equal to self.pdb_db_path

        Modifies self.pdb_db and self.use_pdb_db
        """

        if not self.pdb_db_path:
            self.pdb_db_path = constants.default_pdb_database_path

        ok_if_absent = True if self.pdb_db_path == constants.default_pdb_database_path else False

        if filesnpaths.is_file_exists(self.pdb_db_path, dont_raise=ok_if_absent):
            # The user has a database there! Try and load it
            self.pdb_db = anvio.structureops.PDBDatabase(argparse.Namespace(pdb_database_path=self.pdb_db_path))
            self.pdb_db.check_or_create_db()
            self.pdb_db.get_stored_structure_ids()
            self.use_pdb_db = True
        else:
            self.use_pdb_db = False


    def process(self):
        timer = terminal.Timer()

        try:
            self.load_pdb_db()
            timer.make_checkpoint('PDB DB loaded')

            self.run_fasta_to_pir()
            timer.make_checkpoint('Converted gene FASTA to PIR')

            self.check_database()
            timer.make_checkpoint('Checked databases')

            self.run_search_and_parse_results()
            timer.make_checkpoint('Ran DIAMOND search and parsed hits')

            self.get_structures()
            timer.make_checkpoint('Obtained template structures')

            self.run_align_to_templates(self.list_of_template_code_and_chain_ids)
            timer.make_checkpoint('Sequence aligned to templates')

            self.run_get_model(self.num_models, self.deviation, self.very_fast)
            timer.make_checkpoint('Ran structure predictions')

            self.tidyup()
            self.pick_best_model()
            self.run_add_chain_identifiers_to_best_model()
            timer.make_checkpoint('Picked best model and tidied up')

            self.out["structure_exists"] = True

        except self.EndModeller:
            pass

        except ModellerScriptError as e:
            print(e)

        finally:
            timer.gen_report(title='ID %s Time Report' % str(self.corresponding_gene_call), run=self.run)
            self.abort()

        return self.out


    def get_structures(self):
        """Populate self.template_pdb_dir with template structure PDBs"""

        self.template_pdb_dir = os.path.join(self.directory, "%s_TEMPLATE_PDBS" % str(self.corresponding_gene_call))
        filesnpaths.gen_output_directory(self.template_pdb_dir) # does nothing if already exists

        pdb_paths = {}
        for code, chain in self.list_of_template_code_and_chain_ids:
            five_letter_id = code + chain
            requested_path = J(self.template_pdb_dir, '%s.pdb' % code)

            if self.use_pdb_db and five_letter_id in self.pdb_db.stored_structure_ids:
                # This chain exists in the external database. Export it and get the path
                try:
                    path = self.pdb_db.export_pdb(five_letter_id, requested_path)
                    source = 'Offline DB'
                except ConfigError:
                    # The ID is in the DB, but the PDB content is None
                    path = None
                    source = 'Nowhere'

            elif not self.offline_mode:
                # This chain doesn't exist in an external database, and internet access is assumed.
                # We try and download the protein from the RCSB PDB server. If downloading fails,
                # path is None
                path = download_protein_structure(code, chain=chain, output_path=requested_path, raise_if_fail=False)
                source = 'RCSB PDB Server'

            else:
                # Internet access is not assumed, and the chain wasn't in the external database
                path = None
                source = 'Nowhere'

            self.run.info('%s obtained from' % five_letter_id, source)

            if path:
                pdb_paths[five_letter_id] = path

        # remove templates whose structures are not available
        self.list_of_template_code_and_chain_ids = [
            (code, chain_code)
            for code, chain_code in self.list_of_template_code_and_chain_ids
            if code + chain_code in pdb_paths
        ]

        if not len(self.list_of_template_code_and_chain_ids):
            self.run.warning("No structures of the homologous proteins (templates) were available. Probably something "
                             "is wrong. Stopping here.")
            raise self.EndModeller

        self.run.info("Structures obtained for", ", ".join([code[0]+code[1] for code in self.list_of_template_code_and_chain_ids]))


    def sanity_check(self, skip_warnings=False):
        # the directory files will be dumped into (can exist but must be empty)
        if filesnpaths.is_file_exists(self.directory, dont_raise=True):
            filesnpaths.is_output_dir_writable(self.directory)
            if not filesnpaths.is_dir_empty(self.directory):
                raise ModellerError("You cannot give MODELLER a non-empty directory to work in.")
        else:
            filesnpaths.gen_output_directory(self.directory)

        if not self.lazy_init:
            self.executable = check_MODELLER(self.executable)

        # does target_fasta_path point to a fasta file?
        utils.filesnpaths.is_file_fasta_formatted(self.target_fasta_path)

        # make sure target_fasta is valid
        target_fasta = u.SequenceSource(self.target_fasta_path, lazy_init=False)
        if target_fasta.total_seq != 1:
            raise ConfigError("MODELLER :: The input FASTA file must have exactly one sequence. "
                              "You provided one with {}.".format(target_fasta.total_seq))
        try:
            while next(target_fasta):
                int(target_fasta.id)
        except:
            raise ConfigError("MODELLER :: The defline of this fasta file must be an integer")
        target_fasta.close()

        # parameter consistencies
        if self.deviation < 0.5 or self.deviation > 20:
            self.run.warning("You realize that deviation is given in angstroms, right? You chose {}".format(self.deviation))

        if self.very_fast and self.num_models > 1:
            self.num_models = 1
            self.run.warning("Since you chose --very-fast, there will be little difference, if at all, between models. Anvi'o "
                             "authoritatively sets --num-models to 1 to save you time.")


    def pick_best_model(self):
        """Pick best model based on self.scoring_method and rename to gene_<corresponding_gene_call>.pdb"""

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
        """Gene was not modelled. Return to the starting directory"""

        os.chdir(self.start_dir)


    def tidyup(self):
        """Tidyup operations after running get_model.py

        Some of the files in here are unnecessary, some of the names are disgusting. rename from
        "2.B99990001.pdb" to "gene_2_Model001.pdb" if normal model. Rename from "cluster.opt" to
        "gene_2_ModelAvg.pdb"
        """

        if not "get_model.py" in self.scripts.keys():
            raise ConfigError("You are out of line calling tidyup without running get_model.py")

        # remove all copies of all scripts that were ran
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


    def run_add_chain_identifiers_to_best_model(self):
        """Add chain identifier to best model to appease some third-party services"""

        script_name = "add_chain_identifiers_to_best_model.py"

        # check script exists, then copy the script into the working directory
        self.copy_script_to_directory(script_name)

        dir_name, base_name = os.path.split(self.out['best_model_path'])

        command = [self.executable,
                   script_name,
                   dir_name,
                   base_name]

        self.run_command(command, script_name=script_name)


    def run_get_model(self, num_models, deviation, very_fast):
        """Run get model

        This is the magic of MODELLER. Based on the template alignment file, the structures of the
        templates, and satisfaction of physical constraints, the target protein structure is
        modelled without further user input.
        """

        script_name = "get_model.py"

        # check script exists, then copy the script into the working directory
        self.copy_script_to_directory(script_name)

        # model info
        self.model_info_path = J(self.directory, "gene_{}_ModelInfo.txt".format(self.corresponding_gene_call))

        self.run.info("Number of models", num_models)
        self.run.info("Deviation", str(deviation) + " angstroms")
        self.run.info("Fast optimization", str(very_fast))

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
                         check_output = [self.model_info_path])

        # load the model results information as a dataframe
        self.model_info = pd.read_csv(self.model_info_path, sep="\t", index_col=False)

        self.run.info("Model info", os.path.basename(self.model_info_path))


    def run_align_to_templates(self, templates_info):
        """Align the sequence to the best candidate protein sequences

        After identifying best candidate proteins based on sequence data, this function aligns the
        protein. This alignment file is the main input (besides structures) for the homology
        protein.

        Parameters
        ==========
        templates_info : list of 2-tuples
            The zeroth element is the 4-letter protein code and the first element is the chain
            number. E.g. [('4sda', 'A'), ('iq8p', 'E')]
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
                         check_output = [self.alignment_pir_path,
                                         self.alignment_pap_path,
                                         self.template_family_matrix_path])

        self.run.info("Similarity matrix of templates", os.path.basename(self.template_family_matrix_path))
        self.run.info("Target alignment to templates", ", ".join([os.path.basename(self.alignment_pir_path),
                                                                  os.path.basename(self.alignment_pap_path)]))


    def run_search_and_parse_results(self):
        """Align the protein against the database based on only sequence"""

        # Change to MODELLER working directory
        os.chdir(self.directory)

        columns = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gaps', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore']
        driver = diamond.Diamond(
            query_fasta=self.target_fasta_path,
            target_fasta=J(self.database_dir, self.modeller_database + '.dmnd'),
            outfmt=' '.join(['6'] + columns),
            run=terminal.Run(verbose=False),
            progress=terminal.Progress(verbose=False),
        )
        driver.blastp()

        # Change back to user directory
        os.chdir(self.start_dir)

        search_df = driver.view_as_dataframe(J(self.directory, driver.tabular_output_path))

        matches_found = search_df.shape[0]

        if not matches_found:
            self.run.warning("No proteins with homologous sequence were found for {}. No structure will be modelled".format(self.corresponding_gene_call))
            raise self.EndModeller

        # We need the gene length for pident
        target_fasta = u.SequenceSource(self.target_fasta_path, lazy_init=False)
        while next(target_fasta):
            gene_length = len(target_fasta.seq)

        # add some useful columns
        search_df["code"] = search_df["sseqid"].str[:-1]
        search_df["chain"] = search_df["sseqid"].str[-1]
        search_df["align_fraction"] = (search_df["length"] - search_df["gaps"]) / gene_length
        search_df["proper_pident"] = search_df["pident"] * search_df["align_fraction"]

        # Find best match for align fraction and pident
        code_chain_id_of_best = tuple(search_df.iloc[search_df['proper_pident'].argmax()][['code', 'chain']].values)
        best_hit = search_df.loc[
            (search_df['code'] == code_chain_id_of_best[0]) & \
            (search_df['chain'] == code_chain_id_of_best[1]), ['pident', 'align_fraction']
        ].iloc[0]

        # filter results by self.percent_cutoff and self.alignment_fraction_cutoff
        search_df = search_df[search_df["pident"] >= self.percent_cutoff]
        search_df = search_df[search_df["align_fraction"] >= self.alignment_fraction_cutoff]

        # Rank by the alignment fraction times the percent id
        search_df = search_df.sort_values("proper_pident", ascending=False)

        # If more than 1 template in 1 PDB id, just choose 1
        search_df = search_df.drop_duplicates('code', keep='first')

        matches_after_filter = len(search_df)
        if not matches_after_filter:
            self.run.warning("Gene {} did not have a search result with percent identicalness above or equal "
                             "to {}% and alignment fraction above {}%. The best match was chain {} of https://www.rcsb.org/structure/{}, which had a "
                             "percent identicalness of {:.2f}% and an alignment fraction of {:.3f}. No structure will be modelled.".\
                              format(self.corresponding_gene_call,
                                     self.percent_cutoff,
                                     self.alignment_fraction_cutoff,
                                     code_chain_id_of_best[1],
                                     code_chain_id_of_best[0],
                                     best_hit['pident'],
                                     best_hit['align_fraction']))
            raise self.EndModeller

        # Filter out templates with proper_pident more than 5% less than best match
        # http://merenlab.org/2018/09/04/getting-started-with-anvi-structure/#how-much-do-templates-matter
        search_df = search_df[search_df['proper_pident'] >= (search_df['proper_pident'].max() - 5)]

        # get up to self.modeller.max_number_templates of those with the highest proper_ident scores.
        search_df = search_df.iloc[:min([len(search_df), self.max_number_templates])]

        # Get their chain and 4-letter ids
        self.list_of_template_code_and_chain_ids = list(zip(search_df["code"], search_df["chain"]))

        self.run.info("Max number of templates allowed", self.max_number_templates)
        self.run.info("Number of candidate templates", matches_found)
        self.run.info("After >{}% identical filter".format(self.percent_cutoff), matches_after_filter)
        self.run.info("Number accepted as templates", len(self.list_of_template_code_and_chain_ids))

        # update user on which templates are used, and write the templates to self.out
        for i in range(len(self.list_of_template_code_and_chain_ids)):
            pdb_id, chain_id = self.list_of_template_code_and_chain_ids[i]
            proper_percent_similarity = search_df["proper_pident"].iloc[i]
            percent_similarity = search_df["pident"].iloc[i]
            align_fraction = search_df["align_fraction"].iloc[i]

            self.out["templates"]["pdb_id"].append(pdb_id)
            self.out["templates"]["chain_id"].append(chain_id)
            self.out["templates"]["proper_percent_similarity"].append(proper_percent_similarity)
            self.out["templates"]["percent_similarity"].append(percent_similarity)
            self.out["templates"]["align_fraction"].append(align_fraction)

            self.run.info("Template {}".format(i+1),
                          "Protein ID: {}, Chain {} ({:.1f}% identical, {:.2f} align fraction)".format(pdb_id, chain_id, percent_similarity, align_fraction))


    def check_database(self):
        """Setup the database files

        Downloads the .pir file if it is missing
        Binarizes .pir file if .bin is missing
        Creates the .dmnd file if it is missing
        """

        bin_db_path = J(self.database_dir, self.modeller_database + ".bin")
        pir_db_path = J(self.database_dir, self.modeller_database + ".pir")
        bin_exists = utils.filesnpaths.is_file_exists(bin_db_path, dont_raise=True)
        pir_exists = utils.filesnpaths.is_file_exists(pir_db_path, dont_raise=True)

        if bin_exists and pir_exists:
            # We good
            pass
        else:
            if not pir_exists:
                # Download .pir
                self.run.warning("Anvi'o looked in {} for a database with the name {} and with an extension "
                                 "of either .bin or .pir, but didn't find anything matching that "
                                 "criteria. Anvi'o will try and download the best database it knows of from "
                                 "https://salilab.org/modeller/downloads/pdb_95.pir.gz and use that. "
                                 "You can checkout https://salilab.org/modeller/ for more info about the pdb_95 "
                                 "database".format(self.database_dir, self.modeller_database))

                db_download_path = os.path.join(self.database_dir, "pdb_95.pir.gz")
                download_file("https://salilab.org/modeller/downloads/pdb_95.pir.gz", db_download_path)
                run_command(['gzip', '-d', db_download_path], log_file_path=filesnpaths.get_temp_file_path())

            # Binarize .pir (make .bin)
            self.run.warning("Your database is not in binary format. That means accessing its contents is slower "
                             "than it could be. Anvi'o is going to make a binary format. Just FYI")
            self.run_binarize_database(pir_db_path, bin_db_path)

        dmnd_db_path = J(self.database_dir, self.modeller_database + '.dmnd')

        if os.path.exists(dmnd_db_path):
            return

        self.run.warning("Your diamond database does not exist. It will be created.")

        script_name = "pir_to_fasta.py"

        self.copy_script_to_directory(script_name)

        input_pir_path = J(self.database_dir, self.modeller_database + '.pir')
        fasta_path = J(self.database_dir, self.modeller_database + '.fa')
        dmnd_path = J(self.database_dir, self.modeller_database)

        command = [self.executable,
                   script_name,
                   input_pir_path,
                   fasta_path]

        self.run_command(command,
                         script_name=script_name,
                         rename_log=False)

        temp = u.FastaOutput(filesnpaths.get_temp_file_path())
        fasta = u.SequenceSource(fasta_path)

        while next(fasta):
            temp.write_id(fasta.id)
            temp.write_seq(fasta.seq.replace('-', '').replace('.', 'X'))

        shutil.move(temp.output_file_path, fasta_path)
        fasta.close()
        temp.close()

        driver = diamond.Diamond(
            query_fasta=fasta_path,
            run=terminal.Run(verbose=False),
            progress=terminal.Progress(verbose=False),
        )
        driver.makedb(output_file_path=dmnd_path)

        os.remove(fasta_path)


    def run_binarize_database(self, pir_db_path, bin_db_path):
        """Binarizes a .pir file

        Databases can be read in .pir format, but can be more quickly read in binarized format. This
        does that.

        Parameters
        ==========
        pir_db_path : str
            Path to existing .pir file

        bin_db_path : str
            Path to the will-be-made .bin file
        """

        script_name = "binarize_database.py"

        # check script exists, then copy the script into the working directory
        self.copy_script_to_directory(script_name)

        command = [self.executable,
                   script_name,
                   pir_db_path,
                   bin_db_path]

        self.run_command(command,
                         script_name=script_name,
                         check_output=[bin_db_path],
                         rename_log=False)

        self.run.info("New database", bin_db_path)


    def copy_script_to_directory(self, script_name, add_to_scripts_dict=True, directory=None):
        """Copy source script to working directory

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


    def run_fasta_to_pir(self):
        """Convert a fasta file to a pir format.

        MODELLER uses their own .pir format for search and alignment instead of .fasta. This script
        does the conversion. An example pir formatted sequence shown here:

            >P1;TvLDH
            sequence:TvLDH:::::::0.00: 0.00
            MSEAAHVLITGAAGQIGYILSHWIASGELYGDRQVYLHLLDIPPAMNRLTALTMELEDCAFPHLAGFVATTDPKA
            AFKDIDCAFLVASMPLKPGQVRADLISSNSVIFKNTGEYLSKWAKPSVKVLVIGNPDNTNCEIAMLHAKNLKPEN
            FSSLSMLDQNRAYYEVASKLGVDVKDVHDIIVWGNHGESMVADLTQATFTKEGKTQKVVDVLDHDYVFDTFFKKI
            GHRAWDILEHRGFTSAASPTKAAIQHMKAWLFGTAPGEVLSMGIPVPEGNPYGIKPGVVFSFPCNVDKEGKIHVV
            EGFKVNDWLREKLDFTEKDLFHEKEIALNHLAQGG*

        You can find more details via https://salilab.org/modeller/tutorial/basic.html
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
                         check_output = [self.target_pir_path])

        self.run.info("Target alignment file", os.path.basename(self.target_pir_path))


    def run_command(self, command, script_name, check_output=None, rename_log=True):
        """Base routine for running MODELLER scripts

        Parameters
        ==========
        command : list of strs
            E.g. ['mod921', 'test_script.py', 'input1', 'input2'] corresponds to the command line
            "mod9.21 test_script.py input1 input2"

        script_name : str
            E.g. 'test_script.py'

        check_output : list, None
            Verify that this list of filepaths exist after the command is ran

        rename_log : bool, True
            MODELLER outputs a log that is renamed to reflect the command and gene used
        """

        # first things first, we CD into MODELLER's directory
        os.chdir(self.directory)

        # try and execute the command
        process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        output, error = process.communicate()

        # if MODELLER script gave a traceback, it is caught here and everything is stopped
        if process.returncode:
            error = error.decode('utf-8').strip()

            error = "\n" + "\n".join(error.split('\n'))
            print(terminal.c(error, color='red'))

            if not self.check_db_only:
                self.out["structure_exists"] = False

            raise ModellerScriptError("The MODELLER script {} did not execute properly. Hopefully it is clear "
                                      "from the above error message. No structure is going to be modelled."\
                                       .format(script_name))

        # If we made it this far, the MODELLER script ran to completion. Now check outputs exist
        if check_output:
            for output in check_output:
                utils.filesnpaths.is_file_exists(output)

        # MODELLER outputs a log that we rename right here, right now
        old_log_name = os.path.splitext(script_name)[0] + ".log"
        if rename_log:
            new_log_name = "gene_{}_{}".format(self.corresponding_gene_call, old_log_name)
            os.rename(old_log_name, new_log_name)
        else:
            new_log_name = old_log_name

        # add to logs
        self.logs[script_name] = new_log_name

        self.run.info("Log of {}".format(script_name), new_log_name)

        # last things last, we CD back into the starting directory
        os.chdir(self.start_dir)


    class EndModeller(Exception):
        pass


def check_MODELLER(executable=None):
    """Test if MODELLER is going to work.

    Checks the executable exists, that a license exists, and can produce the expected output of a
    modeller executable. Exists outside of the class MODELLER so it does not have to be checked
    everytime the class is initialized. 

    Parameters
    ==========
    executable : str, None
        The string representation of a binary MODELLER program. E.g "mod9.21". If None,
        up_to_date_modeller_exec is chosen and tested.

    Returns
    =======
    executable : str
        Returns the executable that you _should_ use, which is not necessarily what is input
    """

    executable = executable if executable else up_to_date_modeller_exec

    scripts_folder = J(os.path.dirname(anvio.__file__), 'data/misc/MODELLER/scripts')
    if utils.filesnpaths.is_dir_empty(scripts_folder):
        raise ConfigError("Anvi'o houses all its MODELLER scripts in %s, but your directory "
                          "contains no scripts. Why you did dat?" % scripts_folder)

    try:
        is_program_exists(executable)
    except ConfigError:
        *prefix, sub_version = up_to_date_modeller_exec.split('.')
        prefix, sub_version = ''.join(prefix), int(sub_version)
        for alternate_version in reversed(range(sub_version - 10, sub_version + 10)):
            alternate_program = prefix + '.' + str(alternate_version)
            if is_program_exists(alternate_program, dont_raise=True):
                executable = alternate_program
                break
        else:
            raise ConfigError("Anvi'o needs a MODELLER program to be installed on your system. You didn't specify one "
                              "(which can be done with `--modeller-executable`), so anvi'o tried the most recent version "
                              "it knows about: '%s'. If you are certain you have it on your system (for instance you can run it "
                              "by typing '%s' in your terminal window), you may want to send a detailed bug report. If you "
                              "don't have it on your system, check out these installation instructions on our website: "
                              "http://merenlab.org/2016/06/18/installing-third-party-software/#modeller" % (executable, executable))

    temp_dir = filesnpaths.get_temp_directory_path()
    shutil.copy2(J(scripts_folder, 'fasta_to_pir.py'), temp_dir)

    test_script = J(temp_dir, 'fasta_to_pir.py')
    test_input = J(os.path.dirname(anvio.__file__), 'tests/sandbox/mock_data_for_structure/proteins.fa')
    test_output = J(temp_dir, 'test_out')

    command = [executable,
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
            raise ConfigError("You're making progress and anvi'o is proud of you! You just need to validate your MODELLER "
                              "with a license key (it's free). Please go to https://salilab.org/modeller/registration.html "
                              "to register for a new license. After you receive an e-mail with your key, please open '%s' "
                              "and replace the characters XXXXX with your own key. Save the file and try again. " % license_target_file)

        else:
            error = "\n" + "\n".join(error.split('\n'))
            print(terminal.c(error, color='red'))
            raise ConfigError("The executable you requested is called `%s`, but anvi'o doesn't agree with you that "
                              "it is a working MODELLER program. That was determined by running the command `%s`, which raised the "
                              "error seen above. If you want to specify a specific MODELLER program, you can specify it with "
                              "`--modeller-executable`." % (executable, " ".join(command)))

    # no error was raised. now check if output file exists
    try:
        filesnpaths.is_file_exists(test_output)
    except FilesNPathsError:
        raise ConfigError("The executable you requested is called `%s`, but anvi'o doesn't agree with you that "
                          "it is a working MODELLER program. That was determined by running the command `%s`, which did not "
                          "output the file expected. If you want to specify a specific MODELLER program, you can specify it with "
                          "`--modeller-executable`." % (executable, " ".join(command)))

    return executable



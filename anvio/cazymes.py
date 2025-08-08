#!/usr/bin/env python
# -*- coding: utf-8
"""This file contains CAZyme related classes."""

import os
import glob
import shutil

import anvio
import anvio.dbops as dbops
import anvio.utils as utils
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths

from anvio.errors import ConfigError
from anvio.drivers.hmmer import HMMer
from anvio.dbops import ContigsDatabase
from anvio.parsers import parser_modules
from anvio.tables.genefunctions import TableForGeneFunctions

__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "Matthew Schechter"
__email__ = "mschechter@uchicago.edu"

run = terminal.Run()
progress = terminal.Progress()
pp = terminal.pretty_print

class CAZymeSetup(object):
    def __init__(self, args, run=run, progress=progress):
        """Setup a CAZyme database for anvi'o

        http://www.cazy.org/

        Parameters
        ==========
        args : argparse.Namespace
            See `bin/anvi-setup-cazymes` for available arguments
            - cazyme_data_dir : str, optional
                The directory where the CAZyme data should be stored. If not provided, the data will be stored in the anvi'o data directory.
            - reset : bool, optional
                If True, the data directory will be deleted and recreated. Defaults to False.
            - cazyme_version : str, optional
                Version of CAZyme database that will be installed
        run : terminal.Run, optional
            An object for printing messages to the console.
        progress : terminal.Progress, optional
            An object for printing progress bars to the console.
        """

        self.args = args
        self.run = run
        self.progress = progress
        self.cazyme_data_dir = args.cazyme_data_dir

        filesnpaths.is_program_exists('hmmpress')

        if self.cazyme_data_dir and args.reset:
            raise ConfigError("You are attempting to run CAZyme setup on a non-default data directory (%s) using the --reset flag. "
                              "To avoid automatically deleting a directory that may be important to you, anvi'o refuses to reset "
                              "directories that have been specified with --cazyme-data-dir. If you really want to get rid of this "
                              "directory and regenerate it with CAZyme data inside, then please remove the directory yourself using "
                              "a command like `rm -r %s`. We are sorry to make you go through this extra trouble, but it really is "
                              "the safest way to handle things." % (self.cazyme_data_dir, self.cazyme_data_dir))

        if not self.cazyme_data_dir:
            self.cazyme_data_dir = os.path.join(os.path.dirname(anvio.__file__), 'data/misc/CAZyme')

        filesnpaths.is_output_dir_writable(os.path.dirname(os.path.abspath(self.cazyme_data_dir)))

        self.resolve_db_url()

        if not args.reset and not anvio.DEBUG:
            self.is_database_exists()

        if args.reset:
            filesnpaths.gen_output_directory(self.cazyme_data_dir, delete_if_exists=True, dont_warn=True)
        else:
            filesnpaths.gen_output_directory(self.cazyme_data_dir)


    def resolve_db_url(self):
        """Create path to CAZyme ftp

        Added self values
        =================
        - self.db_version : string
            version of CAZyme database

        """
        if self.args.cazyme_version:
            self.db_version = self.args.cazyme_version.upper()
        else:
            self.db_version = 'V13'

        self.db_url = os.path.join("https://bcb.unl.edu/dbCAN2/download/Databases", f"{self.db_version}", f"dbCAN-HMMdb-{self.db_version}.txt")


    def is_database_exists(self):
        """Determine if CAZyme database has already been downloaded"""
        if os.path.exists(os.path.join(self.cazyme_data_dir, "CAZyme_HMMs.txt")):
            # Try to read the version information if available
            version_file = os.path.join(self.cazyme_data_dir, "version.txt")
            with open(version_file, 'r') as f:
                existing_version = f.readline().strip()

            raise ConfigError(f"It seems you already have CAZyme database {existing_version} installed in the following directory: {self.cazyme_data_dir}. "
                              f"You can use the flag `--reset` to instruct anvi'o to set everything up from scratch at that location, "
                              f"or you can use the parameter `--cazyme-data-dir` to setup the CAZyme databases in a different directory.")


    def download(self, hmmpress_files=True):
        """Download CAZyme database and compress with hmmpress"""

        self.run.info("CAZyme database version", self.db_version)
        self.run.info("Remote database location", self.db_url)
        self.run.info("Local database path", self.cazyme_data_dir)

        try:
            utils.download_file(self.db_url, os.path.join(self.cazyme_data_dir, "CAZyme_HMMs.txt"), progress=self.progress)
        except Exception as e:
            raise ConfigError(f"Anvi'o failed to setup your CAZymes at the stage of downloading the data :/ It is very likely that "
                              f"the version '{self.db_version}' does not exist on the server, or the locations of the files have "
                              f"changed. If you are certain that your internet connection works well otherwise, please let "
                              f"the anvi'o developers know about this issue, and they will find a solution. Here is the error message "
                              f"that came from the depths of the code for your reference: '{e.clear_text()}'.")

        with open(os.path.join(self.cazyme_data_dir, "version.txt"), 'w') as f:
            f.write(f"{self.db_version}\n")

        with open(os.path.join(self.cazyme_data_dir, "db_url.txt"), 'w') as f:
            f.write(f"{self.db_url}\n")

        if hmmpress_files:
            self.hmmpress_files()

        self.run.info_single(f"The CAZyme database {self.db_version} is successfully setup on your computer for anvi'o to use 🎉 Now you "
                             f"can use the program `anvi-run-cazyme` on any contigs-db file to annotate genes in them with CAZymes.",
                             nl_before=1, nl_after=1, mc='green')


    def hmmpress_files(self):
        """Runs hmmpress on CAZyme HMM profiles."""

        file_path = os.path.join(self.cazyme_data_dir, "CAZyme_HMMs.txt")
        cmd_line = ['hmmpress', file_path]
        log_file_path = os.path.join(self.cazyme_data_dir, '00_hmmpress_log.txt')
        ret_val = utils.run_command(cmd_line, log_file_path)

        if ret_val:
            raise ConfigError("Hmm. There was an error while running `hmmpress` on the CAZyme HMM profiles. "
                                "Check out the log file ('%s') to see what went wrong." % (log_file_path))
        else:
            # getting rid of the log file because hmmpress was successful
            os.remove(log_file_path)


class CAZyme(object):
    """Search CAZyme database over contigs-db

    Parameters
    ==========
    args : argparse.Namespace
        See `bin/anvi-run-cazymes` for available arguments
        - contigs_db_path : str
            Path to contigs-db
        - num_threads : str
            Number of threads to run hmmsearch or hmmscan
        - noise_cutoff_terms : str, optional
            Filtering option for HMM search
        - hmm_program : str, optional
            hmmsearch (default) or hmmscan
        - cazyme_data_dir : str, optional
            path to downloaded copy of CAZyme data, default downloaded in anvio
    run : terminal.Run, optional
        An object for printing messages to the console.
    progress : terminal.Progress, optional
        An object for printing progress bars to the console.
    """

    def __init__(self, args, run=run, progress=progress):
        self.run = run
        self.progress = progress

        A = lambda x, t: t(args.__dict__[x]) if x in args.__dict__ else None
        null = lambda x: x
        self.contigs_db_path = A('contigs_db', null)
        self.num_threads = A('num_threads', null)
        self.hmm_program = A('hmmer_program', null) or 'hmmsearch'
        self.noise_cutoff_terms = A('noise_cutoff_terms', null)
        self.cazyme_data_dir = args.cazyme_data_dir
        self.just_do_it = A('just_do_it', null)

        # load_catalog will populate this
        self.function_catalog = {}

        # default noise_cutoff_terms
        if not self.noise_cutoff_terms:
            self.noise_cutoff_terms = "-E 1e-12"

        filesnpaths.is_program_exists(self.hmm_program)
        utils.is_contigs_db(self.contigs_db_path)

        if not self.cazyme_data_dir:
            self.cazyme_data_dir = os.path.join(os.path.dirname(anvio.__file__), 'data/misc/CAZyme')

        self.hmm_file =  os.path.join(self.cazyme_data_dir, "CAZyme_HMMs.txt")
        self.version_txt = os.path.join(self.cazyme_data_dir, "version.txt")
        self.db_url_txt = os.path.join(self.cazyme_data_dir, "db_url.txt")

        self.is_database_exists()

        # Learn about the local database
        self.db_version = open(self.version_txt).readline().strip()
        self.db_url = open(self.db_url_txt).readline().strip()

        # reminder to be a good citizen (Good idea Iva!)
        self.run.warning("Anvi'o will annotate genes in your contigs-db with the dbCAN CAZyme database. "
                         "Please do not forget to cite this database when you publish your science :) "
                         "Here is the best way to cite this resource: http://www.cazy.org/Citing-CAZy",
                         lc='green', header="CITATION")


    def is_database_exists(self):
        """Checks if decompressed database files exist"""

        if not glob.glob(self.cazyme_data_dir):
            raise ConfigError("It seems the CAZyme database is not setup on this system :/ Please run the program "
                              "`anvi-setup-cazymes` to set it up. If you set up CAZyme database at a location that "
                              "is different than the default location using the `--cazyme-data-dir` flag before, "
                              "please makes sure to provide the same path to `anvi-run-cazymes`.")

        # Glob and find what files we have then check if we have them all
        downloaded_files = glob.glob(os.path.join(self.cazyme_data_dir, '*'))

        # here we check if the HMM profile is compressed so we can decompress it for next time
        hmmpress_file_extensions = ["h3f", "h3i", "h3m", "h3p", "txt"]
        extant_extensions = [os.path.basename(file).split(".")[-1] for file in downloaded_files]

        if hmmpress_file_extensions.sort() != extant_extensions.sort() or not os.path.exists(self.version_txt) or not os.path.exists(self.db_url_txt):
            raise ConfigError("Something is not right with the files associated with the CAZyme database. Please run "
                              "'anvi-setup-cazymes --reset' to repeat the setup.")


    def check_hash_in_contigs_db(self):
        """Checks the contigs DB self table to make sure it was not already annotated with CAZymes"""

        contigs_db = ContigsDatabase(self.contigs_db_path)
        current_cazyme_hash_in_contigs_db = contigs_db.db.get_meta_value('cazyme_db_hash', return_none_if_not_in_table=True)
        contigs_db.disconnect()

        if current_cazyme_hash_in_contigs_db and not self.just_do_it:
            raise ConfigError(f"The contigs database {self.contigs_db_path} has already been annotated with CAZyme hits version {self.db_version}. "
                              f"If you really want to overwrite these annotations with new ones, please re-run the command with the flag --just-do-it. "
                              f"For those who need this information, the CAZyme database used to annotate this contigs database previously "
                              f"had the following version/URL hash: {current_cazyme_hash_in_contigs_db}")


    def set_hash_in_contigs_db(self):
        """Modifies the contigs DB self table to indicate which CAZyme database version has been used to annotate it."""

        # Create a hash from the database version and URL for tracking
        hash_content = f"{self.db_version}_{self.db_url}"
        cazyme_hash = utils.get_hash_for_list([hash_content])

        contigs_db = ContigsDatabase(self.contigs_db_path)
        contigs_db.db.set_meta_value('cazyme_db_hash', cazyme_hash)
        contigs_db.disconnect()


    def process(self):
        """Search CAZyme HMMs over contigs-db, parse, and filter results"""

        # safety check for previous annotations so that people don't overwrite those if they don't want to
        self.check_hash_in_contigs_db()

        self.run.info("CAZyme database version", f"{self.db_version} (originally from {self.db_url})" )

        # initialize contigs database
        class Args: pass
        args = Args()
        args.contigs_db = self.contigs_db_path
        run_quiet = terminal.Run(verbose=False)

        contigs_db = dbops.ContigsSuperclass(args, r=run_quiet)
        tmp_directory_path = filesnpaths.get_temp_directory_path()

        # get an instance of gene functions table
        gene_function_calls_table = TableForGeneFunctions(self.contigs_db_path)

        # export AA sequences for genes
        target_files_dict = {'AA:GENE': os.path.join(tmp_directory_path, 'AA_gene_sequences.fa')}
        contigs_db.get_sequences_for_gene_callers_ids(output_file_path=target_files_dict['AA:GENE'],
                                                      simple_headers=True,
                                                      report_aa_sequences=True)

        # run hmmer
        hmmer = HMMer(target_files_dict, num_threads_to_use=self.num_threads, program_to_use=self.hmm_program)
        hmm_hits_file = hmmer.run_hmmer('CAZymes', 'AA', 'GENE', None, None, len(self.function_catalog), self.hmm_file, None, self.noise_cutoff_terms)

        if not hmm_hits_file:
            run.info_single("The HMM search returned no hits :/ So there is nothing to add to the contigs database. But "
                            "now anvi'o will add CAZymes as a functional source with no hits, clean the temporary directories "
                            "and gracefully quit.", nl_before=1, nl_after=1)
            shutil.rmtree(tmp_directory_path)
            hmmer.clean_tmp_dirs()
            gene_function_calls_table.add_empty_sources_to_functional_sources({'CAZyme'})
            return

        # parse hmmer output
        parser = parser_modules['search']['hmmer_table_output'](hmm_hits_file, alphabet='AA', context='GENE', program=self.hmm_program)
        search_results_dict = parser.get_search_results()

        # add functions to database
        functions_dict = {}
        counter = 0
        for hmm_hit in search_results_dict.values():
            accession = hmm_hit['gene_name'].removesuffix('.hmm') # removing the .hmm suffix from CAZyme HMM names
            if hmm_hit['gene_hmm_id'].startswith('PF'): # expaned function string if from PFAM
               function = f"{accession} ({hmm_hit['gene_hmm_id']})"
            else:
               function = hmm_hit['gene_hmm_id']

            functions_dict[counter] = {
                'gene_callers_id': hmm_hit['gene_callers_id'],
                'source': 'CAZyme',
                'accession': accession,
                'function': function,
                'e_value': hmm_hit['e_value']
            }

            counter += 1

        if functions_dict:
            gene_function_calls_table.create(functions_dict)
        else:
            self.run.warning("CAZyme class has no hits to process. Returning empty handed, but still adding CAZyme as "
                             "a functional source.")
            gene_function_calls_table.add_empty_sources_to_functional_sources({'CAZyme'})

        # mark contigs db with hash of cazyme database version for tracking
        self.set_hash_in_contigs_db()

        if anvio.DEBUG:
            run.warning("The temp directories, '%s' and '%s' are kept. Please don't forget to clean those up "
                        "later" % (tmp_directory_path, ', '.join(hmmer.tmp_dirs)), header="Debug")
        else:
            run.info_single('Cleaning up the temp directory (you can use `--debug` if you would '
                            'like to keep it for testing purposes)', nl_before=1, nl_after=1)
            shutil.rmtree(tmp_directory_path)
            hmmer.clean_tmp_dirs()

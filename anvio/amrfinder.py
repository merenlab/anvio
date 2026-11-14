#!/usr/bin/env python
"""Classes to set up and run NCBI's AMRFinderPlus on anvi'o contigs databases."""

import os
import glob
import shutil

import anvio
import anvio.dbops as dbops
import anvio.utils as utils
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths

from anvio.errors import ConfigError
from anvio.tables.genefunctions import TableForGeneFunctions


__copyright__ = "Copyleft 2015-2026, The Anvi'o Project (http://anvio.org/)"
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['jessika-fuessel']


run = terminal.Run()
progress = terminal.Progress()

DEFAULT_AMRFINDER_DATA_DIR = os.path.join(os.path.dirname(anvio.__file__), 'data/misc/AMRFinder')

INSTALL_ADVICE = ("If you want to annotate your genes using AMRFinderPlus, you need to install the program on your "
                  "computer as it is not a part of the standard anvi'o installation. Just run `conda install -c "
                  "conda-forge -c bioconda ncbi-amrfinderplus` in your terminal and try again.")

# each functional annotation source, and the AMRFinderPlus output columns for its accession and function
AMRFINDER_SOURCES = {'AMRFinder': ('Hierarchy node', 'Element name'),
                     'AMRFinder_CLASS': ('Class', 'Class'),
                     'AMRFinder_SUBCLASS': ('Subclass', 'Subclass'),
                     'AMRFinder_SUBTYPE': ('Subtype', 'Subtype'),
                     'AMRFinder_SCOPE': ('Scope', 'Scope'),
                     'AMRFinder_METHOD': ('Method', 'Method')}


class AMRFinderSetup(object):
    """Download the AMRFinderPlus database with NCBI's own `amrfinder_update`.

    Parameters
    ==========
    args : argparse.Namespace
        See `anvio/cli/setup_amrfinder.py` for available arguments
        - amrfinder_data_dir : str, optional
            The directory where the database should be stored. The anvi'o data directory is used if not provided.
        - reset : bool, optional
            If True, the default data directory will be deleted and set up from scratch.
    run : terminal.Run, optional
        An object for printing messages to the console.
    progress : terminal.Progress, optional
        An object for printing progress bars to the console.
    """

    def __init__(self, args, run=run, progress=progress):
        self.run = run
        self.progress = progress

        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        self.amrfinder_data_dir = A('amrfinder_data_dir')
        self.reset = A('reset')

        filesnpaths.is_program_exists('amrfinder_update', advice_if_not_exists=INSTALL_ADVICE)

        if self.amrfinder_data_dir and self.reset:
            raise ConfigError(f"You are attempting to use the `--reset` flag on a data directory you chose yourself "
                              f"('{self.amrfinder_data_dir}'). Anvi'o wants to avoid deleting a directory that may "
                              f"contain important data and the `--reset` flag only works for the default directory. "
                              f"If you want to remove the contents of this directory and re-download AMRFinderPlus's "
                              f"database there you need to remove the directory manually with `rm -r "
                              f"{self.amrfinder_data_dir}`, then you can store the AMRFinderPlus database by creating "
                              f"the directory again with the program `anvi-setup-amrfinder` and the flag "
                              f"`--amrfinder-data-dir`.")

        if not self.amrfinder_data_dir:
            self.amrfinder_data_dir = DEFAULT_AMRFINDER_DATA_DIR


    def setup(self):
        """Run `amrfinder_update` to download the latest AMRFinderPlus database."""

        if os.path.exists(os.path.join(self.amrfinder_data_dir, 'latest')) and not self.reset:
            raise ConfigError(f"It seems there already is an AMRFinderPlus database in '{self.amrfinder_data_dir}'. You "
                              f"can use the flag `--reset` to instruct anvi'o to set everything up from scratch at the "
                              f"default location, or you can use the parameter `--amrfinder-data-dir` to set up the "
                              f"database in a different directory.")

        filesnpaths.gen_output_directory(self.amrfinder_data_dir, delete_if_exists=self.reset, dont_warn=True)

        log_file_path = filesnpaths.get_temp_file_path()

        self.run.info('AMRFinderPlus data directory', self.amrfinder_data_dir)
        self.run.info('Log file', log_file_path)

        self.progress.new('AMRFinderPlus database setup')
        self.progress.update('Running `amrfinder_update` (this may take a few minutes) ...')
        ret_val = utils.run_command(['amrfinder_update', '--database', self.amrfinder_data_dir], log_file_path)
        self.progress.end()

        if ret_val:
            raise ConfigError(f"Something went wrong while `amrfinder_update` was downloading the AMRFinderPlus "
                              f"database :/ The log file may tell you more about what went wrong: '{log_file_path}'.")

        os.remove(log_file_path)

        db_version = open(os.path.join(self.amrfinder_data_dir, 'latest', 'version.txt')).read().strip()

        self.run.info_single(f"Success! AMRFinderPlus database {db_version} is now available on your computer for "
                             f"anvi'o to use. Now you can run `anvi-run-amrfinder` on any contigs-db to annotate "
                             f"genes with AMRFinderPlus.", nl_before=1, nl_after=1, mc='green')


class AMRFinder(object):
    """Annotate genes in a contigs-db with NCBI's AMRFinderPlus.

    Parameters
    ==========
    args : argparse.Namespace
        See `anvio/cli/run_amrfinder.py` for available arguments
        - contigs_db : str
            Path to contigs-db
        - amrfinder_data_dir : str, optional
            Path to the AMRFinderPlus data directory. The anvi'o data directory is used if not provided.
        - num_threads : int, optional
            Number of threads for AMRFinderPlus to use
        - skip_amrfinder_e_values : bool, optional
            If True, do not read e-values from the hmmsearch table of AMRFinderPlus, and store 0 instead.
        - full_amrfinder_report : str, optional
            Output file path for the AMRFinderPlus output table, with HMM e-values, bit scores, and cutoffs
        - just_do_it : bool, optional
            If True, annotate a contigs-db again even if it was annotated with AMRFinderPlus before.
    run : terminal.Run, optional
        An object for printing messages to the console.
    progress : terminal.Progress, optional
        An object for printing progress bars to the console.
    """

    def __init__(self, args, run=run, progress=progress):
        self.args = args
        self.run = run
        self.progress = progress

        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        self.contigs_db_path = A('contigs_db')
        self.amrfinder_data_dir = A('amrfinder_data_dir') or DEFAULT_AMRFINDER_DATA_DIR
        self.num_threads = A('num_threads') or 1
        self.skip_amrfinder_e_values = A('skip_amrfinder_e_values')
        self.full_report_path = A('full_amrfinder_report')
        self.just_do_it = A('just_do_it')

        filesnpaths.is_program_exists('amrfinder', advice_if_not_exists=INSTALL_ADVICE)
        utils.is_contigs_db(self.contigs_db_path)

        if self.full_report_path:
            filesnpaths.is_output_file_writable(self.full_report_path)

        # the data directory can be a single database version, or the directory in which `amrfinder_update` keeps
        # each database version in its own directory, and links the newest one as 'latest'
        if os.path.exists(os.path.join(self.amrfinder_data_dir, 'version.txt')):
            self.database_path = self.amrfinder_data_dir
        else:
            self.database_path = os.path.join(self.amrfinder_data_dir, 'latest')

        version_file_path = os.path.join(self.database_path, 'version.txt')

        if not os.path.exists(version_file_path):
            raise ConfigError(f"There doesn't seem to be an AMRFinderPlus database in '{self.amrfinder_data_dir}' :/ "
                              f"To download the database please run the program `anvi-setup-amrfinder`. If you "
                              f"downloaded the database to another location than the default directory, make sure to "
                              f"provide the correct path to this location using the flag `--amrfinder-data-dir` "
                              f"followed by the directory's path.")

        self.db_version = open(version_file_path).read().strip()

        # reminder to be a good citizen
        self.run.warning("Anvi'o will annotate genes in your contigs-db with NCBI's AMRFinderPlus. Please do not "
                         "forget to cite it when you publish your science :) Here is the best way to cite this "
                         "resource: https://doi.org/10.1038/s41598-021-91456-0", lc='green', header="CITATION")


    def process(self):
        """Run AMRFinderPlus on all genes in the contigs-db, and store its hits as functions."""

        # safety check so that people don't overwrite previous annotations by accident
        contigs_db = dbops.ContigsDatabase(self.contigs_db_path, run=self.run, progress=self.progress)
        previous_db_version = contigs_db.db.get_meta_value('amrfinder_db_version', return_none_if_not_in_table=True)
        contigs_db.disconnect()

        if previous_db_version and not self.just_do_it:
            raise ConfigError(f"Your contigs database '{self.contigs_db_path}' has already been annotated with the "
                              f"AMRFinderPlus database {previous_db_version}. If you would like to replace those "
                              f"annotations with new ones, please run the command again with the flag `--just-do-it`.")

        self.run.info('AMRFinderPlus database', f"{self.db_version} ({self.database_path})")
        self.run.info('Number of threads', self.num_threads)

        temp_dir_path = filesnpaths.get_temp_directory_path()
        amrfinder_output_path = self.run_amrfinder(temp_dir_path)

        expected_fields = ['Protein id', 'HMM accession'] + [field for fields in AMRFINDER_SOURCES.values() for field in fields]
        amrfinder_hits = utils.get_TAB_delimited_file_as_dictionary(amrfinder_output_path, indexing_field=-1, expected_fields=expected_fields)

        hmm_hits = {} if self.skip_amrfinder_e_values else self.get_hmm_hits(temp_dir_path)

        self.store_hits(amrfinder_hits, hmm_hits)

        if self.full_report_path:
            self.write_full_report(amrfinder_output_path, amrfinder_hits, hmm_hits)

        if anvio.DEBUG:
            self.run.warning(f"You included the `--debug` flag and anvi'o will keep the temporary directory "
                             f"'{temp_dir_path}'. Please don't forget to clean it up later.", header="Debug")
        else:
            shutil.rmtree(temp_dir_path)


    def run_amrfinder(self, temp_dir_path):
        """Export the amino acid sequences of all genes in the contigs-db, and run AMRFinderPlus on them.

        Parameters
        ==========
        temp_dir_path : str
            Path to a temporary directory for the input, output, and log files, which also becomes $TMPDIR for AMRFinderPlus

        Returns
        =======
        amrfinder_output_path : str
            Path to the AMRFinderPlus output table
        """

        aa_sequences_path = os.path.join(temp_dir_path, 'AA_gene_sequences.fa')
        amrfinder_output_path = os.path.join(temp_dir_path, 'AMRFinder_output.txt')
        log_file_path = os.path.join(temp_dir_path, 'AMRFinder_log.txt')

        contigs_db = dbops.ContigsSuperclass(self.args, r=self.run, p=self.progress)
        contigs_db.get_sequences_for_gene_callers_ids(gene_caller_ids_list=list(contigs_db.genes_in_contigs_dict.keys()),
                                                      output_file_path=aa_sequences_path,
                                                      report_aa_sequences=True,
                                                      simple_headers=True)

        cmd_line = ['amrfinder',
                    '--protein', aa_sequences_path,
                    '--database', self.database_path,
                    '--plus',
                    '--print_node',
                    '--threads', self.num_threads,
                    '--log', log_file_path,
                    '--output', amrfinder_output_path]

        # AMRFinderPlus creates its temporary directory under $TMPDIR, and keeps it when it runs with `--log`. we
        # point $TMPDIR to our own temporary directory, so we know where to find its hmmsearch table later
        tmpdir_env = os.environ.get('TMPDIR')
        os.environ['TMPDIR'] = temp_dir_path

        self.progress.new('AMRFinderPlus')
        self.progress.update('Searching genes (this may take a while) ...')
        try:
            ret_val = utils.run_command(cmd_line, log_file_path)
        finally:
            self.progress.end()

            if tmpdir_env is None:
                del os.environ['TMPDIR']
            else:
                os.environ['TMPDIR'] = tmpdir_env

        if ret_val:
            raise ConfigError(f"AMRFinderPlus did not finish successfully :/ Its log file may tell you more about "
                              f"what went wrong: '{log_file_path}'.")

        return amrfinder_output_path


    def get_hmm_hits(self, temp_dir_path):
        """Read the e-values and bit scores of HMM hits from the hmmsearch table AMRFinderPlus left behind.

        AMRFinderPlus does not report e-values for HMM hits, but its own hmmsearch table has them.

        Parameters
        ==========
        temp_dir_path : str
            Path to the temporary directory AMRFinderPlus ran in

        Returns
        =======
        hmm_hits : dict
            Maps (gene caller id, HMM accession) to the 'e_value' and 'bit_score' of that hit
        """

        hmmsearch_table_paths = glob.glob(os.path.join(temp_dir_path, 'amrfinder.*', 'hmmsearch'))

        if not hmmsearch_table_paths:
            # the files AMRFinderPlus left behind are of no use without the hmmsearch table
            if not anvio.DEBUG:
                shutil.rmtree(temp_dir_path)

            raise ConfigError("Anvi'o tried to extract the e-values of your HMM hits from the temporary hmmsearch "
                              "table AMRFinderPlus usually leaves behind but could not find it :/ That's not a big "
                              "problem if you can live without these e-values. The AMRFinderPlus database contains "
                              "manually curated thresholds for reliable identification and you shouldn't need to "
                              "worry about that. We don't know why anvi'o couldn't find the file, but maybe there "
                              "was a recent AMRFinderPlus update and the program no longer generates the temporary "
                              "directory with this information. To nevertheless annotate your genes with "
                              "AMRFinderPlus, simply run the program again with the flag `--skip-amrfinder-e-values` to store "
                              "all AMRFinderPlus annotations with an e-value of 0.")

        hmm_hits = {}
        with open(hmmsearch_table_paths[0]) as hmmsearch_table:
            for line in hmmsearch_table:
                if not line.startswith('#'):
                    fields = line.split()
                    hmm_hits[(fields[0], fields[3])] = {'e_value': float(fields[4]), 'bit_score': float(fields[5])}

        return hmm_hits


    def store_hits(self, amrfinder_hits, hmm_hits):
        """Store AMRFinderPlus hits as functions in the contigs-db, and record the database version.

        Parameters
        ==========
        amrfinder_hits : dict
            The AMRFinderPlus output table, one entry per hit
        hmm_hits : dict
            E-values and bit scores of HMM hits (see `get_hmm_hits`)
        """

        functions_dict = {}
        for hit in amrfinder_hits.values():
            hmm_hit = hmm_hits.get((hit['Protein id'], hit['HMM accession']))

            for source, (accession_field, function_field) in AMRFINDER_SOURCES.items():
                # AMRFinderPlus reports 'NA' when there is nothing to report (e.g., there is no class for virulence genes)
                if hit[function_field] == 'NA':
                    continue

                functions_dict[len(functions_dict)] = {'gene_callers_id': int(hit['Protein id']),
                                                       'source': source,
                                                       'accession': hit[accession_field],
                                                       'function': hit[function_field],
                                                       'e_value': hmm_hit['e_value'] if hmm_hit else 0}

        gene_function_calls_table = TableForGeneFunctions(self.contigs_db_path, self.run, self.progress)

        if functions_dict:
            gene_function_calls_table.create(functions_dict)
        else:
            self.run.warning("None of the genes in your contigs-db passed the curated cutoffs of the AMRFinderPlus "
                             "database. Anvi'o will still add the AMRFinderPlus annotation sources to your contigs-db, "
                             "so you know that you used this annotation source already.")

        # sources without hits are registered as empty ones (which also removes their annotations from a previous run)
        sources_without_hits = set(AMRFINDER_SOURCES) - set([entry['source'] for entry in functions_dict.values()])
        if sources_without_hits:
            gene_function_calls_table.add_empty_sources_to_functional_sources(sources_without_hits)

        contigs_db = dbops.ContigsDatabase(self.contigs_db_path, run=self.run, progress=self.progress)
        contigs_db.db.set_meta_value('amrfinder_db_version', self.db_version)
        contigs_db.disconnect()


    def write_full_report(self, amrfinder_output_path, amrfinder_hits, hmm_hits):
        """Write the AMRFinderPlus output table with the HMM e-value, bit score, and curated cutoff of each hit.

        Parameters
        ==========
        amrfinder_output_path : str
            Path to the AMRFinderPlus output table, which gives the column order
        amrfinder_hits : dict
            The AMRFinderPlus output table, one entry per hit
        hmm_hits : dict
            E-values and bit scores of HMM hits (see `get_hmm_hits`)
        """

        # the curated HMM cutoffs (the full sequence bit score thresholds, TC1) are in the AMRFinderPlus database
        families = utils.get_TAB_delimited_file_as_dictionary(os.path.join(self.database_path, 'fam.tsv'))
        hmm_cutoffs = {family['hmm_id']: family['hmm_tc1'] for family in families.values()}

        for hit in amrfinder_hits.values():
            hmm_hit = hmm_hits.get((hit['Protein id'], hit['HMM accession']))
            hit['HMM e-value'] = hmm_hit['e_value'] if hmm_hit else 'NA'
            hit['HMM bit score'] = hmm_hit['bit_score'] if hmm_hit else 'NA'
            hit['HMM cutoff'] = hmm_cutoffs.get(hit['HMM accession'], 'NA')

        headers = utils.get_columns_of_TAB_delim_file(amrfinder_output_path, include_first_column=True) + ['HMM e-value', 'HMM bit score', 'HMM cutoff']
        utils.store_dict_as_TAB_delimited_file(amrfinder_hits, self.full_report_path, headers=['entry'] + headers, do_not_write_key_column=True)
        self.run.info('Full AMRFinderPlus report', self.full_report_path)

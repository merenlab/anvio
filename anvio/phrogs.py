#!/usr/bin/env python
# pylint: disable=line-too-long
"""This file contains PHROGs setup and annotation classes."""

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
from anvio.parsers import parser_modules
from anvio.tables.genefunctions import TableForGeneFunctions


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "The Anvi'o Project"
__email__ = "info@anvio.org"


run = terminal.Run()
progress = terminal.Progress()


class PHROGsSetup(object):
    def __init__(self, args, run=run, progress=progress):
        self.args = args
        self.run = run
        self.progress = progress
        self.phrogs_data_dir = args.phrogs_data_dir
        self.phrogs_version = args.phrogs_version or 'v4'

        filesnpaths.is_program_exists('hmmpress')

        if self.phrogs_data_dir and args.reset:
            raise ConfigError("You are attempting to run PHROGs setup on a non-default data directory (%s) using the --reset flag. "
                              "To avoid automatically deleting a directory that may be important to you, anvi'o refuses to reset "
                              "directories that have been specified with --phrogs-data-dir. If you really want to get rid of this "
                              "directory and regenerate it with PHROGs data inside, then please remove the directory yourself using "
                              "a command like `rm -r %s`. We are sorry to make you go through this extra trouble, but it really is "
                              "the safest way to handle things." % (self.phrogs_data_dir, self.phrogs_data_dir))

        if not self.phrogs_data_dir:
            self.phrogs_data_dir = os.path.join(os.path.dirname(anvio.__file__), 'data/misc/PHROGs')

        filesnpaths.is_output_dir_writable(os.path.dirname(os.path.abspath(self.phrogs_data_dir)))

        self.db_url = "https://phrogs.lmge.uca.fr/downloads_from_website"
        self.hmm_archive_name = "HMM_phrog.tar.gz"
        self.annotation_file_name = f"phrog_annot_{self.phrogs_version}.tsv"
        self.local_hmm_archive_path = os.path.join(self.phrogs_data_dir, self.hmm_archive_name)
        self.local_annotation_path = os.path.join(self.phrogs_data_dir, "phrog_annot.tsv")
        self.local_hmm_file = os.path.join(self.phrogs_data_dir, "PHROGs.hmm")

        if not args.reset and not anvio.DEBUG:
            self.is_database_exists()

        if args.reset:
            filesnpaths.gen_output_directory(self.phrogs_data_dir, delete_if_exists=True, dont_warn=True)
        else:
            filesnpaths.gen_output_directory(self.phrogs_data_dir)


    def is_database_exists(self):
        if os.path.exists(self.local_hmm_file) and os.path.exists(self.local_annotation_path):
            raise ConfigError("It seems you already have PHROGs database installed in '%s', please use --reset flag if you want to re-download it." % self.phrogs_data_dir)


    def download(self):
        hmm_download_url = f"{self.db_url}/{self.hmm_archive_name}"
        annotation_download_url = f"{self.db_url}/{self.annotation_file_name}"
        self.run.info("PHROGs HMM archive URL", hmm_download_url)
        self.run.info("PHROGs annotation URL", annotation_download_url)
        self.run.info("PHROGs data directory", self.phrogs_data_dir)

        try:
            utils.download_file(hmm_download_url, self.local_hmm_archive_path, progress=self.progress, run=self.run)
            utils.download_file(annotation_download_url, self.local_annotation_path, progress=self.progress, run=self.run)
        except Exception as e:
            raise ConfigError("Anvi'o failed to download PHROGs files. If your internet connection is healthy, this may indicate the "
                              "database URLs have changed. Please report this issue to the anvi'o developers. Original error: '%s'." % e)

        extracted_hmms_dir = os.path.join(self.phrogs_data_dir, "HMM_phrog")
        if os.path.exists(extracted_hmms_dir):
            shutil.rmtree(extracted_hmms_dir)

        utils.tar_extract_file(self.local_hmm_archive_path, output_file_path=self.phrogs_data_dir, keep_original=True)

        hmm_files = sorted(glob.glob(os.path.join(extracted_hmms_dir, "*.hmm")))
        if not hmm_files:
            raise ConfigError("Anvi'o expected to find .hmm files under '%s' after extracting the PHROGs archive, but found none. "
                              "The archive format may have changed." % extracted_hmms_dir)

        utils.concatenate_files(self.local_hmm_file, hmm_files, remove_concatenated_files=False)
        self.hmmpress_profiles()

        with open(os.path.join(self.phrogs_data_dir, "version.txt"), 'w') as f:
            f.write(f"{self.phrogs_version}\n")

        with open(os.path.join(self.phrogs_data_dir, "db_url.txt"), 'w') as f:
            f.write(f"{self.db_url}\n")

        self.run.info_single("The PHROGs database is successfully setup on your computer for anvi'o to use 🎉 Now you can run "
                             "`anvi-run-phrogs` on contigs databases to annotate genes with PHROGs functions.",
                             nl_before=1, nl_after=1, mc='green')


    def hmmpress_profiles(self):
        cmd_line = ['hmmpress', self.local_hmm_file]
        log_file_path = os.path.join(self.phrogs_data_dir, '00_hmmpress_log.txt')
        ret_val = utils.run_command(cmd_line, log_file_path)

        if ret_val:
            raise ConfigError("There was an error while running `hmmpress` on PHROGs HMM profiles. "
                              "Check out the log file ('%s') to see what went wrong." % log_file_path)

        os.remove(log_file_path)


class PHROGs(object):
    def __init__(self, args, run=run, progress=progress):
        self.run = run
        self.progress = progress

        A = lambda x, t: t(args.__dict__[x]) if x in args.__dict__ else None
        null = lambda x: x
        self.contigs_db_path = A('contigs_db', null)
        self.num_threads = A('num_threads', null)
        self.hmm_program = A('hmmer_program', null) or 'hmmsearch'
        self.noise_cutoff_terms = A('noise_cutoff_terms', null) or '--cut_ga'
        self.phrogs_data_dir = A('phrogs_data_dir', null)
        self.just_do_it = A('just_do_it', null)

        filesnpaths.is_program_exists(self.hmm_program)
        utils.is_contigs_db(self.contigs_db_path)

        if not self.phrogs_data_dir:
            self.phrogs_data_dir = os.path.join(os.path.dirname(anvio.__file__), 'data/misc/PHROGs')

        self.hmm_file = os.path.join(self.phrogs_data_dir, "PHROGs.hmm")
        self.annotation_file = os.path.join(self.phrogs_data_dir, "phrog_annot.tsv")
        self.version_file = os.path.join(self.phrogs_data_dir, "version.txt")
        self.db_url_file = os.path.join(self.phrogs_data_dir, "db_url.txt")
        self.function_catalog = {}

        self.is_database_exists()
        self.load_catalog()


    def is_database_exists(self):
        expected_files = [self.hmm_file, self.annotation_file]
        for file_path in expected_files:
            if not os.path.exists(file_path):
                raise ConfigError("It seems the PHROGs database is not setup on this system :/ Please run the program "
                                  "`anvi-setup-phrogs` to set it up. If you set up PHROGs at a location that is different "
                                  "than the default location using the `--phrogs-data-dir` flag before, please provide "
                                  "the same path to `anvi-run-phrogs`.")


    def canonicalize_phrog_id(self, raw_id):
        if raw_id is None:
            return None

        phrog_id = str(raw_id).strip()
        if not phrog_id:
            return None

        if phrog_id.lower().startswith('phrog_'):
            phrog_suffix = phrog_id.split('_', maxsplit=1)[1]
        elif phrog_id.lower().startswith('phrog'):
            phrog_suffix = phrog_id[5:]
        else:
            phrog_suffix = phrog_id

        phrog_suffix = phrog_suffix.strip()
        if phrog_suffix.isdigit():
            return f"phrog_{int(phrog_suffix):04d}"
        else:
            return f"phrog_{phrog_suffix.lower()}"


    def load_catalog(self):
        rows_dict = utils.get_TAB_delimited_file_as_dictionary(self.annotation_file)

        for row in rows_dict.values():
            phrog_id = row.get('#phrog') or row.get('phrog') or row.get('PHROG') or row.get('phrog_id')
            annotation = row.get('Annotation') or row.get('annotation') or row.get('function') or row.get('description')
            category = row.get('Category') or row.get('category')

            canonical_phrog_id = self.canonicalize_phrog_id(phrog_id)
            if not canonical_phrog_id:
                continue

            if annotation and category:
                function = f"{annotation} [{category}]"
            elif annotation:
                function = annotation
            elif category:
                function = f"PHROG function category: {category}"
            else:
                function = f"Unannotated PHROG {canonical_phrog_id}"

            self.function_catalog[canonical_phrog_id] = function

        if not self.function_catalog:
            raise ConfigError("The PHROGs annotation file '%s' could not be parsed into a usable catalog. "
                              "Please make sure this file is a valid PHROGs annotation table." % self.annotation_file)

        self.run.info("PHROGs annotation entries loaded", len(self.function_catalog))
        if os.path.exists(self.version_file):
            self.run.info("PHROGs version", open(self.version_file).readline().strip())


    def process(self):
        if not self.just_do_it:
            self.run.warning("Anvi'o will annotate genes in your contigs-db with PHROGs. If your contigs-db already "
                             "contains PHROGs hits and you intend to overwrite them, rerun this command with --just-do-it.",
                             header="Heads up", lc='green')

        class Args:
            pass

        args = Args()
        args.contigs_db = self.contigs_db_path
        run_quiet = terminal.Run(verbose=False)

        contigs_db = dbops.ContigsSuperclass(args, r=run_quiet)
        tmp_directory_path = filesnpaths.get_temp_directory_path()
        gene_function_calls_table = TableForGeneFunctions(self.contigs_db_path)

        target_files_dict = {'AA:GENE': os.path.join(tmp_directory_path, 'AA_gene_sequences.fa')}
        contigs_db.get_sequences_for_gene_callers_ids(output_file_path=target_files_dict['AA:GENE'],
                                                      simple_headers=True,
                                                      report_aa_sequences=True)

        hmmer = HMMer(target_files_dict, num_threads_to_use=self.num_threads, program_to_use=self.hmm_program)
        hmm_hits_file = hmmer.run_hmmer('PHROGs', 'AA', 'GENE', None, None, len(self.function_catalog), self.hmm_file, None, self.noise_cutoff_terms)

        if not hmm_hits_file:
            self.run.info_single("The HMM search returned no hits. Anvi'o will add PHROGs as a functional source with no hits "
                                 "and exit gracefully.", nl_before=1, nl_after=1)
            shutil.rmtree(tmp_directory_path)
            hmmer.clean_tmp_dirs()
            gene_function_calls_table.add_empty_sources_to_functional_sources({'PHROGs'})
            return

        parser = parser_modules['search']['hmmer_table_output'](hmm_hits_file, alphabet='AA', context='GENE', program=self.hmm_program)
        search_results_dict = parser.get_search_results()

        functions_dict = {}
        counter = 0
        for hmm_hit in search_results_dict.values():
            canonical_id = self.canonicalize_phrog_id(hmm_hit.get('gene_hmm_id') or hmm_hit.get('gene_name'))
            functions_dict[counter] = {
                'gene_callers_id': hmm_hit['gene_callers_id'],
                'source': 'PHROGs',
                'accession': canonical_id,
                'function': self.function_catalog.get(canonical_id, f"Unknown PHROG hit ({canonical_id})"),
                'e_value': hmm_hit['e_value']
            }
            counter += 1

        if functions_dict:
            gene_function_calls_table.create(functions_dict)
        else:
            self.run.warning("PHROGs has no hits to process. Returning empty handed, but still adding PHROGs as a functional source.")
            gene_function_calls_table.add_empty_sources_to_functional_sources({'PHROGs'})

        if anvio.DEBUG:
            self.run.warning("The temp directories, '%s' and '%s' are kept. Please don't forget to clean those up "
                             "later" % (tmp_directory_path, ', '.join(hmmer.tmp_dirs)), header="Debug")
        else:
            self.run.info_single('Cleaning up the temp directory (you can use `--debug` if you would '
                                 'like to keep it for testing purposes)', nl_before=1, nl_after=1)
            shutil.rmtree(tmp_directory_path)
            hmmer.clean_tmp_dirs()

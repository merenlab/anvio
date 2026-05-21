#!/usr/bin/env python
# pylint: disable=line-too-long
"""This file contains PHROGs setup and annotation classes."""

import os
import glob
import shutil
import argparse

import anvio
import anvio.dbops as dbops
import anvio.utils as utils
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths

from anvio.errors import ConfigError
from anvio.drivers.hmmer import HMMer
from anvio.parsers import parser_modules
from anvio.tables.genefunctions import TableForGeneFunctions

__copyright__ = "Copyleft 2015-2026, The Anvi'o Project (http://anvio.org/)"
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "Georges Kanaan"
__email__ = "georges@gkanaan.com"


class PHROGsSetup(object):
    def __init__(self, args, run=terminal.Run(), progress=terminal.Progress()):
        self.args = args
        self.run = run
        self.progress = progress
        self.phrogs_data_dir = args.phrogs_data_dir
        self.phrogs_version = args.phrogs_version or 'v4'

        filesnpaths.is_program_exists('hmmpress')
        filesnpaths.is_program_exists('hmmbuild')

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
        self.msa_archive_name = "MSA_phrogs.tar.gz"
        self.annotation_file_name = f"phrog_annot_{self.phrogs_version}.tsv"

        self.local_msa_archive_path = os.path.join(self.phrogs_data_dir, self.msa_archive_name)
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
        msa_download_url = f"{self.db_url}/{self.msa_archive_name}"
        annotation_download_url = f"{self.db_url}/{self.annotation_file_name}"

        self.run.info("PHROGs MSA archive URL", msa_download_url)
        self.run.info("PHROGs annotation URL", annotation_download_url)
        self.run.info("PHROGs data directory", self.phrogs_data_dir)

        try:
            utils.download_file(
                msa_download_url,
                self.local_msa_archive_path,
                check_certificate=False,
                progress=self.progress,
                run=self.run,
            )
            utils.download_file(
                annotation_download_url,
                self.local_annotation_path,
                check_certificate=False,
                progress=self.progress,
                run=self.run,
            )
        except Exception as e:
            raise ConfigError("Anvi'o failed to download PHROGs files. If your internet connection is healthy, this may indicate the "
                              "database URLs have changed. Please report this issue to the anvi'o developers. Original error: '%s'." % e)

        extracted_fma_dir = os.path.join(self.phrogs_data_dir, "MSA_Phrogs_M50_FMA")
        if os.path.exists(extracted_fma_dir):
            shutil.rmtree(extracted_fma_dir)

        utils.tar_extract_file(self.local_msa_archive_path, output_file_path=self.phrogs_data_dir, keep_original=True)

        fma_files = sorted(glob.glob(os.path.join(extracted_fma_dir, "*.fma")))

        if not fma_files:
            raise ConfigError("Anvi'o expected to find FMA files under '%s' after extracting the PHROGs archive, but found none. "
                              "The archive format may have changed." % extracted_fma_dir)

        hmm_files = self.build_hmms_from_fmas(fma_files)

        utils.concatenate_files(self.local_hmm_file, hmm_files, remove_concatenated_files=True)
        self.hmmpress_profiles()
        
        with open(os.path.join(self.phrogs_data_dir, "version.txt"), 'w') as f:
            f.write(f"{self.phrogs_version}\n")

        with open(os.path.join(self.phrogs_data_dir, "db_url.txt"), 'w') as f:
            f.write(f"{self.db_url}\n")

        self.run.info_single("The PHROGs database is successfully setup on your computer for anvi'o to use 🎉 Now you can run "
                             "`anvi-run-phrogs` on contigs databases to annotate genes with PHROGs functions.",
                             nl_before=1, nl_after=1, mc='green')

    def build_hmms_from_fmas(self, fma_files):
        log_file_path = os.path.join(self.phrogs_data_dir, '00_hmmbuild_log.txt')
        fma_file = ""
        hmm_file = ""
        hmm_files = []

        try:
            self.progress.new("Building PHROGs .hmm files from FMA families", progress_total_items=len(fma_files))

            for fma_file in fma_files:
                hmm_file = os.path.splitext(fma_file)[0] + '.hmm'
                utils.run_command_and_get_output(
                    ['hmmbuild', hmm_file, fma_file],
                    log_file_path=log_file_path
                )
                hmm_files.append(hmm_file)
                self.progress.increment()

            self.progress.end()

        except Exception:
            if os.path.exists(hmm_file):
                os.remove(hmm_file)
            self.progress.end()
            raise ConfigError("Anvi'o expected to be able to build PHROGs .hmm files from FMA families using `hmmbuild`, "
                              "but something went wrong. Please check the log file ('%s') to see what happened." % log_file_path)

        # Delete FMAs
        for fma_file in fma_files:
            os.remove(fma_file)
        
        return hmm_files

    def hmmpress_profiles(self):
        cmd_line = ['hmmpress', self.local_hmm_file]
        log_file_path = os.path.join(self.phrogs_data_dir, '00_hmmpress_log.txt')
        ret_val = utils.run_command(cmd_line, log_file_path)

        if ret_val:
            raise ConfigError("There was an error while running `hmmpress` on PHROGs HMM profiles. "
                              "Check out the log file ('%s') to see what went wrong." % log_file_path)

        os.remove(log_file_path)


class PHROGs(object):
    def __init__(self, args, run=terminal.Run(), progress=terminal.Progress()):
        self.args = args
        self.run = run
        self.progress = progress

        A = lambda x, t: t(args.__dict__[x]) if x in args.__dict__ else None
        null = lambda x: x

        self.contigs_db_path = A('contigs_db', null)
        self.num_threads = A('num_threads', null) or 1
        self.hmm_program = A('hmmer_program', null) or 'hmmsearch'
        self.phrogs_data_dir = A('phrogs_data_dir', null)
        self.noise_cutoff_terms = A('noise_cutoff_terms', null) or '--cut_ga'
        self.just_do_it = A('just_do_it', null)

        filesnpaths.is_program_exists(self.hmm_program)
        utils.is_contigs_db(self.contigs_db_path)

        if not self.phrogs_data_dir:
            self.phrogs_data_dir = os.path.join(os.path.dirname(anvio.__file__), 'data/misc/PHROGs')

        self.local_hmm_file = os.path.join(self.phrogs_data_dir, 'PHROGs.hmm')
        self.local_annotation_path = os.path.join(self.phrogs_data_dir, 'phrog_annot.tsv')

        self.is_database_exists()
        self.load_catalog()

    def is_database_exists(self):
        if not os.path.exists(self.local_hmm_file):
            raise ConfigError("It seems you do not have PHROGs HMM profiles installed in '%s'. Please run 'anvi-setup-phrogs' first." %
                              self.phrogs_data_dir)

        if not os.path.exists(self.local_annotation_path):
            raise ConfigError("It seems you do not have PHROGs annotations installed in '%s'. Please run 'anvi-setup-phrogs' first." %
                              self.phrogs_data_dir)

    def load_catalog(self):
        self.function_catalog = utils.get_TAB_delimited_file_as_dictionary(
            self.local_annotation_path,
            column_names=['phrog', 'color', 'annot', 'category'],
        )

    def get_function_from_catalog(self, accession, ok_if_missing_from_catalog=False):
        accession = str(accession)

        if accession not in self.function_catalog:
            if ok_if_missing_from_catalog:
                return "Unknown function with PHROGs accession %s" % accession
            raise ConfigError("It seems hmmscan/hmmsearch found an accession id that does not exist in the PHROGs catalog: %s" % accession)

        return self.function_catalog[accession]['annot']

    def process(self):
        args = argparse.Namespace(contigs_db=self.contigs_db_path)
        contigs_db = dbops.ContigsSuperclass(args)
        tmp_directory_path = filesnpaths.get_temp_directory_path()

        gene_function_calls_table = TableForGeneFunctions(self.contigs_db_path, self.run, self.progress)

        target_files_dict = {'AA:GENE': os.path.join(tmp_directory_path, 'AA_gene_sequences.fa')}
        contigs_db.get_sequences_for_gene_callers_ids(output_file_path=target_files_dict['AA:GENE'],
                                                      simple_headers=True,
                                                      report_aa_sequences=True)

        hmmer = HMMer(target_files_dict, num_threads_to_use=self.num_threads, program_to_use=self.hmm_program)
        hmm_hits_file = hmmer.run_hmmer('PHROGs', 'AA', 'GENE', None, None, len(self.function_catalog), self.local_hmm_file,
                                        None, self.noise_cutoff_terms)

        if not hmm_hits_file:
            self.run.info_single("The HMM search returned no hits :/ So there is nothing to add to the contigs database. But "
                                 "now anvi'o will add PHROGs as a functional source with no hits, clean the temporary directories "
                                 "and gracefully quit.", nl_before=1, nl_after=1)
            shutil.rmtree(tmp_directory_path)
            hmmer.clean_tmp_dirs()
            gene_function_calls_table.add_empty_sources_to_functional_sources({'PHROGs'})
            return

        parser = parser_modules['search']['hmmer_table_output'](hmm_hits_file, alphabet='AA', context='GENE', program=self.hmm_program)
        search_results_dict = parser.get_search_results()

        functions_dict = {}
        counter = 0
        for hmm_hit in search_results_dict.values():
            functions_dict[counter] = {
                'gene_callers_id': hmm_hit['gene_callers_id'],
                'source': 'PHROGs',
                'accession': hmm_hit['gene_hmm_id'],
                'function': self.get_function_from_catalog(hmm_hit['gene_hmm_id'], ok_if_missing_from_catalog=True),
                'e_value': hmm_hit['e_value'],
            }
            counter += 1

        if functions_dict:
            gene_function_calls_table.create(functions_dict, drop_previous_annotations_first=bool(self.just_do_it))
        else:
            self.run.warning("PHROGs class has no hits to process. Returning empty handed, but still adding PHROGs as a functional source.")
            gene_function_calls_table.add_empty_sources_to_functional_sources({'PHROGs'})

        if anvio.DEBUG:
            self.run.warning("The temp directories, '%s' and '%s' are kept. Please don't forget to clean those up later" %
                             (tmp_directory_path, ', '.join(hmmer.tmp_dirs)), header="Debug")
        else:
            self.run.info_single('Cleaning up the temp directory (you can use `--debug` if you would like to keep it for testing purposes)',
                                 nl_before=1, nl_after=1)
            shutil.rmtree(tmp_directory_path)
            hmmer.clean_tmp_dirs()

      
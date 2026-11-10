#!/usr/bin/env python
# pylint: disable=line-too-long
"""This file contains PHROGs setup and annotation classes."""

import os
import glob
import shutil
import argparse
import multiprocess as multiprocessing

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

        # reminder to be a good citizen
        self.run.warning("Anvi'o will set up the PHROGs database, as described in Terzian, Olo Ndela et al. (2021) "
                         "(doi:10.1093/nargab/lqab067). When you publish your findings, please do not forget to "
                         "properly credit this work.", lc='green', header="CITATION")

        self.phrogs_data_dir = args.phrogs_data_dir
        self.phrogs_version = args.phrogs_version or 'v4'
        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        self.num_threads = A('num_threads') or 1

        filesnpaths.is_program_exists('hmmpress')
        filesnpaths.is_program_exists('hmmbuild')

        if self.phrogs_data_dir and args.reset and os.path.exists(self.phrogs_data_dir):
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
        self.msa_extracted_dir_name = "MSA_Phrogs_M50_FASTA"
        self.annotation_file_name = f"phrog_annot_{self.phrogs_version}.tsv"

        self.local_msa_archive_path = os.path.join(self.phrogs_data_dir, self.msa_archive_name)
        self.local_extracted_fma_dir = os.path.join(self.phrogs_data_dir, self.msa_extracted_dir_name)
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

        if os.path.exists(self.local_extracted_fma_dir):
            shutil.rmtree(self.local_extracted_fma_dir)

        utils.tar_extract_file(self.local_msa_archive_path, output_file_path=self.phrogs_data_dir, keep_original=True)

        fma_files = sorted(glob.glob(os.path.join(self.local_extracted_fma_dir, "*.fma")))

        if not fma_files:
            raise ConfigError("Anvi'o expected to find FMA files under '%s' after extracting the PHROGs archive, but found none. "
                              "The archive format may have changed." % self.local_extracted_fma_dir)

        hmm_files = self.build_hmms_from_fmas(fma_files)

        utils.concatenate_files(self.local_hmm_file, hmm_files, remove_concatenated_files=not anvio.DEBUG)
        self.hmmpress_profiles()

        # clean up after successful setup
        if not anvio.DEBUG:
            os.remove(self.local_msa_archive_path)
            shutil.rmtree(self.local_extracted_fma_dir)
        else:
            self.run.warning(f"Since you used the `--debug` flag, anvi'o kept all the intermediate files from this setup. "
                             f"These include the downloaded archive at '{self.local_msa_archive_path}', and the folder "
                             f"'{self.local_extracted_fma_dir}' with all {len(fma_files)} FMA alignments and the individual "
                             f".hmm files built from them, which together take up quite a bit of space. Please don't forget "
                             f"to clean those up later.", header="Debug")

        with open(os.path.join(self.phrogs_data_dir, "version.txt"), 'w') as f:
            f.write(f"{self.phrogs_version}\n")

        with open(os.path.join(self.phrogs_data_dir, "db_url.txt"), 'w') as f:
            f.write(f"{self.db_url}\n")

        self.run.info_single("The PHROGs database is successfully setup on your computer for anvi'o to use 🎉 Now you can run "
                             "`anvi-run-phrogs` on contigs databases to annotate genes with PHROGs functions.",
                             nl_before=1, nl_after=1, mc='green')

    @staticmethod
    def _hmmbuild_worker(fma_files, log_file_path, output_queue):
        """Multiprocessing worker that builds one .hmm file per FMA alignment it is given.

        Parameters
        ==========
        fma_files : list of str
            Paths to the FMA alignment files this worker is responsible for. Each .hmm file is
            written next to its FMA file.
        log_file_path : str
            Path to this worker's own log file, where `hmmbuild` errors are written.
        output_queue : multiprocess.Manager.Queue
            Queue where the worker reports each file as `(True, fma_file, None)` on success, or
            `(False, fma_file, error_message)` on failure. The worker stops after its first failure.
        """

        for fma_file in fma_files:
            hmm_file = os.path.splitext(fma_file)[0] + '.hmm'

            try:
                utils.run_command_and_get_output(['hmmbuild', '--cpu', '1', hmm_file, fma_file], log_file_path=log_file_path)
            except Exception as e:
                if os.path.exists(hmm_file):
                    os.remove(hmm_file)

                output_queue.put((False, fma_file, getattr(e, 'e', None) or str(e)))
                return

            output_queue.put((True, fma_file, None))

    def build_hmms_from_fmas(self, fma_files):
        num_fma_files = len(fma_files)

        # deal the FMA files round-robin into one chunk per worker
        num_workers = min(self.num_threads, num_fma_files)
        chunks = [fma_files[i::num_workers] for i in range(num_workers)]
        log_file_paths = [os.path.join(self.phrogs_data_dir, f'00_hmmbuild_log_worker_{i}.txt') for i in range(num_workers)]

        self.run.info('Number of FMA families', num_fma_files)
        self.run.info('Number of hmmbuild workers', num_workers)

        manager = multiprocessing.Manager()
        output_queue = manager.Queue()

        workers = []
        for chunk, log_file_path in zip(chunks, log_file_paths):
            worker = multiprocessing.Process(target=PHROGsSetup._hmmbuild_worker, args=(chunk, log_file_path, output_queue))
            workers.append(worker)
            worker.start()

        self.progress.new("Building PHROGs .hmm files from FMA families", progress_total_items=num_fma_files)
        self.progress.update(f"0 of {num_fma_files} families")

        num_built = 0
        failure = None
        try:
            while num_built < num_fma_files:
                success, fma_file, error_message = output_queue.get()

                if not success:
                    failure = (fma_file, error_message)
                    break

                num_built += 1
                self.progress.update(f"{num_built} of {num_fma_files} families")
                self.progress.increment(increment_to=num_built)
        finally:
            # always end progress and terminate workers, even on error
            self.progress.end()
            for worker in workers:
                worker.terminate()
            for worker in workers:
                worker.join()

        if failure:
            fma_file, error_message = failure

            # other workers may have been terminated mid-build, so we remove all .hmm files to avoid
            # leaving partial models behind
            for f in fma_files:
                hmm_file = os.path.splitext(f)[0] + '.hmm'
                if os.path.exists(hmm_file):
                    os.remove(hmm_file)

            raise ConfigError(f"Anvi'o expected to be able to build PHROGs .hmm files from FMA families using `hmmbuild`, "
                              f"but something went wrong with '{os.path.basename(fma_file)}'. Here is the error message: "
                              f"\"{error_message}\". You may find more details in the hmmbuild log files at "
                              f"'{os.path.join(self.phrogs_data_dir, '00_hmmbuild_log_worker_*.txt')}'.")

        # we build this list from the input order so the final PHROGs.hmm file is always the same regardless of the number of threads
        hmm_files = [os.path.splitext(f)[0] + '.hmm' for f in fma_files]

        # Delete FMAs
        if not anvio.DEBUG:
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

        if not anvio.DEBUG:
            os.remove(log_file_path)


class PHROGs(object):
    def __init__(self, args, run=terminal.Run(), progress=terminal.Progress()):
        self.args = args
        self.run = run
        self.progress = progress

        # reminder to be a good citizen
        self.run.warning("Anvi'o will annotate your database with the PHROGs database, as described in Terzian, "
                         "Olo Ndela et al. (2021) (doi:10.1093/nargab/lqab067). When you publish your findings, "
                         "please do not forget to properly credit this work.", lc='green', header="CITATION")

        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        self.contigs_db_path = A('contigs_db')
        self.num_threads = A('num_threads') or 1
        self.hmm_program = A('hmmer_program') or 'hmmsearch'
        self.phrogs_data_dir = A('phrogs_data_dir')
        self.noise_cutoff_terms = A('noise_cutoff_terms')

        if self.noise_cutoff_terms is None:
            self.noise_cutoff_terms = '-E 1e-10'
            self.run.info_single(f"No `--noise-cutoff-terms` was provided, so anvi'o will use the default "
                                 f"'{self.noise_cutoff_terms}'. You may wish to choose a different cutoff.",
                                 nl_before=1, nl_after=1, mc='yellow')

        filesnpaths.is_program_exists(self.hmm_program)
        utils.is_contigs_db(self.contigs_db_path)

        if not self.phrogs_data_dir:
            self.phrogs_data_dir = os.path.join(os.path.dirname(anvio.__file__), 'data/misc/PHROGs')

        self.local_hmm_file = os.path.join(self.phrogs_data_dir, 'PHROGs.hmm')
        self.local_annotation_path = os.path.join(self.phrogs_data_dir, 'phrog_annot.tsv')
        self.local_version_file = os.path.join(self.phrogs_data_dir, 'version.txt')

        # names of the functional sources this class adds to the contigs database
        self.function_source = 'PHROG_FUNCTION'
        self.category_source = 'PHROG_CATEGORY'

        self.is_database_exists()
        self.get_version()
        self.load_catalog()

    def get_version(self):
        if os.path.exists(self.local_version_file):
            with open(self.local_version_file) as f:
                version = f.read().strip()
        else:
            version = 'unknown'

        self.run.info('PHROGs version', version)

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
        """Returns a `(function, category)` tuple for a PHROGs accession, from the 'annot' and 'category' columns

        It is not uncommon for PHROGs to have annotation 'NA' and category 'unknown function' in the
        annotation table, and these values will be stored in the contigs database. But when the family
        ID is not in the annotation table at all, this functon returns longer 'Unknown' strings if
        ok_if_missing_from_catalog is True.
        """

        accession = str(accession)

        # the annotation table is keyed by integer PHROG IDs, while model names look like 'phrog_1'
        catalog_key = accession.removeprefix('phrog_')

        if catalog_key not in self.function_catalog:
            if ok_if_missing_from_catalog:
                return (f"Unknown function with PHROGs accession {accession}",
                        f"Unknown category with PHROGs accession {accession}")
            raise ConfigError(f"It seems hmmscan/hmmsearch found an accession id that does not exist in the PHROGs catalog: {accession}")

        entry = self.function_catalog[catalog_key]
        return entry['annot'], entry['category']

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
            self.run.info_single(f"The HMM search returned no hits :/ So there is nothing to add to the contigs database. But "
                                 f"now anvi'o will add {self.function_source} and {self.category_source} as functional sources "
                                 f"with no hits, clean the temporary directories and gracefully quit.", nl_before=1, nl_after=1)
            shutil.rmtree(tmp_directory_path)
            hmmer.clean_tmp_dirs()
            gene_function_calls_table.add_empty_sources_to_functional_sources({self.function_source, self.category_source})
            return

        parser = parser_modules['search']['hmmer_table_output'](hmm_hits_file, alphabet='AA', context='GENE', program=self.hmm_program)
        search_results_dict = parser.get_search_results()

        functions_dict = {}
        counter = 0
        for hmm_hit in search_results_dict.values():
            # models built from PHROGs FMA files have no ACC line, so we use the model NAME (e.g., 'phrog_1')
            accession = hmm_hit['gene_name']
            function, category = self.get_function_from_catalog(accession, ok_if_missing_from_catalog=True)

            for source, value in [(self.function_source, function), (self.category_source, category)]:
                functions_dict[counter] = {
                    'gene_callers_id': hmm_hit['gene_callers_id'],
                    'source': source,
                    'accession': accession,
                    'function': value,
                    'e_value': hmm_hit['e_value'],
                }
                counter += 1

        if functions_dict:
            gene_function_calls_table.create(functions_dict)
        else:
            self.run.warning(f"PHROGs class has no hits to process. Returning empty handed, but still adding "
                             f"{self.function_source} and {self.category_source} as functional sources.")
            gene_function_calls_table.add_empty_sources_to_functional_sources({self.function_source, self.category_source})

        if anvio.DEBUG:
            self.run.warning("The temp directories, '%s' and '%s' are kept. Please don't forget to clean those up later" %
                             (tmp_directory_path, ', '.join(hmmer.tmp_dirs)), header="Debug")
        else:
            self.run.info_single('Cleaning up the temp directory (you can use `--debug` if you would like to keep it for testing purposes)',
                                 nl_before=1, nl_after=1)
            shutil.rmtree(tmp_directory_path)
            hmmer.clean_tmp_dirs()

import os
import re
import glob
import json
import time
import shutil
import pandas as pd
import multiprocessing as mp

import anvio
import anvio.utils as utils
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths

from anvio.errors import ConfigError
from anvio.parsers import parser_modules
from anvio.drivers.hmmer import HMMer
from anvio.drivers.muscle import Muscle

from anvio.metabolism.setup import KeggSetup
from anvio.metabolism.constants import STRAY_KO_ANVIO_SUFFIX, GLOBAL_MAP_ID_PATTERN


run_quiet = terminal.Run(log_file_path=None, verbose=False)
progress_quiet = terminal.Progress(verbose=False)


class KOfamDownload(KeggSetup):
    """Class for setting up KOfam HMM profiles.

    Parameters
    ==========
    args: Namespace object
        All the arguments supplied by user to command-line programs relying on this
        class, such as `anvi-setup-kegg-data`. If using this class through the API, please
        provide a Namespace object with the Boolean 'reset' parameter.
    skip_init: Boolean
        Developers can use this flag to skip the sanity checks and creation of directories
        when testing this class.
    """

    def __init__(self, args, run=terminal.Run(), progress=terminal.Progress(), skip_init=False):
        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        self.args = args
        self.run = run
        self.progress = progress
        self.skip_init = skip_init
        self.include_stray_kos = True if A('include_stray_KOs') else False

        self.run.info_single("Info from KOfam Download")
        self.run.info("Stray KOs will be processed (`--include-stray-KOs` flag)", self.include_stray_kos)

        KeggSetup.__init__(self, self.args, skip_init=self.skip_init)

        filesnpaths.is_program_exists('hmmpress')
        if self.include_stray_kos:
            filesnpaths.is_program_exists('hmmbuild')
            filesnpaths.is_program_exists('muscle')

        # ftp path for HMM profiles and KO list
            # for ko list, add /ko_list.gz to end of url
            # for profiles, add /profiles.tar.gz  to end of url
        self.database_url = "ftp://ftp.genome.jp/pub/db/kofam"
        # dictionary mapping downloaded file name to final decompressed file name or folder location
        self.kofam_files = {'ko_list.gz': self.ko_list_file_path, 'profiles.tar.gz': self.kegg_data_dir}

        expected_files_for_kofams = [self.ko_list_file_path]
        if self.only_processing:
            expected_files_for_kofams.append(os.path.join(self.kegg_data_dir, 'profiles.tar.gz'))
        else:
            expected_files_for_kofams.append(self.kofam_hmm_file_path)

        if not args.reset and not anvio.DEBUG and not self.skip_init:
            self.is_database_exists(expected_files_for_kofams, fail_if_exists=(not self.only_processing))

        if self.download_from_kegg and not self.only_processing and not self.kegg_archive_path and not self.skip_init:
            filesnpaths.gen_output_directory(self.kegg_hmm_data_dir, delete_if_exists=args.reset)
            filesnpaths.gen_output_directory(self.orphan_data_dir, delete_if_exists=args.reset)


    def download_profiles(self):
        """This function downloads the Kofam profiles."""

        self.run.info("Kofam Profile Database URL", self.database_url)

        try:
            for file_name in self.kofam_files.keys():
                utils.download_file(self.database_url + '/' + file_name,
                    os.path.join(self.kegg_data_dir, file_name), progress=self.progress, run=self.run)
        except Exception as e:
            print(e)
            raise ConfigError("Anvi'o failed to download KEGG KOfam profiles from the KEGG website. Something "
                              "likely changed on the KEGG end. Please contact the developers to see if this is "
                              "a fixable issue. If it isn't, we may be able to provide you with a legacy KEGG "
                              "data archive that you can use to setup KEGG with the --kegg-archive flag.")


    def decompress_profiles(self):
        """This function decompresses the Kofam profiles."""

        self.progress.new('Decompressing files')
        for file_name in self.kofam_files.keys():
            self.progress.update('Decompressing file %s' % file_name)
            full_path = os.path.join(self.kegg_data_dir, file_name)

            if full_path.endswith("tar.gz"):
                utils.tar_extract_file(full_path, output_file_path=self.kofam_files[file_name], keep_original=False)
            else:
                utils.gzip_decompress_file(full_path, output_file_path=self.kofam_files[file_name], keep_original=False)

            self.progress.update("File decompressed. Yay.")
        self.progress.end()


    def confirm_downloaded_profiles(self):
        """This function verifies that all Kofam profiles have been properly downloaded.

        It is intended to be run after the files have been decompressed. The profiles directory should contain hmm files (ie, K00001.hmm);
        all KO numbers from the ko_list file (except those in ko_skip_list) should be included.

        This function must be called after setup_ko_dict() so that the self.ko_dict attribute is established.
        """

        ko_nums = self.ko_dict.keys()
        for k in ko_nums:
            if k not in self.ko_skip_list:
                hmm_path = os.path.join(self.kegg_data_dir, f"profiles/{k}.hmm")
                if not os.path.exists(hmm_path):
                    raise ConfigError(f"The KOfam HMM profile at {hmm_path} does not exist. This probably means that something went wrong "
                                      f"while downloading the KOfam database. Please run `anvi-setup-kegg-data` with the --reset "
                                      f"flag. If that still doesn't work, please contact the developers to see if the issue is fixable. "
                                      f"If it isn't, we may be able to provide you with a legacy KEGG data archive that you can use to "
                                      f"setup KEGG with the --kegg-archive flag.")


    def move_orphan_files(self):
        """This function moves the following to the orphan files directory:

            - profiles that do not have ko_list entries
            - profiles whose ko_list entries have no scoring threshold (in ko_no_threshold_list)

        And, the following profiles should not have been downloaded, but if they were then we move them, too:
            - profiles whose ko_list entries have no data at all (in ko_skip_list)
        """

        if not os.path.exists(self.orphan_data_dir): # should not happen but we check just in case
            raise ConfigError(f"Hmm. Something is out of order. The orphan data directory {self.orphan_data_dir} does not exist "
                              "yet, but it needs to in order for the move_orphan_files() function to work.")

        no_kofam_path = os.path.join(self.orphan_data_dir, "hmm_profiles_with_no_kofams.hmm")
        no_kofam_file_list = []
        no_threshold_file_list = []
        no_data_path = os.path.join(self.orphan_data_dir, "hmm_profiles_with_kofams_with_no_data.hmm")
        no_data_file_list = []

        hmm_list = [k for k in glob.glob(os.path.join(self.kegg_data_dir, 'profiles/*.hmm'))]
        for hmm_file in hmm_list:
            ko = re.search('profiles/(K\d{5})\.hmm', hmm_file).group(1)
            if ko not in self.ko_dict.keys():
                if ko in self.ko_no_threshold_list:
                    no_threshold_file_list.append(hmm_file)
                elif ko in self.ko_skip_list: # these should not have been downloaded, but if they were we will move them
                    no_data_file_list.append(hmm_file)
                else:
                    no_kofam_file_list.append(hmm_file)

        # now we concatenate the orphan and stray KO hmms into the orphan data directory
        if no_kofam_file_list:
            utils.concatenate_files(no_kofam_path, no_kofam_file_list, remove_concatenated_files=True)
            self.progress.reset()
            self.run.warning(f"Please note that while anvi'o was building your databases, she found {len(no_kofam_file_list)} "
                             f"HMM profiles that did not have any matching KOfam entries. We have removed those HMM "
                             f"profiles from the final database. You can find them under the directory '{self.orphan_data_dir}'.")

        if no_threshold_file_list:
            utils.concatenate_files(self.stray_ko_hmms_from_kegg, no_threshold_file_list, remove_concatenated_files=False)
            filesnpaths.gen_output_directory(os.path.join(self.orphan_data_dir, "profiles"), delete_if_exists=True)
            for k_path in no_threshold_file_list:
                k = os.path.basename(k_path)
                # move individual profiles temporarily to the orphan data dir, so they don't get combined with the regular KOs
                # but we can still use them later if necessary for --include-stray-KOs
                os.rename(k_path, os.path.join(self.orphan_data_dir, f"profiles/{k}"))
            self.progress.reset()
            self.run.warning(f"Please note that while anvi'o was building your databases, she found {len(no_threshold_file_list)} "
                             f"KOfam entries that did not have any threshold to remove weak hits. We have removed those HMM "
                             f"profiles from the final database. You can find them under the directory '{self.orphan_data_dir}'. "
                             f"If you used the flag --include-stray-KOs, we will estimate their bit score thresholds using KEGG GENES "
                             f"data so that you can annotate these KOs downstream if you wish.")

        if no_data_file_list:
            utils.concatenate_files(no_data_path, no_data_file_list, remove_concatenated_files=True)
            self.progress.reset()
            self.run.warning(f"Please note that while anvi'o was building your databases, she found {len(no_data_file_list)} "
                             f"HMM profiles that did not have any associated data (besides an annotation) in their KOfam entries. "
                             f"We have removed those HMM profiles from the final database. You can find them under the directory "
                             f"'{self.orphan_data_dir}'.")


    def exec_hmmpress_command_on_ko_file(self, hmm_file_path, log_file_path):
        """Given a path to a set of KO HMMs and a log file path, this function executes the appropriate
        `hmmpress` command and deletes the log file afterwards if it was successful.
        """

        cmd_line = ['hmmpress', hmm_file_path]
        ret_val = utils.run_command(cmd_line, log_file_path)

        if ret_val:
            raise ConfigError("Hmm. There was an error while running `hmmpress` on the Kofam HMM profiles. "
                              "Check out the log file ('%s') to see what went wrong." % (log_file_path))
        else:
            # getting rid of the log file because hmmpress was successful
            os.remove(log_file_path)


    def run_hmmpress(self):
        """This function concatenates the Kofam profiles and runs hmmpress on them."""

        self.progress.new('Preparing Kofam HMM Profiles')

        self.progress.update('Verifying the Kofam directory %s contains all HMM profiles' % self.kegg_data_dir)
        self.confirm_downloaded_profiles()

        self.progress.update('Handling orphan files')
        self.move_orphan_files()

        self.progress.update('Concatenating HMM profiles into one file...')
        hmm_list = [k for k in glob.glob(os.path.join(self.kegg_data_dir, 'profiles/*.hmm'))]
        utils.concatenate_files(self.kofam_hmm_file_path, hmm_list, remove_concatenated_files=False)

        self.progress.update('Running hmmpress on KOs...')
        self.exec_hmmpress_command_on_ko_file(self.kofam_hmm_file_path, os.path.join(self.kegg_hmm_data_dir, '00_hmmpress_log.txt'))

        self.progress.end()


    def download_ko_files(self, kos_to_download, destination_dir, dont_raise=True):
        """Multi-threaded download of KEGG Orthology files.

        Parameters
        ==========
        kos_to_download: list of str
            List of KOs to download Orthology files for
        destination_dir: file path
            Where to download the files to
        dont_raise : Boolean
            If True (default), this function won't raise an error if some files failed to download.
        Returns
        =======
        undownloaded : list of str
            List of KOs that failed to download (will be empty if all were successful).
        """

        num_kos = len(kos_to_download)
        self.progress.new('Downloading KO files', progress_total_items=num_kos)

        manager = mp.Manager()
        input_queue = manager.Queue()
        output_queue = manager.Queue()
        for ko in kos_to_download:
            ko_file_path = os.path.join(destination_dir, ko)
            url = self.kegg_rest_api_get + '/' + ko
            input_queue.put((url, ko_file_path))
        workers: List[mp.Process] = []
        for _ in range(self.num_threads):
            worker = mp.Process(target=_download_worker, args=(input_queue, output_queue))
            workers.append(worker)
            worker.start()

        downloaded_count = 0
        undownloaded_count = 0
        undownloaded = []
        while downloaded_count + undownloaded_count < num_kos:
            output = output_queue.get()
            if output is True:
                downloaded_count += 1
                self.progress.update(f"{downloaded_count} / {num_kos} KO files downloaded")
                self.progress.increment(increment_to=downloaded_count)
            else:
                undownloaded_count += 1
                undownloaded.append(os.path.splitext(os.path.basename(output))[0])

        self.progress.end()
        for worker in workers:
            worker.terminate()
        if undownloaded:
            if dont_raise:
                self.run.warning(f"Files for the following KOs failed to download despite multiple attempts: "
                                f"{', '.join(undownloaded)}. If this is unacceptable to you, you can try to "
                                f"re-run this program to see if things will work on the next try.")
            else:
                raise ConfigError(f"Files for the following KOs failed to download despite multiple attempts: "
                                  f"{', '.join(undownloaded)}. Since the function responsible for handling this was "
                                  f"told to quit should this happen, well, here we are. If skipping these failed KOs "
                                  f"is okay, you could always run this function with `dont_raise=True`.")

        return undownloaded


    def download_kegg_genes_files(self, genes_to_download, destination_dir, dont_raise=True):
        """Multi-threaded download of KEGG GENES files.

        Parameters
        ==========
        genes_to_download: list of str
            List of KEGG GENES accessions to download. Example format: "ctc:CTC_p60" (lowercase organism code, colon, gene accession)
        destination_dir: file path
            Where to download the files to
        dont_raise : Boolean
            If True (default), this function won't raise an error if some files failed to download.
        Returns
        =======
        undownloaded : List
            List of files that failed to download (will be empty if all were successful).
        """

        num_genes = len(genes_to_download)
        self.progress.new('Downloading KEGG GENES files', progress_total_items=num_genes)

        manager = mp.Manager()
        input_queue = manager.Queue()
        output_queue = manager.Queue()
        for g in genes_to_download:
            genes_file_path = os.path.join(destination_dir, g)
            url = self.kegg_rest_api_get + '/' + g
            input_queue.put((url, genes_file_path))
        workers: List[mp.Process] = []
        for _ in range(self.num_threads):
            worker = mp.Process(target=_download_worker, args=(input_queue, output_queue))
            workers.append(worker)
            worker.start()

        downloaded_count = 0
        undownloaded_count = 0
        undownloaded = []
        while downloaded_count + undownloaded_count < num_genes:
            output = output_queue.get()
            if output is True:
                downloaded_count += 1
                self.progress.update(f"{downloaded_count} / {num_genes} KEGG GENES files downloaded")
                self.progress.increment(increment_to=downloaded_count)
            else:
                undownloaded_count += 1
                undownloaded.append(os.path.splitext(os.path.basename(output))[0])

        self.progress.end()
        for worker in workers:
            worker.terminate()
        if undownloaded:
            if dont_raise:
                self.run.warning(f"Files for the following KEGG GENES failed to download despite multiple attempts: "
                                f"{', '.join(undownloaded)}. If this is unacceptable to you, you can try to "
                                f"re-run this program to see if things will work on the next try.")
            else:
                raise ConfigError(f"Files for the following KEGG GENES failed to download despite multiple attempts: "
                                  f"{', '.join(undownloaded)}. Since the function responsible for handling this was "
                                  f"told to quit should this happen, well, here we are. If skipping these failed KOs "
                                  f"is okay, you could always run this function with `dont_raise=True`.")

        return undownloaded


    def get_kegg_gene_accessions_from_ko_files(self, ko_list, ko_file_dir):
        """Extracts KEGG GENES accessions from KO files and returns a dictionary mapping KO to its GENES.

        Parameters
        ==========
        ko_list: list of str
            List of KEGG accessions to process.
        ko_file_dir: file path
            Where the KO files are located
        Returns
        =======
        ko_to_genes : dict
            Dictionary with KOs as keys and list of KEGG GENES accessions as values
        """

        ko_to_genes = {k : [] for k in ko_list}

        for ko in ko_list:
            ko_file_path = os.path.join(ko_file_dir, ko)
            genes_acc_list = self.extract_data_field_from_kegg_file(ko_file_path, target_field="GENES")

            kegg_genes_code_list = []
            for i, acc in enumerate(genes_acc_list):
                acc_fields = acc.split(": ")            # example accession is "CTC: CTC_p60(tetX)"
                org_code = acc_fields[0].lower()        # the organism code (before the colon) needs to be converted to lowercase
                gene_name = acc_fields[1]

                # sometimes we have multiple genes per organism, like this: "PSOM: 113322169 113340172"
                all_genes = gene_name.split(' ')
                for g in all_genes:
                    g = g.split('(')[0] # the gene name needs to have anything in parentheses removed. ex. CTC_p60(tetX) becomes CTC_p60
                    kegg_genes_code = f"{org_code}:{g}"
                    kegg_genes_code_list.append(kegg_genes_code)

            ko_to_genes[ko] = kegg_genes_code_list

        return ko_to_genes


    def kegg_gene_sequences_to_fasta_file(self, kegg_genes_files, target_fasta_file):
        """This function extracts the amino acid sequences for a list of KEGG GENES and prints them to a FASTA file.

        Parameters
        ==========
        kegg_genes_files : List of str
            List of paths to KEGG GENES file to extract sequences from
        target_fasta_file : list of str
            Path to FASTA file in which to store the sequences

        Returns
        =======
        seq_tuples : List of tuples
            Each sequence added to the FASTA file is also returned in this list, where each tuple contains
            (KEGG GENES name, amino acid sequence). Note that the seq name is taken from the name of the KEGG GENES file.
        """

        seq_tuples = []
        for i, gene_file_path in enumerate(kegg_genes_files):
            seq_name = os.path.basename(gene_file_path)
            # obtain the amino acid sequence and save it to the fasta file
            aa_sequence_data = self.extract_data_field_from_kegg_file(gene_file_path, target_field="AASEQ")

            aaseq = ""
            with open(target_fasta_file, 'a') as fasta:
                fasta.write(f">{i}\n") # we label the gene with its index because the HMMER parser expects an int, not a string, as the gene name
                for seq in aa_sequence_data[1:]: # we skip the first element, which is the sequence length
                    fasta.write(f"{seq}\n")
                    aaseq += seq
            seq_tuples.append((seq_name, aaseq))

        return seq_tuples


    def build_HMM_from_seqs(self, hmm_name, tuple_of_seqs, hmm_output_file, log_file_path):
        """This function aligns sequences and builds an HMM from them using `muscle` and `hmmbuild`.

        Parameters
        ==========
        hmm_name : str
            What to name the model (ie 'NAME' field in the .hmm file)
        tuple_of_seqs : List of (sequence name, sequence) tuples
            The sequences to align with 'muscle' to create the `hmmbuild` input.
            See anvio.drivers.muscle for example format
        hmm_output_file : str
            File path where to store the new HMM model
        log_file_path : str
            File path for the log file of `hmmbuild`
        """

        if len(tuple_of_seqs) < 2:
            raise ConfigError(f"The function build_HMM_from_seqs() can't build an alignment from less than "
                              f"2 sequences, but that is what it got. Here is the sequence (if any) passed to this "
                              f"function: {tuple_of_seqs}. No alignment, no HMM. Sorry!")

        m = Muscle(progress=progress_quiet, run=run_quiet)
        clw_alignment = m.run_stdin(tuple_of_seqs, debug=anvio.DEBUG, clustalw_format=True)

        hmmbuild_cmd_line = ['hmmbuild', '-n', hmm_name, '--informat', 'clustallike', hmm_output_file, '-'] # sending '-' in place of an alignment file so it reads from stdin
        utils.run_command_STDIN(hmmbuild_cmd_line, log_file_path, clw_alignment)

        if not os.path.exists(hmm_output_file):
            raise ConfigError(f"It seems that the `hmmbuild` command failed because there is no output model at {hmm_output_file}. "
                              f"Perhaps the log file {log_file_path} will hold some answers for you.")


    def estimate_bitscore_for_ko(self, ko, kegg_genes_for_ko, kegg_genes_fasta, ko_model_file):
        """This function estimates the bitscore of a single KEGG Ortholog.

        It runs `hmmscan` of the KO model against the provided list of its KEGG GENE
        sequences, and then computes the minimum bit score to use as a threshold for
        annotating this protein family.

        Parameters
        ==========
        ko : str
            KEGG identifier for the KO
        kegg_genes_for_ko : list of str
            List of KEGG GENE accessions that were used to generate the KO model (for sanity check and
            number of sequences)
        kegg_genes_fasta : str
            Path to FASTA file where the KEGG GENES sequences for this KO are stored
        ko_model_file : str
            File path of the .hmm file containg the KO model (doesn't need to contain only this model,
            but must be hmmpressed already)
        Returns
        =======
        threshold : float
            estimated bit score threshold for the KO's HMM. Will be None if kegg_genes_for_ko is empty
        """

        # sanity check for empty KEGG GENES list
        if not kegg_genes_for_ko:
            self.run.warning(f"The function estimate_bitscore_for_ko() received an empty list of KEGG GENES "
                             f"for {ko}, so it cannot estimate a bit score threshold. The function will return "
                             f"a threshold of `None` for this KO.")
            return None

        # we run hmmscan of the KO against its associated GENES sequences and process the hits
        target_file_dict = {'AA:GENE': kegg_genes_fasta}
        hmmer = HMMer(target_file_dict, num_threads_to_use=self.num_threads, progress=progress_quiet, run=run_quiet)
        hmm_hits_file = hmmer.run_hmmer('KO {ko}', 'AA', 'GENE', None, None, len(kegg_genes_for_ko), ko_model_file, None, None)

        if not hmm_hits_file:
            raise ConfigError(f"No HMM hits were found for the KO model {ko}. This is seriously concerning, because we were running it against "
                              f"gene sequences that were used to generate the model.")

        parser = parser_modules['search']['hmmer_table_output'](hmm_hits_file, alphabet='AA', context='GENE', run=run_quiet)
        search_results_dict = parser.get_search_results()

        # take the minimum of hits from current KO model as bit score threshold
        all_relevant_bitscores = []
        for hit, hit_info_dict in search_results_dict.items():
            if hit_info_dict['gene_name'] == ko or hit_info_dict['gene_name'] == f"{ko}{STRAY_KO_ANVIO_SUFFIX}":
                all_relevant_bitscores.append(hit_info_dict['bit_score'])

        threshold = min(all_relevant_bitscores)
        return threshold


    def process_all_stray_kos(self):
        """This driver function processes each stray KO and creates a file of bit score thresholds for them.

        The following steps are run for each stray KO:
        1. download of its KO file
        2. identification and download of the KEGG GENES sequences in this family
        3. alignment with `muscle` and `hmmbuild` to create a new model (since KEGG GENES updates faster than KOfam models do)
        4. `hmmscan` of the new KO model against these sequences to get bit scores
        5. computing the minimum bit score to use as a threshold for annotating this family
        """

        num_strays = len(self.ko_no_threshold_list)
        self.run.info("Number of Stray KOs to process", num_strays)

        self.stray_ko_file_dir = os.path.join(self.orphan_data_dir, "00_STRAY_KO_FILES")
        self.stray_ko_genes_dir = os.path.join(self.orphan_data_dir, "01_STRAY_GENES_FILES")
        self.stray_ko_seqs_dir = os.path.join(self.orphan_data_dir, "02_STRAY_GENES_FASTA")
        self.stray_ko_hmms_dir = os.path.join(self.orphan_data_dir, "03_STRAY_KO_HMMS")
        filesnpaths.gen_output_directory(self.stray_ko_file_dir, delete_if_exists=True)
        filesnpaths.gen_output_directory(self.stray_ko_genes_dir, delete_if_exists=True)
        filesnpaths.gen_output_directory(self.stray_ko_seqs_dir, delete_if_exists=True)
        filesnpaths.gen_output_directory(self.stray_ko_hmms_dir, delete_if_exists=True)

        ko_files_not_downloaded = self.download_ko_files(self.ko_no_threshold_list, self.stray_ko_file_dir)

        ko_files_to_process = list(set(self.ko_no_threshold_list) - set(ko_files_not_downloaded))
        self.run.info("Number of Stray KO files successfully downloaded", len(ko_files_to_process))
        if ko_files_not_downloaded:
            self.run.warning(f"FYI, some stray KOs failed to download from KEGG. We will not estimate their bit score thresholds "
                             f"and you will not be able to annotate them later. Here they are: {', '.join(ko_files_not_downloaded)}")

        ko_to_gene_accessions = self.get_kegg_gene_accessions_from_ko_files(ko_files_to_process, self.stray_ko_file_dir)
        kegg_genes_to_download = []
        for acc_list in ko_to_gene_accessions.values():
            kegg_genes_to_download.extend(acc_list)

        kegg_genes_not_downloaded = self.download_kegg_genes_files(kegg_genes_to_download, self.stray_ko_genes_dir)
        kegg_genes_downloaded = list(set(kegg_genes_to_download) - set(kegg_genes_not_downloaded))
        self.run.info("Number of KEGG GENES files successfully downloaded", len(kegg_genes_downloaded))
        if kegg_genes_not_downloaded:
            self.run.warning(f"FYI, some KEGG GENES files failed to download from KEGG. They will not be used for "
                             f"estimating bit score thresholds. Here they are: {', '.join(kegg_genes_not_downloaded)}")

        self.progress.new("Extracting amino acid sequences for Stray KOs", progress_total_items=len(ko_files_to_process))
        ko_to_gene_seqs_list = {} # we'll store the sequences to align with muscle here. yes, we just stored them in a file.
        cur_num = 0
        for k in ko_files_to_process:
            self.progress.update(f"Working on {k} [{cur_num} of {len(ko_files_to_process)}]")
            self.progress.increment(increment_to=cur_num)
            downloaded_genes_list = [a for a in ko_to_gene_accessions[k] if a in kegg_genes_downloaded]
            gene_file_paths = [os.path.join(self.stray_ko_genes_dir, code) for code in downloaded_genes_list]
            ko_to_gene_seqs_list[k] = self.kegg_gene_sequences_to_fasta_file(gene_file_paths, os.path.join(self.stray_ko_seqs_dir, f"GENES_FOR_{k}.fa"))
            cur_num += 1
        self.progress.end()

        self.progress.new("Aligning genes and creating new HMMs for Stray KOs", progress_total_items=len(ko_files_to_process))
        list_of_new_HMMs = []
        hmmbuild_log = os.path.join(self.orphan_data_dir, "hmmbuild.log")
        cur_num = 0
        new_models = 0
        old_models = 0
        models_without_genes = []
        kos_with_one_gene = []
        models_with_anvio_version = []
        for k in ko_files_to_process:
            self.progress.update(f"Working on {k} [{cur_num} of {len(ko_files_to_process)}]")
            self.progress.increment(increment_to=cur_num)
            if not len(ko_to_gene_seqs_list[k]):
                models_without_genes.append(k)
            elif len(ko_to_gene_seqs_list[k]) == 1:
                # with only 1 sequence, we can't build a new model. We can try to use KEGG's since it is guaranteed to fit this sequence
                kegg_model_file = os.path.join(self.orphan_data_dir, f"profiles/{k}.hmm")
                if not os.path.exists(kegg_model_file):
                    kos_with_one_gene.append(k) # if we don't have the OG model, we just skip this one
                else:
                    list_of_new_HMMs.append(kegg_model_file)
                    old_models += 1
            else:
                hmm_model_file = os.path.join(self.stray_ko_hmms_dir, f"{k}_anvio.hmm")
                self.build_HMM_from_seqs(f"{k}{STRAY_KO_ANVIO_SUFFIX}", ko_to_gene_seqs_list[k], hmm_model_file, hmmbuild_log)
                list_of_new_HMMs.append(hmm_model_file)
                new_models += 1
                models_with_anvio_version.append(k)
            cur_num += 1
        self.progress.end()

        self.run.info("Number of Stray KOs with new HMMs built by anvi'o to incorporate potentially new KEGG GENES", new_models)
        self.run.info("Number of Stray KOs using KEGG's original HMM because the family includes only one gene sequence", old_models)
        if models_without_genes:
            self.run.warning(f"We weren't able to download any KEGG GENE sequences for some stray KOs, and therefore will not "
                             f"be able to estimate bit score threshold for these KOs. Here they are: {', '.join(models_without_genes)}")
            self.run.info("Number of KOs without downloaded gene sequences", len(models_without_genes))
            ko_files_to_process = list(set(ko_files_to_process) - set(models_without_genes))
        if kos_with_one_gene:
            self.run.warning(f"The following stray KOs had exactly one KEGG GENE sequence, so we couldn't build a new HMM for them, "
                             f"but we also couldn't find their models from KEGG, so we won't estimate bit score thresholds for them: "
                             f"{', '.join(kos_with_one_gene)}")
            self.run.info("Number of KOs without HMMs", len(kos_with_one_gene))
            ko_files_to_process = list(set(ko_files_to_process) - set(kos_with_one_gene))

        self.progress.new('Concatenating new Stray HMM files...')
        utils.concatenate_files(self.stray_ko_hmm_file_path, list_of_new_HMMs, remove_concatenated_files=False)
        self.progress.update('Running hmmpress on new Stray KO HMMs...')
        self.exec_hmmpress_command_on_ko_file(self.stray_ko_hmm_file_path, os.path.join(self.orphan_data_dir, '00_hmmpress_log.txt'))
        self.progress.end()
        self.run.info("File storing all new HMMs generated for Stray KOs", self.stray_ko_hmm_file_path)

        self.progress.new("Estimating bit score threshold for Stray KOs", progress_total_items=len(ko_files_to_process))
        threshold_dict = {}
        cur_num = 0
        for k in ko_files_to_process:
            self.progress.update(f"Working on {k} [{cur_num} of {len(ko_files_to_process)}]")
            self.progress.increment(increment_to=cur_num)
            downloaded_genes_list = [a for a in ko_to_gene_accessions[k] if a in kegg_genes_downloaded]
            threshold_dict[k] = self.estimate_bitscore_for_ko(k, kegg_genes_for_ko=downloaded_genes_list,
                                        kegg_genes_fasta=os.path.join(self.stray_ko_seqs_dir, f"GENES_FOR_{k}.fa"),
                                        ko_model_file=self.stray_ko_hmm_file_path)
            cur_num += 1
        self.progress.end()

        # we need to re-load the ko dictionary so that we have access to the definitions of the stray KOs
        # cannot do this before this point because the absence of an stray KO from this dict controls whether it is moved to the
        # stray data directory (and we want to keep the strays separate since we process them specially)
        self.setup_ko_dict(exclude_threshold=(not self.include_stray_kos), suppress_warnings=True)

        # write the thresholds to a file
        thresholds_not_none = 0
        with open(self.stray_ko_thresholds_file, 'w') as out:
            out.write("knum\tthreshold\tscore_type\tdefinition\n")
            for k, thr in threshold_dict.items():
                if thr:
                    model_name = k
                    if k in models_with_anvio_version:
                        model_name = f"{k}{STRAY_KO_ANVIO_SUFFIX}"
                    ko_definition = self.ko_dict[k]['definition']
                    out.write(f"{model_name}\t{thr}\tfull\t{ko_definition}\n")
                    thresholds_not_none += 1
        self.run.info("File with estimated bit score thresholds", self.stray_ko_thresholds_file)
        self.run.info("Number of estimated thresholds", thresholds_not_none)

        # clean up downloaded files
        if not anvio.DEBUG:
            os.remove(hmmbuild_log)
            for d in [self.stray_ko_file_dir, self.stray_ko_genes_dir, self.stray_ko_seqs_dir, self.stray_ko_hmms_dir]:
                shutil.rmtree(d)
            self.run.warning("The KO and GENES files downloaded from KEGG for processing the stray KOs, as well as the "
                             "individual new HMM files that anvi'o generated for these KOs, are now deleted to save space. "
                             "If you want to keep them, next time run the program with `--debug`.")


    def setup_kofams(self):
        """This function downloads, decompresses, and runs `hmmpress` on KOfam profiles."""

        if not self.only_processing:
            self.download_profiles()

        if not self.only_download:
            self.decompress_profiles()
            self.setup_ko_dict() # get ko dict attribute
            self.run_hmmpress()

            if self.include_stray_kos:
                self.process_all_stray_kos()

            # there is no reason to keep the original HMM profiles around, unless we are debugging
            if not anvio.DEBUG:
                shutil.rmtree(os.path.join(self.kegg_data_dir, "profiles"))
                shutil.rmtree(os.path.join(self.orphan_data_dir, "profiles"))


class ModulesDownload(KeggSetup):
    """Class for setting up all KEGG data related to pathway prediction, namely KOfam profiles and KEGG MODULES;
    reaction networks, which require MODULES, BRITE, and binary relation files;
    and pathway map images and reference KO, EC, and RN KGML files.

    Parameters
    ==========
    args: Namespace object
        All the arguments supplied by user to command-line programs relying on this
        class, such as `anvi-setup-kegg-data`. If using this class through the API, please
        provide a Namespace object with the Boolean 'reset' parameter.
    skip_init: Boolean
        Developers can use this flag to skip the sanity checks and creation of directories
        when testing this class.
    """

    def __init__(self, args, run=terminal.Run(), progress=terminal.Progress(), skip_init=False):
        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        self.args = args
        self.run = run
        self.progress = progress
        self.skip_init = skip_init
        self.overwrite_modules_db = A('overwrite_output_destinations')

        self.run.info_single("Info from MODULES Download")

        # we also need the init of the superclass
        KeggSetup.__init__(self, self.args, skip_init=self.skip_init)

        # download from KEGG option: module/pathway map htext files and API link
        self.kegg_module_download_path = "https://www.genome.jp/kegg-bin/download_htext?htext=ko00002.keg&format=htext&filedir="
        self.kegg_pathway_download_path = "https://www.genome.jp/kegg-bin/download_htext?htext=br08901.keg&format=htext&filedir="
        self.kegg_rest_api_get = "http://rest.kegg.jp/get"

        # check if the data is already downloaded
        expected_files_for_modules = [self.kegg_module_file,
                                      self.kegg_module_data_dir]

        if not args.reset and not anvio.DEBUG and not self.skip_init:
            self.is_database_exists(expected_files_for_modules, fail_if_exists=(not self.only_processing))

        # generate subfolders if necessary
        if self.download_from_kegg and not self.only_processing and not self.kegg_archive_path and not self.skip_init:
            filesnpaths.gen_output_directory(self.kegg_module_data_dir, delete_if_exists=args.reset)

    def download_kegg_module_file(self):
        """This function downloads the KEGG module file, which tells us which module files to download."""

        # download the kegg module file, which lists all modules
        try:
            utils.download_file(self.kegg_module_download_path, self.kegg_module_file, progress=self.progress, run=self.run)
        except Exception as e:
            print(e)
            raise ConfigError("Anvi'o failed to download the KEGG Module htext file from the KEGG website. Something "
                              "likely changed on the KEGG end. Please contact the developers to see if this is "
                              "a fixable issue. If it isn't, we may be able to provide you with a legacy KEGG "
                              "data archive that you can use to setup KEGG with the --kegg-archive flag.")


    def process_module_file(self):
        """This function reads the kegg module file into a dictionary. It should be called during setup to get the KEGG module numbers so that KEGG modules can be downloaded.

        The structure of this file is like this:

        +D    Module
        #<h2><a href="/kegg/kegg2.html"><img src="/Fig/bget/kegg3.gif" align="middle" border=0></a>&nbsp; KEGG Modules</h2>
        !
        A<b>Pathway modules</b>
        B
        B  <b>Carbohydrate metabolism</b>
        C    Central carbohydrate metabolism
        D      M00001  Glycolysis (Embden-Meyerhof pathway), glucose => pyruvate [PATH:map00010 map01200 map01100]
        D      M00002  Glycolysis, core module involving three-carbon compounds [PATH:map00010 map01200 map01230 map01100]
        D      M00003  Gluconeogenesis, oxaloacetate => fructose-6P [PATH:map00010 map00020 map01100]

        In other words, a bunch of initial lines to be ignored, and thereafter the line's information can be determined by the one-letter code at the start.
        A = Pathway modules (metabolic pathways) or signature modules (gene sets that indicate a phenotypic trait, ie toxins).
        B = Category of module (a type of metabolism for pathway modules. For signature modules, either Gene Set or Module Set)
        C = Sub-category of module
        D = Module

        """
        self.module_dict = {}

        filesnpaths.is_file_exists(self.kegg_module_file)
        filesnpaths.is_file_plain_text(self.kegg_module_file)

        f = open(self.kegg_module_file, 'r')
        self.progress.new("Parsing KEGG Module file")

        current_module_type = None
        current_category = None
        current_subcategory = None

        for line in f.readlines():
            line = line.strip('\n')
            first_char = line[0]

            # garbage lines
            if first_char in ["+", "#", "!"]:
                continue
            else:
                # module type
                if first_char == "A":
                    fields = re.split('<[^>]*>', line) # we split by the html tag here
                    current_module_type = fields[1]
                # Category
                elif first_char == "B":
                    fields = re.split('<[^>]*>', line) # we split by the html tag here
                    if len(fields) == 1: # sometimes this level has lines with only a B
                        continue
                    current_category = fields[1]
                # Sub-category
                elif first_char == "C":
                    fields = re.split('\s{2,}', line) # don't want to split the subcategory name, so we have to split at least 2 spaces
                    current_subcategory = fields[1]
                # module
                elif first_char == "D":
                    fields = re.split('\s{2,}', line)
                    mnum = fields[1]
                    self.module_dict[mnum] = {"name" : fields[2], "type" : current_module_type, "category" : current_category, "subcategory" : current_subcategory}
                # unknown code
                else:
                    raise ConfigError("While parsing the KEGG file %s, we found an unknown line code %s. This has "
                                      "made the file unparseable. It is likely that an update to KEGG has broken "
                                      "things such that anvi'o doesn't know what is going on anymore. Sad, we know. :( "
                                      "Please contact the developers to see if this is a fixable issue, and in the "
                                      "meantime use an older version of the KEGG data directory (if you have one). "
                                      "If we cannot fix it, we may be able to provide you with a legacy KEGG "
                                      "data archive that you can use to setup KEGG with the --kegg-archive flag." % (self.kegg_module_file, first_char))
        self.progress.end()


    def download_modules(self):
        """This function downloads the KEGG modules."""

        total = len(self.module_dict.keys())
        self.run.info("KEGG Module Database URL", self.kegg_rest_api_get)
        self.run.info("Number of KEGG Modules to download", total)
        self.run.info("Number of threads used for download", self.num_threads)

        self.progress.new("Downloading KEGG Module files", progress_total_items=total)
        manager = mp.Manager()
        input_queue = manager.Queue()
        output_queue = manager.Queue()
        for mnum in self.module_dict.keys():
            file_path = os.path.join(self.kegg_module_data_dir, mnum)
            url = self.kegg_rest_api_get + '/' + mnum
            input_queue.put((url, file_path))
        workers: List[mp.Process] = []
        for _ in range(self.num_threads):
            worker = mp.Process(target=_download_worker, args=(input_queue, output_queue))
            workers.append(worker)
            worker.start()

        downloaded_count = 0
        undownloaded_count = 0
        undownloaded = []
        while downloaded_count + undownloaded_count < total:
            output = output_queue.get()
            if output is True:
                downloaded_count += 1
                self.progress.update(f"{downloaded_count} / {total} module files downloaded")
                self.progress.increment(increment_to=downloaded_count)
            else:
                undownloaded_count += 1
                undownloaded.append(os.path.splitext(os.path.basename(output))[0])

        for worker in workers:
            worker.terminate()
        self.progress.end()

        if undownloaded:
            raise ConfigError(
                "Unfortunately, files for the following modules failed to download despite multiple attempts, "
                f"and so the database needs to be set up again: {', '.join(undownloaded)}"
            )


    def confirm_downloaded_modules(self):
        """This function verifies that all module files have been downloaded.

        It checks that there is a module file for every module in the self.module_dict dictionary;
        for that reason, it must be called after the function that creates that attribute,
        process_module_file(), has already been called. To verify that each file has been downloaded
        properly, we check that the last line is '///'.
        """

        for mnum in self.module_dict.keys():
            file_path = os.path.join(self.kegg_module_data_dir, mnum)
            if not os.path.exists(file_path):
                raise ConfigError(f"The module file for {mnum} does not exist at its expected location, {file_path}. "
                                  f"This probably means that something is wrong with your downloaded data, since this "
                                  f"module is present in the KEGG MODULE file that lists all modules you *should* have "
                                  f"on your computer. Very sorry to tell you this, but you need to re-download the KEGG "
                                  f"data. We recommend the --reset flag.")
            # verify entire file has been downloaded
            f = open(file_path, 'r')
            f.seek(0, os.SEEK_END)
            f.seek(f.tell() - 4, os.SEEK_SET)
            last_line = f.readline().strip('\n')
            if not last_line == '///':
                raise ConfigError("The KEGG module file %s was not downloaded properly. We were expecting the last line in the file "
                                  "to be '///', but instead it was %s. Formatting of these files may have changed on the KEGG website. "
                                  "Please contact the developers to see if this is a fixable issue. If it isn't, we may be able to "
                                  "provide you with a legacy KEGG data archive that you can use to setup KEGG with the --kegg-archive flag."
                                  % (file_path, last_line))
        self.run.info("Number of module files found", len(self.module_dict))


    def setup_modules_data(self):
        """This is a driver function which executes the setup process for pathway prediction and reaction network data from KEGG."""

        # FIXME: we will have to move user setup to a completely separate program at some point
        # PS. user setup related functions belong to the superclass for now
        if self.user_input_dir:
            self.setup_user_data()
        else:
            # download the data first
            # unless user requested only processing (mostly for developers and the adventurous)
            if not self.only_processing:
                self.download_kegg_module_file()
                self.process_module_file() # get module dict attribute
                self.download_modules()
                self.confirm_downloaded_modules()

            else:
                # get required attributes for database setup and make sure all expected files were downloaded
                self.process_module_file()
                self.confirm_downloaded_modules()

            # process the modules file into a database
            if not self.only_download:
                self.setup_modules_db(db_path=self.kegg_modules_db_path, module_data_directory=self.kegg_module_data_dir)

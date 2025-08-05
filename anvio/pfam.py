#!/usr/bin/env python
# -*- coding: utf-8
"""This file contains PfamSetup and Pfam classes."""

import anvio
import anvio.dbops as dbops
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths

from anvio.drivers.hmmer import HMMer
from anvio.parsers import parser_modules
from anvio.errors import ConfigError
from anvio.tables.genefunctions import TableForGeneFunctions

import os
import gzip
import numpy as np
import pandas as pd
import shutil
import requests
import glob

from io import BytesIO
from scipy.stats import entropy
from anvio.dbinfo import is_contigs_db
from anvio.utils.commandline import run_command
from anvio.utils.files import (
    get_TAB_delimited_file_as_dictionary,
    get_chunk,
    get_file_md5,
    gzip_decompress_file
)
from anvio.utils.network import download_file


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "Özcan Esen"
__email__ = "ozcanesen@gmail.com"


run = terminal.Run()
progress = terminal.Progress()
pp = terminal.pretty_print


def read_remote_file(url, is_gzip=True):
    remote_file = requests.get(url)

    if remote_file.status_code == 404:
        raise Exception("'%s' returned 404 Not Found. " % url)

    if is_gzip:
        buf = BytesIO(remote_file.content)
        fg = gzip.GzipFile(fileobj=buf)
        return fg.read().decode('utf-8')

    return remote_file.content.decode('utf-8')


class PfamSetup(object):
    def __init__(self, args, run=run, progress=progress):
        """Setup a Pfam database for anvi'o

        Parameters
        ==========
        args : argparse.Namespace
            See `bin/anvi-setup-pfams` for available arguments
        """

        self.args = args
        self.run = run
        self.progress = progress
        self.pfam_data_dir = args.pfam_data_dir

        filesnpaths.is_program_exists('hmmpress')

        if self.pfam_data_dir and args.reset:
            raise ConfigError("You are attempting to run Pfam setup on a non-default data directory (%s) using the --reset flag. "
                              "To avoid automatically deleting a directory that may be important to you, anvi'o refuses to reset "
                              "directories that have been specified with --pfam-data-dir. If you really want to get rid of this "
                              "directory and regenerate it with Pfam data inside, then please remove the directory yourself using "
                              "a command like `rm -r %s`. We are sorry to make you go through this extra trouble, but it really is "
                              "the safest way to handle things." % (self.pfam_data_dir, self.pfam_data_dir))

        if not self.pfam_data_dir:
            self.pfam_data_dir = os.path.join(os.path.dirname(anvio.__file__), 'data/misc/Pfam')

        filesnpaths.is_output_dir_writable(os.path.dirname(os.path.abspath(self.pfam_data_dir)))

        if not args.reset and not anvio.DEBUG:
            self.is_database_exists()

        if args.reset:
            filesnpaths.gen_output_directory(self.pfam_data_dir, delete_if_exists=True, dont_warn=True)
        else:
            filesnpaths.gen_output_directory(self.pfam_data_dir)


        self.resolve_database_url()
        self.files = ['Pfam-A.hmm.gz', 'Pfam.version.gz', 'Pfam-A.clans.tsv.gz']


    def resolve_database_url(self):
        if self.args.pfam_version:
            page_index = 'releases/Pfam%s' % self.args.pfam_version
            self.run.info('Attempting to use version', self.args.pfam_version)
        else:
            page_index = 'current_release'
            self.run.info_single('No Pfam version specified. Using current release.')

        self.database_url = "http://ftp.ebi.ac.uk/pub/databases/Pfam/%s" % page_index


    def is_database_exists(self):
        if os.path.exists(os.path.join(self.pfam_data_dir, 'Pfam-A.hmm') or os.path.exists(os.path.join(self.pfam_data_dir, 'Pfam-A.hmm.gz'))):
            raise ConfigError("It seems you already have Pfam database installed in '%s', please use --reset flag if you want to re-download it." % self.pfam_data_dir)


    def get_remote_version(self):
        input_file = os.path.join(self.pfam_data_dir, '/Pfam.version.gz')
        if os.path.exists(input_file):
            try:
                with gzip.open(input_file,'rt') as f:
                    content = f.read()
            except gzip.BadGzipFile as err:
                raise Exception ("BadGzipFile", input_file)
        else:
            content = read_remote_file(self.database_url + '/Pfam.version.gz')

        # below we are parsing this, not so elegant.
        # Pfam release       : 31.0
        # Pfam-A families    : 16712
        # Date               : 2017-02
        # Based on UniProtKB : 2016_10

        version = content.strip().split('\n')[0].split(':')[1].strip()
        release_date = content.strip().split('\n')[2].split(':')[1].strip()

        self.run.info("Found Pfam version", "%s (%s)" % (version, release_date))

    def download(self, hmmpress_files=True):
        self.run.info("Database URL", self.database_url)

        for file_name in self.files:
            local_file = os.path.join(self.pfam_data_dir, file_name)
            if not os.path.exists(local_file):
                download_file(self.database_url + '/' + file_name,
                    os.path.join(self.pfam_data_dir, file_name), progress=self.progress, run=self.run)
            else:
                self.run.info("found local file", file_name)

        self.confirm_downloaded_files()
        self.decompress_files()
        if hmmpress_files:
            self.hmmpress_files()


    def confirm_downloaded_files(self):
        chksums_file = os.path.join(self.pfam_data_dir, 'md5_checksums')
        if os.path.exists(chksums_file): 
            self.run.info("found local file", 'md5_checksums')
            with open(chksums_file,'r') as fh:
                checksums_file = fh.read() 
        else:
            try:
                checksums_file = read_remote_file(self.database_url + '/md5_checksums', is_gzip=False).strip()
            except:
                self.run.warning("Checksum file '%s' is not available in FTP, Anvi'o won't be able to verify downloaded files." % (self.database_url + '/md5_checksums'))
                return

        checksums = {}
        for line in checksums_file.split('\n'):
            if not line: continue
            checksum, file_name = [item.strip() for item in line.strip().split()]
            checksums[file_name] = checksum

        for file_name in self.files:
            if not filesnpaths.is_file_exists(os.path.join(self.pfam_data_dir, file_name), dont_raise=True):
                raise ConfigError(f"Unfortunately, we failed to download the file {file_name}, please re-run setup "
                                  "with the --reset flag.")

            hash_on_disk = get_file_md5(os.path.join(self.pfam_data_dir, file_name))
            expected_hash = checksums[file_name]

            if not expected_hash == hash_on_disk:
                raise ConfigError(f"Please re-run setup with --reset, the file hash for {file_name} doesn't match to the hash "
                                  "we expected. If you continue to get this error after doing that, try removing the entire "
                                  f"Pfams data directory ({self.pfam_data_dir}) manually and running setup again (without the --reset flag).")


    def decompress_files(self):
        """Decompresses Pfam HMM profiles."""

        for file_name in self.files:
            full_path = os.path.join(self.pfam_data_dir, file_name)

            if full_path.endswith('.gz'):
                if not os.path.exists(full_path) and os.path.exists(full_path[:-3]):
                    self.run.warning("It seems the file at %s is already decompressed. You are probably seeing "
                                     "this message because Pfams was set up previously on this computer. Hakuna Matata. Anvi'o will "
                                     "simply skip decompressing this file at this time. But if you think there is an issue, you can "
                                     "re-do the Pfam setup by running `anvi-setup-pfams` again and using the --reset flag."
                                     % (full_path[:-3]))
                    continue
                elif not os.path.exists(full_path):
                    raise ConfigError("Oh no. The file at %s does not exist. Something is terribly wrong. :( Anvi'o suggests re-running "
                                      "`anvi-setup-pfams` using the --reset flag." % (full_path))
                gzip_decompress_file(full_path)
                os.remove(full_path)


    def hmmpress_files(self):
        """Runs hmmpress on Pfam HMM profiles."""

        for file_path in glob.glob(os.path.join(self.pfam_data_dir, '*.hmm')):
            cmd_line = ['hmmpress', file_path]
            log_file_path = os.path.join(self.pfam_data_dir, '00_hmmpress_log.txt')
            ret_val = run_command(cmd_line, log_file_path)

            if ret_val:
                raise ConfigError("Hmm. There was an error while running `hmmpress` on the Pfam HMM profiles. "
                                  "Check out the log file ('%s') to see what went wrong." % (log_file_path))
            else:
                # getting rid of the log file because hmmpress was successful
                os.remove(log_file_path)


class Pfam(object):
    def __init__(self, args, run=run, progress=progress):

        self.run = run
        self.progress = progress

        A = lambda x, t: t(args.__dict__[x]) if x in args.__dict__ else None
        null = lambda x: x
        self.contigs_db_path = A('contigs_db', null)
        self.num_threads = A('num_threads', null)
        self.hmm_program = A('hmmer_program', null) or 'hmmsearch'
        self.pfam_data_dir = A('pfam_data_dir', null)

        # load_catalog will populate this
        self.function_catalog = {}

        filesnpaths.is_program_exists(self.hmm_program)
        is_contigs_db(self.contigs_db_path)

        if not self.pfam_data_dir:
            self.pfam_data_dir = os.path.join(os.path.dirname(anvio.__file__), 'data/misc/Pfam')

        # here, in the process of checking whether Pfam has been downloaded into the pfam_data_dir,
        # we also decompress and hmmpress the profile if it is currently gzipped
        self.is_database_exists()

        self.run.info('Pfam database directory', self.pfam_data_dir)

        self.get_version()
        self.load_catalog()


    def is_database_exists(self):
        """Checks if database files exist and decompresses them if compressed

        This function verifies that pfam_data_dir contains the Pfam hmm profiles and checks whether
        they are compressed or not. If they are compressed, we decompress them and run hmmpress.
        """

        if not (os.path.exists(os.path.join(self.pfam_data_dir, 'Pfam-A.hmm.gz')) or os.path.exists(os.path.join(self.pfam_data_dir, 'Pfam-A.hmm'))):
            raise ConfigError("It seems you do not have Pfam database installed, please run 'anvi-setup-pfams' to download it.")

        # here we check if the HMM profile is compressed so we can decompress it for next time
        if os.path.exists(os.path.join(self.pfam_data_dir, 'Pfam-A.hmm.gz')):
            self.run.warning("Anvi'o has detected that your Pfam database is currently compressed. It will now be unpacked before "
                             "running HMMs.")
            gzip_decompress_file(os.path.join(self.pfam_data_dir, 'Pfam-A.hmm.gz'), keep_original=False)

            cmd_line = ['hmmpress', os.path.join(self.pfam_data_dir, 'Pfam-A.hmm')]
            log_file_path = os.path.join(self.pfam_data_dir, '00_hmmpress_log.txt')
            ret_val = run_command(cmd_line, log_file_path)

            if ret_val:
                raise ConfigError("Hmm. There was an error while running `hmmpress` on the Pfam HMM profiles. "
                                  "Check out the log file ('%s') to see what went wrong." % (log_file_path))
            else:
                # getting rid of the log file because hmmpress was successful
                os.remove(log_file_path)


    def get_version(self):
        with open(os.path.join(self.pfam_data_dir, 'Pfam.version')) as f:
            content = f.read()

        # below we are parsing this, not so elegant.
        # Pfam release       : 31.0
        # Pfam-A families    : 16712
        # Date               : 2017-02
        # Based on UniProtKB : 2016_10

        version = content.strip().split('\n')[0].split(':')[1].strip()
        release_date = content.strip().split('\n')[2].split(':')[1].strip()

        self.run.info("Pfam database version", "%s (%s)" % (version, release_date))


    def load_catalog(self):
        catalog_path = os.path.join(self.pfam_data_dir, 'Pfam-A.clans.tsv')
        self.function_catalog = get_TAB_delimited_file_as_dictionary(
            catalog_path,
            column_names=['accession', 'clan', 'unknown_column1', 'unknown_column2', 'function'],
            no_header=True
        )


    def get_function_from_catalog(self, accession, ok_if_missing_from_catalog=False):
        if '.' in accession:
            accession = accession.split('.')[0]

        if not accession in self.function_catalog:
            if ok_if_missing_from_catalog:
                return "Unknown function with PFAM accession %s" % accession
            else:
                raise ConfigError("It seems hmmscan/hmmsearch found an accession id that does not exists "
                                  "in Pfam catalog: %s" % accession)

        return self.function_catalog[accession]['function'] # maybe merge other columns too?


    def process(self):
        hmm_file = os.path.join(self.pfam_data_dir, 'Pfam-A.hmm')

        # initialize contigs database
        class Args: pass
        args = Args()
        args.contigs_db = self.contigs_db_path
        contigs_db = dbops.ContigsSuperclass(args)
        tmp_directory_path = filesnpaths.get_temp_directory_path()

        # get an instance of gene functions table
        gene_function_calls_table = TableForGeneFunctions(self.contigs_db_path, self.run, self.progress)

        # export AA sequences for genes
        target_files_dict = {'AA:GENE': os.path.join(tmp_directory_path, 'AA_gene_sequences.fa')}
        contigs_db.get_sequences_for_gene_callers_ids(output_file_path=target_files_dict['AA:GENE'],
                                                      simple_headers=True,
                                                      report_aa_sequences=True)

        # run hmmer
        hmmer = HMMer(target_files_dict, num_threads_to_use=self.num_threads, program_to_use=self.hmm_program)
        hmm_hits_file = hmmer.run_hmmer('Pfam', 'AA', 'GENE', None, None, len(self.function_catalog), hmm_file, None, '--cut_ga')

        if not hmm_hits_file:
            run.info_single("The HMM search returned no hits :/ So there is nothing to add to the contigs database. But "
                            "now anvi'o will add PFAMs as a functional source with no hits, clean the temporary directories "
                            "and gracefully quit.", nl_before=1, nl_after=1)
            shutil.rmtree(tmp_directory_path)
            hmmer.clean_tmp_dirs()
            gene_function_calls_table.add_empty_sources_to_functional_sources({'Pfam'})
            return

        # parse hmmer output
        parser = parser_modules['search']['hmmer_table_output'](hmm_hits_file, alphabet='AA', context='GENE', program=self.hmm_program)
        search_results_dict = parser.get_search_results()

        # add functions to database
        functions_dict = {}
        counter = 0
        for hmm_hit in search_results_dict.values():
            functions_dict[counter] = {
                'gene_callers_id': hmm_hit['gene_callers_id'],
                'source': 'Pfam',
                'accession': hmm_hit['gene_hmm_id'],
                'function': self.get_function_from_catalog(hmm_hit['gene_hmm_id'], ok_if_missing_from_catalog=True),
                'e_value': hmm_hit['e_value'],
            }

            counter += 1

        if functions_dict:
            gene_function_calls_table.create(functions_dict)
        else:
            self.run.warning("Pfam class has no hits to process. Returning empty handed, but still adding Pfam as "
                             "a functional source.")
            gene_function_calls_table.add_empty_sources_to_functional_sources({'Pfam'})

        if anvio.DEBUG:
            run.warning("The temp directories, '%s' and '%s' are kept. Please don't forget to clean those up "
                        "later" % (tmp_directory_path, ', '.join(hmmer.tmp_dirs)), header="Debug")
        else:
            run.info_single('Cleaning up the temp directory (you can use `--debug` if you would '
                            'like to keep it for testing purposes)', nl_before=1, nl_after=1)
            shutil.rmtree(tmp_directory_path)
            hmmer.clean_tmp_dirs()


class HMMProfile(object):
    """Reads, modifies, and writes hmm profiles (.hmm)"""

    def __init__(self, filepath, run=terminal.Run(), progress=terminal.Progress()):
        self.run = run
        self.progress = progress
        self.filepath = filepath

        filesnpaths.is_file_exists(self.filepath)

        self.data = {}
        self.load()


    def load(self):
        self.progress.new('HMM profile')
        self.progress.update('Loading %s' % self.filepath)

        with open(self.filepath) as f:
            for i, raw_profile in enumerate(get_chunk(f, separator='//\n', read_size=32768)):
                if not raw_profile.strip():
                    continue

                if i % 100 == 0:
                    self.progress.update('%d profiles loaded' % i)
                    self.progress.increment(increment_to=i)

                profile = self.process_raw_profile(raw_profile)
                self.data[profile['ACC']] = profile

        self.progress.end()
        self.run.info('HMM profiles loaded from', self.filepath)
        self.run.info('Number of HMM profiles loaded', '%d profiles' % len(self.data))


    def filter(self, by, subset):
        """Create a new .hmm that is a subset of the instantiated .hmm

        Parameters
        ==========
        by : str
            Which key would you like to filter by? For example, 'NAME', 'ACC', 'LENG'

        subset : set
            For each profile, if its corresponding value of `by` matches any of the elements in this
            set, it is included in the new .hmm file.
        """

        accessions_to_keep = set()
        num_profiles = len(self.data)

        self.progress.new('Filtering HMMs by %s' % by)

        i = 0
        for acc, profile in self.data.items():

            if i % 100 == 0:
                self.progress.update('%d/%d processed' % (i, num_profiles))
                self.progress.increment(increment_to=i)

            if profile[by] in subset:
                accessions_to_keep.add(acc)

            i += 1

        self.data = {acc: self.data[acc] for acc in accessions_to_keep}

        self.progress.end()
        self.run.info('Filtered by', by)
        self.run.info('Num subset values', len(subset))
        self.run.info('Profiles kept', len(self.data))


    def write(self, filepath):
        """Write self.data to an .hmm file.

        Parameters
        ==========
        filepath : str
            If explictly passed None, the source .hmm file will be overwritten.
        """

        if filepath is None:
            filepath = self.filepath

        with open(filepath, 'w') as f:
            for profile in self.data.values():
                f.write(profile['RAW'] + '//\n')

        self.run.info('.hmm written to', filepath)


    def __getitem__(self, key):
        return self.data[key]


    def process_raw_profile(self, raw_profile):
        """Process a profile.

        Processes a string that looks like this:

        HMMER3/f [3.2.1 | June 2018]
        NAME  Ubie_methyltran
        ACC   PF01209.17
        DESC  ubiE/COQ5 methyltransferase family
        LENG  233
        ALPH  amino
        RF    no

        ...

        STATS LOCAL FORWARD   -5.2458  0.70364
        HMM          A        C        D        E        F        G        H        I        K        L        M        N        P        Q        R        S        T        V        W        Y
                    m->m     m->i     m->d     i->m     i->i     d->m     d->d
          COMPO   2.54716  4.35070  2.86386  2.72305  3.17454  2.77327  3.72537  2.80280  2.64902  2.43657  3.49755  3.12309  3.54723  3.17117  3.00291  2.68310  2.88876  2.58381  4.57982  3.36830
                  2.68618  4.42225  2.77519  2.73123  3.46354  2.40513  3.72494  3.29354  2.67741  2.69355  4.24690  2.90347  2.73739  3.18146  2.89801  2.37887  2.77519  2.98518  4.58477  3.61503
                  0.24616  4.06597  1.60415  0.61958  0.77255  0.00000        *
              1   2.70153  4.53236  3.34187  2.81919  3.62674  3.53680  3.76783  2.88731  2.03395  2.51546  2.26844  3.23461  3.95064  3.01094  2.73876  2.86019  2.94221  2.70788  5.06922  3.81779      1 k - - X
                  2.68618  4.42225  2.77519  2.73123  3.46354  2.40513  3.72494  3.29354  2.67741  2.69355  4.24690  2.90347  2.73739  3.18146  2.89801  2.37887  2.77519  2.98518  4.58477  3.61503
                  0.03206  3.85187  4.57421  0.61958  0.77255  0.38135  1.14866
              2   2.61757  4.49839  3.28139  2.75711  3.69124  3.48770  3.79249  3.00826  2.63457  2.66714  2.81777  3.19513  3.91494  2.57012  2.99323  2.76655  1.92096  2.77781  5.10625  3.84713      2 t - - X
                  2.68618  4.42225  2.77519  2.73123  3.46354  2.40513  3.72494  3.29354  2.67741  2.69355  4.24690  2.90347  2.73739  3.18146  2.89801  2.37887  2.77519  2.98518  4.58477  3.61503
                  0.02580  4.06597  4.78831  0.61958  0.77255  0.48576  0.95510

                    ....

            232   2.59256  4.07223  4.02476  3.44210  3.10552  3.69012  3.96905  2.14938  3.31845  2.21385  3.16613  3.65707  4.05752  3.54613  3.50650  2.97719  2.45795  2.02283  3.32232  2.83766    259 v - - E
                  2.68618  4.42225  2.77519  2.73123  3.46354  2.40513  3.72494  3.29354  2.67741  2.69355  4.24690  2.90347  2.73739  3.18146  2.89801  2.37887  2.77519  2.98518  4.58477  3.61503
                  0.02580  4.06597  4.78831  0.61958  0.77255  0.48576  0.95510
            233   3.15895  5.12095  3.34717  2.93780  4.55191  3.62263  3.88116  3.99448  0.80088  3.51942  4.47451  3.35216  4.11266  3.06460  2.48782  3.19066  3.39982  3.69284  5.53001  4.37266    260 k - - E
                  2.68618  4.42225  2.77519  2.73123  3.46354  2.40513  3.72494  3.29354  2.67741  2.69355  4.24690  2.90347  2.73739  3.18146  2.89801  2.37887  2.77519  2.98518  4.58477  3.61503
                  0.01744  4.05761        *  0.61958  0.77255  0.00000        *
        //
        """

        # Init the profile, we will populate this
        profile = {
            'RAW': raw_profile,
            'MATCH_STATES': {},
            'NAME': None,
            'ACC': None,
            'DESC': None,
            'LENG': None,
            'ALPH': None,
            'MAP': None,
            'CONS': None,
        }

        # Parse the header
        FIELDS = {
            'NAME': lambda x: x,
            'ACC': lambda x: x,
            'DESC': lambda x: x,
            'LENG': lambda x: int(x),
            'ALPH': lambda x: x,
            'MAP': lambda x: True if x == 'yes' else False,
            'CONS': lambda x: True if x == 'yes' else False,
        }

        profile_split = raw_profile.split()
        for i, TAG in enumerate(profile_split):
            if TAG in FIELDS:
                profile[TAG] = FIELDS[TAG](profile_split[i+1])
                FIELDS.pop(TAG)

            if TAG == 'HMM':
                # We are done parsing the header
                break

        if len(FIELDS):
            raise ConfigError("Anvi'o was expecting each profile HMM to have the following fields, however this "
                              "profile did not: %s. Here is the raw profile: %s" % (', '.join(FIELDS), raw_profile))

        if profile['ACC'].find('.') != -1:
            profile['ACC'], profile['VERSION'] = profile['ACC'].split('.')

        # Now we go line by line
        raw_lines = raw_profile.split('\n')
        line_no = 0

        # Get the alphabet and markov state
        for line in raw_lines:
            if not line.startswith('HMM '):
                line_no += 1
                continue

            alphabet = line.split()[1:]
            states = raw_lines[line_no+1].split()[1:]

            profile['alphabet'] = alphabet
            profile['states'] = states

            line_no += 2
            break

        # Find the start
        for line in raw_lines[line_no:]:
            if not line.split()[0] == '1':
                line_no += 1
                continue
            break

        profile['MATCH_STATES'] = {
            'MATCH_STATE': [],
            'CS': [],
            'MM': [],
            'RF': [],
            'CONS': [],
            'MAP': [],
            'IC': [], # information content (https://en.wikipedia.org/wiki/Sequence_logo)
        }

        for i in range(line_no, len(raw_lines[line_no:]), 3):
            emission_line, insertion_line, state_line = raw_lines[i:i+3]
            emission = emission_line.split()

            # These are not used currently but maybe someday will be
            insertion = insertion_line.split()
            state = state_line.split()

            assign_type = lambda x, t: t(x) if x != '-' else '-'
            profile['MATCH_STATES']['MATCH_STATE'].append(int(emission.pop(0)))
            profile['MATCH_STATES']['CS'].append(assign_type(emission.pop(-1), str))
            profile['MATCH_STATES']['MM'].append(assign_type(emission.pop(-1), str))
            profile['MATCH_STATES']['RF'].append(assign_type(emission.pop(-1), str))
            profile['MATCH_STATES']['CONS'].append(assign_type(emission.pop(-1).upper(), str))
            profile['MATCH_STATES']['MAP'].append(assign_type(emission.pop(-1), int))

            # Store the information content of the match state
            emission_probs = np.exp(-np.array([0 if x == '*' else float(x) for x in emission]))
            information_content = np.log2(len(profile['alphabet'])) - entropy(emission_probs, base=2)
            profile['MATCH_STATES']['IC'].append(information_content)

        # Cast match state info as a dataframe
        profile['MATCH_STATES'] = pd.DataFrame(profile['MATCH_STATES'])

        # We 0-index the match states because we are rebels
        profile['MATCH_STATES']['MATCH_STATE'] -= 1
        profile['MATCH_STATES'].set_index('MATCH_STATE', drop=True, inplace=True)

        return profile

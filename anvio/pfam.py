#!/usr/bin/env python
# -*- coding: utf-8
"""
    This file contains PfamSetup and Pfam classes.

"""
import os
import gzip
import shutil
import requests
from io import BytesIO

import anvio
import anvio.dbops as dbops
import anvio.utils as utils
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths

from anvio.drivers.hmmer import HMMer
from anvio.parsers import parser_modules
from anvio.errors import ConfigError
from anvio.tables.genefunctions import TableForGeneFunctions


__author__ = "Developers of anvi'o (see AUTHORS.txt)"
__copyright__ = "Copyleft 2015-2018, the Meren Lab (http://merenlab.org/)"
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "Özcan Esen"
__email__ = "ozcanesen@gmail.com"


run = terminal.Run()
progress = terminal.Progress()
pp = terminal.pretty_print


def read_remote_file(url, is_gzip=True):
    remote_file = requests.get(url)

    if is_gzip:
        buf = BytesIO(remote_file.content)
        fg = gzip.GzipFile(fileobj=buf)
        return fg.read().decode('utf-8')

    return remote_file.content.decode('utf-8')


class PfamSetup(object):
    def __init__(self, args, run=run, progress=progress):
        self.args = args
        self.run = run
        self.progress = progress
        self.pfam_data_dir = args.pfam_data_dir

        filesnpaths.is_program_exists('hmmpress')

        if not self.pfam_data_dir:
            self.pfam_data_dir = os.path.join(os.path.dirname(anvio.__file__), 'data/misc/Pfam')

        if not args.reset:
            self.is_database_exists()

        filesnpaths.gen_output_directory(self.pfam_data_dir, delete_if_exists=args.reset)

        self.database_url = "http://ftp.ebi.ac.uk/pub/databases/Pfam/current_release"
        self.files = ['Pfam-A.hmm.gz', 'Pfam.version.gz', 'Pfam-A.clans.tsv.gz']


    def is_database_exists(self):
        if os.path.exists(os.path.join(self.pfam_data_dir, 'Pfam-A.hmm.gz')):
            raise ConfigError("It seems you already have Pfam database installed in '%s', please use --reset flag if you want to re-download it." % self.pfam_data_dir)


    def get_remote_version(self):
        content = read_remote_file(self.database_url + '/Pfam.version.gz')

        # below we are parsing this, not so elegant.
        # Pfam release       : 31.0
        # Pfam-A families    : 16712
        # Date               : 2017-02
        # Based on UniProtKB : 2016_10

        version = content.strip().split('\n')[0].split(':')[1].strip()
        release_date = content.strip().split('\n')[2].split(':')[1].strip()

        self.run.info("Current Pfam version on EBI", "%s (%s)" % (version, release_date))


    def download(self):
        self.run.info("Database URL", self.database_url)

        for file_name in self.files:
            utils.download_file(self.database_url + '/' + file_name,
                os.path.join(self.pfam_data_dir, file_name), progress=self.progress, run=self.run)

        self.confirm_downloaded_files()
        self.decompress_files()


    def confirm_downloaded_files(self):
        checksums_file = read_remote_file(self.database_url + '/md5_checksums', is_gzip=False).strip()
        checksums = {}

        for line in checksums_file.split('\n'):
            checksum, file_name = [item.strip() for item in line.strip().split()]
            checksums[file_name] = checksum

        for file_name in self.files:
            if not filesnpaths.is_file_exists(os.path.join(self.pfam_data_dir, file_name), dont_raise=True):
                 # TO DO: Fix messages :(
                raise ConfigError("Have missing file %s, please run --reset" % file_name)

            hash_on_disk = utils.get_file_md5(os.path.join(self.pfam_data_dir, file_name))
            expected_hash = checksums[file_name]

            if not expected_hash == hash_on_disk:
                # TO DO: Fix messages :(
                raise ConfigError("Please run with --reset, one file hash doesn't match. %s" % file_name)


    def decompress_files(self):
        # Decompressing Pfam-A.hmm.gz is not necessary, HMMer class works with .gz

        for file_name in ['Pfam.version.gz', 'Pfam-A.clans.tsv.gz']:
            full_path = os.path.join(self.pfam_data_dir, file_name)

            utils.gzip_decompress_file(full_path)
            os.remove(full_path)


class Pfam(object):
    def __init__(self, args, run=run, progress=progress):
        self.args = args
        self.run = run
        self.progress = progress
        self.contigs_db_path = args.contigs_db
        self.num_threads = args.num_threads
        self.pfam_data_dir = args.pfam_data_dir

        # load_catalog will populate this
        self.function_catalog = {}

        filesnpaths.is_program_exists('hmmscan')
        utils.is_contigs_db(self.contigs_db_path)

        if not self.pfam_data_dir:
            self.pfam_data_dir = os.path.join(os.path.dirname(anvio.__file__), 'data/misc/Pfam')

        self.is_database_exists()

        self.run.info('Pfam database directory', self.pfam_data_dir)

        self.get_version()
        self.load_catalog()


    def is_database_exists(self):
        if not os.path.exists(os.path.join(self.pfam_data_dir, 'Pfam-A.hmm.gz')):
            raise ConfigError("It seems you do not have Pfam database installed, please run 'anvi-setup-pfams' to download it.")


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
        self.function_catalog = utils.get_TAB_delimited_file_as_dictionary(catalog_path,
            column_names=['accession', 'clan', 'unknown_column1', 'unknown_column2', 'function'])


    def get_function_from_catalog(self, accession):
        if '.' in accession:
            accession = accession.split('.')[0]

        if not accession in self.function_catalog:
            # TO DO: messsages
            raise ConfigError("It seems hmmscan found a accession id that does not exists in Pfam catalog, Id: %s" % accession)

        return self.function_catalog[accession]['function'] # maybe merge other columns too?


    def process(self):
        hmm_file = os.path.join(self.pfam_data_dir, 'Pfam-A.hmm.gz')

        # initialize contigs database
        class Args: pass
        args = Args()
        args.contigs_db = self.contigs_db_path
        contigs_db = dbops.ContigsSuperclass(args)
        tmp_directory_path = filesnpaths.get_temp_directory_path()

        # export AA sequences for genes
        target_files_dict = {'AA:GENE': os.path.join(tmp_directory_path, 'AA_gene_sequences.fa')}
        contigs_db.gen_FASTA_file_of_sequences_for_gene_caller_ids(output_file_path=target_files_dict['AA:GENE'],
                                                                   simple_headers=True,
                                                                   rna_alphabet=False,
                                                                   report_aa_sequences=True)

        # run hmmscan
        hmmer = HMMer(target_files_dict, num_threads_to_use=self.num_threads)
        hmm_hits_file = hmmer.run_hmmscan('Pfam', 'AA', 'GENE', None, None, len(self.function_catalog), hmm_file, None, '--cut_ga')

        # parse hmmscan output
        parser = parser_modules['search']['hmmscan'](hmm_hits_file, alphabet='AA', context='GENE')
        search_results_dict = parser.get_search_results()

        # add functions to database
        functions_dict = {}
        counter = 0
        for hmm_hit in search_results_dict.values():
            functions_dict[counter] = {
                'gene_callers_id': hmm_hit['gene_callers_id'],
                'source': 'Pfam',
                'accession': hmm_hit['gene_hmm_id'],
                'function': self.get_function_from_catalog(hmm_hit['gene_hmm_id']),
                'e_value': hmm_hit['e_value'],
            }

            counter += 1

        gene_function_calls_table = TableForGeneFunctions(self.contigs_db_path, self.run, self.progress)
        gene_function_calls_table.create(functions_dict)

        if anvio.DEBUG:
            run.warning("The temp directories, '%s' and '%s' are kept. Please don't forget to clean those up\
                         later" % (tmp_directory_path, ', '.join(hmmer.tmp_dirs)), header="Debug")
        else:
            run.info_single('Cleaning up the temp directory (you can use `--debug` if you would\
                             like to keep it for testing purposes)', nl_before=1, nl_after=1)
            shutil.rmtree(tmp_directory_path)
            hmmer.clean_tmp_dirs()

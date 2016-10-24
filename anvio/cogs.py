# -*- coding: utf-8
# pylint: disable=line-too-long
"""
    Making sense of COGs.
"""

import os
import sys
import gzip
import shutil

import anvio
import anvio.utils as utils
import anvio.dictio as dictio
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths

from anvio.drivers.diamond import Diamond
from anvio.drivers.blast import BLAST
from anvio.errors import ConfigError

# just to make sure things don't break too far when they do:
COGs_DATA_VERSION='1'


__author__ = "A. Murat Eren"
__copyright__ = "Copyright 2016, The anvio Project"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "A. Murat Eren"
__email__ = "a.murat.eren@gmail.com"


run = terminal.Run()
progress = terminal.Progress()
pp = terminal.pretty_print

J = lambda x, y: os.path.join(x, y)

COGs_PATH = os.path.join(os.path.dirname(anvio.__file__), 'data/misc/COGs/COGs.txt')
CATEGORIES_PATH = os.path.join(os.path.dirname(anvio.__file__), 'data/misc/COGs/CATEGORIES.txt')


class COGsData:
    """A class to make sense of COG ids and categories"""
    def __init__(self, COGs_path=COGs_PATH, categories_path=CATEGORIES_PATH):
        filesnpaths.is_file_tab_delimited(COGs_path)
        filesnpaths.is_file_tab_delimited(categories_path)

        self.cogs = utils.get_TAB_delimited_file_as_dictionary(COGs_path, no_header=True, column_names=['COG', 'categories', 'annotation'])
        self.categories = utils.get_TAB_delimited_file_as_dictionary(categories_path, no_header=True, column_names=['category', 'description'])

        for cog in self.cogs:
            self.cogs[cog]['categories'] = [c.strip() for c in self.cogs[cog]['categories'].split(',')]

        for cat in self.categories:
            self.categories[cat] = self.categories[cat]['description']


class SetupCOGs:
    """A class to download and setup the COG data from NCBI."""
    def __init__(self, args, run=run, progress=progress):
        self.COG_data_dir = J(os.path.dirname(anvio.__file__), 'data/misc/COGs')
        self.COG_data_dir_version = J(self.COG_data_dir, '.VERSION')
        self.raw_NCBI_files_dir = J(self.COG_data_dir, 'RAW_DATA_FROM_NCBI')

        self.run = run
        self.progress = progress

        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        self.num_threads = A('num_threads') or 1
        self.reset = A('reset')

        if not os.path.exists(self.COG_data_dir):
            try:
                os.mkdir(self.COG_data_dir)
                open(self.COG_data_dir_version, 'w').write(COGs_DATA_VERSION)
            except Exception, e:
                raise ConfigError, "So the COGs data directory is not there, and anvi'o wants to create one. But it didn't\
                                    go that well. It could be due to permissions (which may require you to run this with sudo\
                                    or may need to ask your sys admin to do it for you since this is a one time operation), or\
                                    it could be due to something totally irrelevant. Here is the error message: '%s'" % e

        filesnpaths.is_output_dir_writable(self.COG_data_dir)

        self.raw_files = {
                #'cog2003-2014.csv': {
                #    'url': 'ftp://ftp.ncbi.nih.gov/pub/COG/COG2014/data/cog2003-2014.csv',
                #    'func': self.format_p_id_to_cog_id_cPickle,
                #    'formatted_file_name': 'PID-TO-CID.cPickle'},
                #'cognames2003-2014.tab': {
                #    'url': 'ftp://ftp.ncbi.nih.gov/pub/COG/COG2014/data/cognames2003-2014.tab',
                #    'func': self.format_cog_names,
                #    'formatted_file_name': 'COG.txt'},
                #'fun2003-2014.tab': {
                #    'url': 'ftp://ftp.ncbi.nih.gov/pub/COG/COG2014/data/fun2003-2014.tab',
                #    'func': self.format_categories,
                #    'formatted_file_name': 'CATEGORIES.txt'},
                'prot2003-2014.fa.gz': {
                    'url': 'ftp://ftp.ncbi.nih.gov/pub/COG/COG2014/data/prot2003-2014.fa.gz',
                    'func': self.format_protein_db,
                    'formatted_file_name': 'IGNORE_THIS_AND_SEE_THE_FUNCTION'}
        }


    def process(self):
        run.info('COG data dir', self.COG_data_dir)

        if self.reset:
            run.warning('This program will remove everything in the COGs data directory, then download and reformat\
                         everything from scratch. Press ENTER to continue, or press CTRL + C to cancel if you are\
                         not OK with that.')
            self.wait_for_the_user()

            # OK. reset the crap out of it.
            shutil.rmtree(self.COG_data_dir)
            os.mkdir(self.COG_data_dir)
            open(self.COG_data_dir_version, 'w').write(COGs_DATA_VERSION)
        else:
            run.warning("This program will first check whether you have all the raw files, and then will attempt to\
                         regenerate everything that is necessary from them. Press ENTER to continue, or\
                         press CTRL + C to cancel if it doesn't sound right to you.")
            self.wait_for_the_user()

        if not os.path.exists(self.COG_data_dir_version) or open(self.COG_data_dir_version).read() != COGs_DATA_VERSION:
            raise ConfigError, "The version of your COGs data directory is different than what anvi'o hoping to see.\
                                It seems you need to (re)run anvi'o script to download and format COG data from NCBI."

        # get raw files
        self.get_raw_data()

        # format raw files
        self.setup_raw_data()


    def format_p_id_to_cog_id_cPickle(self, input_file_path, output_file_path):
        progress.new('Formatting protein ids to COG ids file')
        progress.update('...')

        p_id_to_cog_id = {}

        for line in open(input_file_path, 'rU').readlines():
            fields = line.strip('\n').split(',')
            p_id = fields[0]
            COG = fields[6]

            if p_id in p_id_to_cog_id:
                p_id_to_cog_id[p_id].append(COG)
            else:
                p_id_to_cog_id[p_id] = [COG]

        dictio.write_serialized_object(p_id_to_cog_id, output_file_path)

        progress.end()


    def format_cog_names(self, input_file_path, output_file_path):
        progress.new('Formatting COG names file')
        progress.update('...')
        output = open(output_file_path, 'w')
        for line in open(input_file_path, 'rU').readlines():
            if line.startswith('#'):
                continue
            COG, function, name = line.strip('\n').split('\t')

            # get rid of non-ascii chars:
            name = ''.join([i if ord(i) < 128 else '' for i in name])

            output.write('\t'.join([COG, ', '.join(list(function)), name]) + '\n')
            progress.end()


    def format_categories(self, input_file_path, output_file_path):
        progress.new('Formatting COG categories file')
        progress.update('...')

        output = open(output_file_path, 'w')
        for line in open(input_file_path, 'rU').readlines():
            if line.startswith('#'):
                continue

            category, description = line.strip('\n').split('\t')

            # get rid of non-ascii chars:
            description = ''.join([i if ord(i) < 128 else '' for i in description])

            output.write('\t'.join([category, '[%s] %s' % (category, description)]) + '\n')

        progress.end()


    def format_protein_db(self, input_file_path, output_file_path):
        progress.new('Formatting raw files')
        progress.update('Decompressing protein sequences')

        # poor man's uncompress
        temp_fasta_path = filesnpaths.get_temp_file_path()
        with open(temp_fasta_path, 'w') as f_out, gzip.open(input_file_path, 'rb') as f_in:
            f_out.write(f_in.read())

        progress.end()

        if utils.is_program_exists('diamond', dont_raise=True):
            output_dir = J(self.COG_data_dir, 'DB_DIAMOND')
            if os.path.exists(output_dir):
                shutil.rmtree(output_dir)

            os.mkdir(output_dir)

            output_db_path = J(output_dir, 'COGs')
            log_file_path = J(output_dir, 'log.txt')

            self.run.info('Diamond log', log_file_path)

            diamond = Diamond(temp_fasta_path)
            diamond.num_threads = self.num_threads
            diamond.run.log_file_path = log_file_path
            diamond.target_db_path = output_db_path
            diamond.makedb()
        else:
            self.run.warning("Diamond does not seem to be installed on this system, so anvi'o is not going to\
                              generate a search database for it. Remember this when/if things go South.")

        if utils.is_program_exists('makeblastdb', dont_raise=True) and utils.is_program_exists('blastp', dont_raise=True):
            output_dir = J(self.COG_data_dir, 'DB_BLAST')
            if os.path.exists(output_dir):
                shutil.rmtree(output_dir)

            os.mkdir(output_dir)

            output_db_path = J(output_dir, 'COGs')
            log_file_path = J(output_dir, 'log.txt')

            self.run.info('BLAST log', log_file_path)

            blast = BLAST(temp_fasta_path)
            blast.target_db_path = output_db_path
            blast.run.log_file_path = log_file_path
            blast.num_threads = self.num_threads
            blast.makedb()
        else:
            self.run.warning("BLAST tools do not seem to be installed on this system, so anvi'o is not going to\
                              generate a search database for them to be used. Keep this in mind for later.")

        os.remove(temp_fasta_path)


    def get_raw_data(self):
        if not os.path.exists(self.raw_NCBI_files_dir):
            os.mkdir(self.raw_NCBI_files_dir)
            open(self.COG_data_dir_version, 'w').write(COGs_DATA_VERSION)

        for file_name in self.raw_files:
            file_path = J(self.raw_NCBI_files_dir, file_name)
            if not os.path.exists(file_path):
                utils.download_file(self.raw_files[file_name]['url'], file_path, progress=progress, run=run)


    def setup_raw_data(self):
        for file_name in self.raw_files:
            file_path = J(self.raw_NCBI_files_dir, file_name)

            if not os.path.exists(file_path):
                raise ConfigError, "Something is wrong :/ Raw files are not in place..."

            self.raw_files[file_name]['func'](file_path, J(self.COG_data_dir, self.raw_files[file_name]['formatted_file_name']))


    def wait_for_the_user(self):
        try:
            raw_input("Press ENTER to continue...\n")
        except:
            sys.exit()



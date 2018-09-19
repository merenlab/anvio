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
import anvio.dbops as dbops
import anvio.dictio as dictio
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths

from anvio.drivers.diamond import Diamond
from anvio.drivers.blast import BLAST
from anvio.errors import ConfigError
from anvio.tables.genefunctions import TableForGeneFunctions

# just to make sure things don't break too far when they do:
COG_DATA_VERSION='2'


__author__ = "Developers of anvi'o (see AUTHORS.txt)"
__copyright__ = "Copyleft 2015-2018, the Meren Lab (http://merenlab.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "A. Murat Eren"
__email__ = "a.murat.eren@gmail.com"


run = terminal.Run()
progress = terminal.Progress()
pp = terminal.pretty_print

J = lambda x, y: os.path.join(x, y)


class Args():
    pass


class COGs:
    """A class to run COGs"""
    def __init__(self, args=Args(), run=run, progress=progress):
        self.args = args
        self.run = run
        self.progress = progress

        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        self.num_threads = A('num_threads')
        self.contigs_db_path = A('contigs_db')
        self.search_with = A('search_with') or 'diamond'
        self.temp_dir_path = A('temporary_dir_path')
        self.sensitive = A('sensitive')

        self.log_file_path = None

        self.default_search_method = 'diamond'
        self.search_methods_factory = {'diamond': self.search_with_diamond,
                                       'blastp': self.search_with_ncbi_blast}
        self.available_search_methods = [p for p in self.search_methods_factory.keys() if utils.is_program_exists(p, dont_raise=True)]

        if not len(self.available_search_methods):
            raise ConfigError("None of the serach methods this class could use, which include '%s', seem to be\
                               available on your system :/" % (', '.join(list(self.search_methods_factory.keys()))))

        if self.default_search_method not in self.available_search_methods:
            self.default_search_method = self.available_search_methods[0]


        self.hits = None # the search function will take care of this one.

        if len(args.__dict__):
            self.COG_setup = COGsSetup(args)
            self.COG_data_dir = self.COG_setup.COG_data_dir
            self.available_db_search_program_targets = self.COG_setup.get_formatted_db_paths()
            self.essential_files = self.COG_setup.get_essential_file_paths()


    def process(self, aa_sequences_file_path=None):
        if self.search_with not in self.available_search_methods:
            raise ConfigError("Let us start by making it clear that we probably like '%s' as much as you do, but it doesn't\
                               seem to be available on your system OR recognized by the COGs class since anvi'o couldn't\
                               find it among the available search methods. You probably need to try something else :/" \
                                                                                                    % self.search_with)

        if self.search_with not in self.available_db_search_program_targets:
            raise ConfigError("Anvi'o understands that you want to use '%s' to search for COGs, however, there is no\
                               database formatted under the COGs data directory for that program :/ You may need to\
                               re-run the COGs setup, UNLESS, you set up your COG data directory somewhere else than what\
                               anvi'o attempts to use at the moment ('%s'). If that is the case, this may be the best\
                               time to point the right directory using the --cog-data-dir parameter, or the environmental\
                               variable 'ANVIO_COG_DATA_DIR'." % (self.search_with, self.COG_data_dir))

        if not aa_sequences_file_path and not self.contigs_db_path:
            raise ConfigError("You either need to provide an anvi'o contigs database path, or a FASTA file for AA\
                               sequences")

        if aa_sequences_file_path and self.contigs_db_path:
            raise ConfigError("You can't provide both an AA sequences file and a contigs database. Choose one!")

        if self.contigs_db_path:
            utils.is_contigs_db(self.contigs_db_path)

        if not self.temp_dir_path:
            self.temp_dir_path = filesnpaths.get_temp_directory_path()
            self.remove_temp_dir_path = True
        else:
            filesnpaths.is_file_exists(self.temp_dir_path)
            filesnpaths.is_output_dir_writable(self.temp_dir_path)

            self.run.warning("Because you set the temporary directory path by hand, anvi'o will not remove its content\
                              when it is done. But she certainly hopes that you will clean those files later.")

            self.remove_temp_dir_path = False

        self.run.info('COG data directory', self.COG_data_dir)
        self.run.info('Searching with', self.search_with)
        self.run.info('Directory to store temporary files', self.temp_dir_path)
        self.run.info('Directory will be removed after the run', self.remove_temp_dir_path)

        if not aa_sequences_file_path:
            aa_sequences_file_path = dbops.export_aa_sequences_from_contigs_db(self.contigs_db_path, J(self.temp_dir_path, 'aa_sequences.fa'))

        # do the search
        search_results_tabular = self.search_methods_factory[self.search_with](aa_sequences_file_path)

        # convert the output to a hits dict
        self.hits = utils.get_BLAST_tabular_output_as_dict(search_results_tabular, target_id_parser_func=lambda x: x.split('|')[1])

        # store hits into the contigs database
        self.store_hits_into_contigs_db()

        if self.remove_temp_dir_path:
            shutil.rmtree(self.temp_dir_path)


    def store_hits_into_contigs_db(self):
        if not self.hits:
            self.run.warning("COGs class has no hits to process. Returning empty handed.")
            return

        cogs_data = COGsData(self.args)
        cogs_data.init_p_id_to_cog_id_dict()

        functions_dict = {}
        self.__entry_id = 0


        def add_entry(gene_callers_id, source, accession, function, e_value):
            functions_dict[self.__entry_id] = {'gene_callers_id': int(gene_callers_id),
                                               'source': source,
                                               'accession': accession,
                                               'function': function,
                                               'e_value': float(e_value)}
            self.__entry_id += 1

        # let's keep track of hits that match to missing COGs
        hits_for_missing_cogs = 0
        missing_cogs_found = set([])

        for gene_callers_id in self.hits:
            ncbi_protein_id = self.hits[gene_callers_id]['hit']

            in_proteins_FASTA_not_in_cogs_CSV = []
            if ncbi_protein_id not in cogs_data.p_id_to_cog_id:
                in_proteins_FASTA_not_in_cogs_CSV.append((ncbi_protein_id, gene_callers_id),)
            else:
                COG_ids = cogs_data.p_id_to_cog_id[ncbi_protein_id]

            annotations = []
            categories = set([])
            for COG_id in COG_ids:
                # is missing?
                if COG_id in cogs_data.missing_cogs:
                    missing_cogs_found.add(COG_id)
                    hits_for_missing_cogs += 1
                    continue

                # resolve categories
                for category in cogs_data.cogs[COG_id]['categories']:
                    categories.add(category)

                # append annotation
                annotations.append(cogs_data.cogs[COG_id]['annotation'])

            # all these shitty heuristics... If there are multiple COG ids or categories, separate them from each other by '!!!' so parsing
            # them later is possible. Am I embarrassed? Yes. Is there a better way of doing this efficiently? Absolutely. What time is it?
            # 9pm. Where am I? In the lab. Is it OK for me to let this slip away if it means for me to go home sooner? Yes, probably. Am I
            # gonna remember this crap in the code for the next two months at random times in the shower and feel bad about myself? Fuck yes.
            add_entry(gene_callers_id, 'COG_FUNCTION', '!!!'.join(COG_ids), '!!!'.join(annotations), self.hits[gene_callers_id]['evalue'])
            add_entry(gene_callers_id, 'COG_CATEGORY', '!!!'.join(categories), '!!!'.join(categories), 0.0)

        # store hits in contigs db.
        gene_function_calls_table = TableForGeneFunctions(self.contigs_db_path, self.run, self.progress)
        gene_function_calls_table.create(functions_dict)

        if len(missing_cogs_found):
            self.run.warning('Although your COGs are successfully added to the database, there were some COG IDs your genes hit\
                              were among the ones that were not described in the raw data. Here is the list of %d COG IDs that\
                              were hit %d times: %s.' % (len(missing_cogs_found), hits_for_missing_cogs, ', '.join(missing_cogs_found)))

        if len(in_proteins_FASTA_not_in_cogs_CSV):
            # so some of the hits represented in the FASTA file from the NCBI were not put in the
            # CSV file from NCBI to associate them with COGs
            report_output_file_path = filesnpaths.get_temp_file_path()
            report_output = open(report_output_file_path, 'w')
            report_output.write('anvio_gene_callers_id\tNCBI_protein_id\n')

            for protein_id, gene_callers_id in in_proteins_FASTA_not_in_cogs_CSV:
                report_output.write('%s\t%s\n' % (gene_callers_id, protein_id))

            report_output.close()

            self.run.warning("This is important. %s hits for your genes that appeared in the proteins FASTA file from the NCBI had protein\
                              IDs that were not described in the CSV file from the NCBI that associates each protein ID with a COG function.\
                              That's OK if you don't care. But if you would like to take a look, anvi'o stored a report\
                              file for you at %s" \
                        % (len(in_proteins_FASTA_not_in_cogs_CSV), report_output_file_path))


    def search_with_diamond(self, aa_sequences_file_path):
        diamond = Diamond(aa_sequences_file_path, run=self.run, progress=self.progress, num_threads=self.num_threads)

        diamond.target_fasta = self.available_db_search_program_targets['diamond']
        self.run.log_file_path = self.log_file_path or J(self.temp_dir_path, 'log.txt')
        diamond.search_output_path = J(self.temp_dir_path, 'diamond-search-results')
        diamond.tabular_output_path = J(self.temp_dir_path, 'diamond-search-results.txt')

        diamond.sensitive = self.sensitive
        diamond.max_target_seqs = 1

        diamond.blastp()
        diamond.view()

        return diamond.tabular_output_path


    def search_with_ncbi_blast(self, aa_sequences_file_path):
        blast = BLAST(aa_sequences_file_path, run=self.run, progress=self.progress, num_threads=self.num_threads)

        blast.target_fasta = self.available_db_search_program_targets['blastp']
        self.run.log_file_path = self.log_file_path or J(self.temp_dir_path, 'log.txt')
        blast.search_output_path = J(self.temp_dir_path, 'blast-search-results.txt')
        blast.max_target_seqs = 1

        blast.blast()

        return blast.search_output_path


class COGsData:
    """A class to make sense of COG ids and categories"""
    def __init__(self, args=Args(), cog_data_dir = None, run=run, progress=progress, panic_on_failure_to_init=False):
        self.run = run
        self.progress = progress

        if cog_data_dir:
            args.cog_data_dir = cog_data_dir

        self.setup = COGsSetup(args)
        self.essential_files = self.setup.get_essential_file_paths()

        self.cogs = None
        self.p_id_to_cog_id = None
        self.missing_cogs = None
        self.categories = None
        self.initialized = False

        if self.essential_files:
            self.init()
        elif panic_on_failure_to_init:
            raise ConfigError("It seems you don't have your COG data set up on this system. Whatever you were\
                                trying to do is not going to continue being done :( Did you setup your COGs? If\
                                not, you can take a look at the program `anvi-setup-ncbi-cogs`. Maybe you did\
                                setup into another directory than the default destination? If that is the case\
                                maybe you can use the `--cog-data-dir` parameter if it is applicable? No? None\
                                of these work? Well. Anvi'o hates it as much as you do when things come to this.")


    def init(self):
        self.progress.new('Initializing COGs Data')
        self.progress.update('Reading COG functions ...')
        self.cogs = utils.get_TAB_delimited_file_as_dictionary(self.essential_files['COG.txt'], no_header=True, column_names=['COG', 'categories', 'annotation'])

        self.progress.update('Reading COG categories ...')
        self.categories = utils.get_TAB_delimited_file_as_dictionary(self.essential_files['CATEGORIES.txt'], no_header=True, column_names=['category', 'description'])

        self.progress.update('Reading missing COG IDs ...')
        self.missing_cogs = dictio.read_serialized_object(self.essential_files['MISSING_COG_IDs.cPickle'])

        self.progress.end()

        for cog in self.cogs:
            self.cogs[cog]['categories'] = [c.strip() for c in self.cogs[cog]['categories'].split(',')]

        for cat in self.categories:
            self.categories[cat] = self.categories[cat]['description']

        self.initialized = True


    def init_p_id_to_cog_id_dict(self):
        self.progress.new('Initializing COGs Data')
        self.progress.update('Reading NCBI Protein ID to COG id converter ...')

        self.p_id_to_cog_id = dictio.read_serialized_object(self.essential_files['PID-TO-CID.cPickle'])

        self.progress.end()


class COGsSetup:
    """A class to download and setup the COG data from NCBI."""
    def __init__(self, args=Args(), cog_data_dir = None, run=run, progress=progress):
        self.run = run
        self.progress = progress

        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        self.num_threads = A('num_threads') or 1
        self.just_do_it = A('just_do_it')
        self.reset = A('reset')
        self.cog_data_source = 'unknown'

        if cog_data_dir:
            self.COG_data_dir = cog_data_dir
            self.cog_data_source = 'The function call.'
        elif A('cog_data_dir'):
            self.COG_data_dir = A('cog_data_dir')
            self.cog_data_source = 'The command line parameter.'
        elif 'ANVIO_COG_DATA_DIR' in os.environ:
            self.COG_data_dir = os.environ['ANVIO_COG_DATA_DIR']
            self.cog_data_source = 'The environmental variable.'
        else:
            self.COG_data_dir = J(os.path.dirname(anvio.__file__), 'data/misc/COG')
            self.cog_data_source = "The anvi'o default."

        self.COG_data_dir = os.path.abspath(os.path.expanduser(self.COG_data_dir))

        self.COG_data_dir_version = J(self.COG_data_dir, '.VERSION')
        self.raw_NCBI_files_dir = J(self.COG_data_dir, 'RAW_DATA_FROM_NCBI')

        self.run.info('COG data directory', self.COG_data_dir)
        self.run.info('COG data source', self.cog_data_source)

        self.files = {
                'cog2003-2014.csv': {
                    'url': 'ftp://ftp.ncbi.nih.gov/pub/COG/COG2014/data/cog2003-2014.csv',
                    'func': self.format_p_id_to_cog_id_cPickle,
                    'type': 'essential',
                    'formatted_file_name': 'PID-TO-CID.cPickle'},
                'cognames2003-2014.tab': {
                    'url': 'ftp://ftp.ncbi.nih.gov/pub/COG/COG2014/data/cognames2003-2014.tab',
                    'func': self.format_cog_names,
                    'type': 'essential',
                    'formatted_file_name': 'COG.txt'},
                'fun2003-2014.tab': {
                    'url': 'ftp://ftp.ncbi.nih.gov/pub/COG/COG2014/data/fun2003-2014.tab',
                    'func': self.format_categories,
                    'type': 'essential',
                    'formatted_file_name': 'CATEGORIES.txt'},
                'prot2003-2014.fa.gz': {
                    'url': 'ftp://ftp.ncbi.nih.gov/pub/COG/COG2014/data/prot2003-2014.fa.gz',
                    'func': self.format_protein_db,
                    'type': 'database',
                    'formatted_file_name': 'IGNORE_THIS_AND_SEE_THE_FUNCTION'}
        }

        self.cogs_found_in_proteins_fasta = set([])
        self.cogs_found_in_cog_names_file = set([])


    def get_formatted_db_paths(self):
        formatted_db_paths = {}

        diamond_db_path = J(self.COG_data_dir, 'DB_DIAMOND')
        if os.path.exists(diamond_db_path):
            formatted_db_paths['diamond'] = J(diamond_db_path, 'COG')

        blast_db_path = J(self.COG_data_dir, 'DB_BLAST')
        if os.path.exists(blast_db_path):
            formatted_db_paths['blastp'] = J(blast_db_path, 'COG/COG.fa')

        return formatted_db_paths


    def get_essential_file_paths(self):
        if not os.path.exists(self.COG_data_dir):
            # the COG_data_dir is not there
            return None

        essential_files = {}
        for v in list(self.files.values()):
            if v['type'] == 'essential':
                essential_files[v['formatted_file_name']] = J(self.COG_data_dir, v['formatted_file_name'])

        # add the missing COG IDs file into the list:
        essential_files['MISSING_COG_IDs.cPickle'] = J(self.COG_data_dir, 'MISSING_COG_IDs.cPickle')

        for file_name in essential_files:
            if not os.path.exists(essential_files[file_name]):
                raise ConfigError("At least one essential formatted file that is necesary for COG operations is not where it should\
                                    be ('%s'). You should run COG setup, with the flag `--reset` if necessary, to make sure things\
                                    are in order." % essential_files[file_name])

        return essential_files


    def create(self):
        run.info('COG data dir', self.COG_data_dir)

        if not os.path.exists(self.COG_data_dir):
            try:
                os.mkdir(self.COG_data_dir)
                open(self.COG_data_dir_version, 'w').write(COG_DATA_VERSION)
            except Exception as e:
                raise ConfigError("So the COG data directory is not there, and anvi'o wants to create one. But it didn't\
                                    go that well. It could be due to permissions (which may require you to run this with sudo\
                                    or may need to ask your sys admin to do it for you since this is a one time operation), or\
                                    it could be due to something totally irrelevant. Here is the error message: '%s'" % e)

        filesnpaths.is_output_dir_writable(self.COG_data_dir)

        if self.reset:
            run.warning('This program will remove everything in the COG data directory, then download and reformat\
                         everything from scratch.')
            self.wait_for_the_user()

            # OK. reset the crap out of it.
            shutil.rmtree(self.COG_data_dir)
            os.mkdir(self.COG_data_dir)
            open(self.COG_data_dir_version, 'w').write(COG_DATA_VERSION)
        else:
            run.warning("This program will first check whether you have all the raw files, and then will attempt to\
                         regenerate everything that is necessary from them.")
            self.wait_for_the_user()

        if not os.path.exists(self.COG_data_dir_version) or open(self.COG_data_dir_version).read().strip() != COG_DATA_VERSION:
            raise ConfigError("The version of your COG data directory is different than what anvi'o hoping to see.\
                                It seems you need to (re)run anvi'o script to download and format COG data from NCBI.")

        # get raw files
        self.get_raw_data()

        # format raw files
        self.setup_raw_data()

        # identify missing COGs
        self.generate_missing_cog_ids_file()


    def generate_missing_cog_ids_file(self):
        missing_cog_ids = self.cogs_found_in_proteins_fasta.difference(self.cogs_found_in_cog_names_file)

        if len(missing_cog_ids):
            self.run.warning("%d of %d COG IDs that appear in the list of orthology domains file (which links protein IDs\
                              to COG names), are missing from the COG names file (which links COG IDs to function names and\
                              categoires). Because clearly even the files that are distributed together should not be expected to\
                              be fully compatible. Anvi'o thanks everyone for their contributions." % \
                                                        (len(missing_cog_ids), len(self.cogs_found_in_proteins_fasta)))

        dictio.write_serialized_object(missing_cog_ids, J(self.COG_data_dir, 'MISSING_COG_IDs.cPickle'))


    def format_p_id_to_cog_id_cPickle(self, input_file_path, output_file_path):
        progress.new('Formatting protein ids to COG ids file')
        progress.update('...')

        p_id_to_cog_id = {}

        for line in open(input_file_path, 'rU').readlines():
            fields = line.strip('\n').split(',')
            p_id = fields[0]
            COG = fields[6]

            self.cogs_found_in_proteins_fasta.add(COG)

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

        try:
            lines = open(input_file_path).readlines()
        except UnicodeDecodeError:
            lines = open(input_file_path, encoding='ISO-8859-1').readlines()

        for line in lines:
            if line.startswith('#'):
                continue
            COG, function, name = line.strip('\n').split('\t')
            self.cogs_found_in_cog_names_file.add(COG)

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

            output.write('\t'.join([category, description]) + '\n')

        progress.end()


    def format_protein_db(self, input_file_path, output_file_path):
        progress.new('Formatting raw files')
        progress.update('Decompressing protein sequences')

        # poor man's uncompress
        temp_fasta_path = filesnpaths.get_temp_file_path()
        with open(temp_fasta_path, 'wb') as f_out, gzip.open(input_file_path, 'rb') as f_in:
            f_out.write(f_in.read())

        progress.end()

        if utils.is_program_exists('diamond', dont_raise=True):
            output_dir = J(self.COG_data_dir, 'DB_DIAMOND')
            if os.path.exists(output_dir):
                shutil.rmtree(output_dir)

            os.mkdir(output_dir)

            output_db_path = J(output_dir, 'COG')
            log_file_path = J(output_dir, 'log.txt')

            self.run.info('Diamond log', log_file_path)

            diamond = Diamond(temp_fasta_path)
            diamond.num_threads = self.num_threads
            diamond.run.log_file_path = log_file_path
            diamond.makedb(output_db_path)
        else:
            self.run.warning("Diamond does not seem to be installed on this system, so anvi'o is not going to\
                              generate a search database for it. Remember this when/if things go South.")

        if utils.is_program_exists('makeblastdb', dont_raise=True) and utils.is_program_exists('blastp', dont_raise=True):
            output_dir = J(self.COG_data_dir, 'DB_BLAST')
            if os.path.exists(output_dir):
                shutil.rmtree(output_dir)

            os.mkdir(output_dir)

            output_db_path = J(output_dir, 'COG')
            log_file_path = J(output_dir, 'log.txt')

            self.run.info('BLAST log', log_file_path)

            blast = BLAST(temp_fasta_path)
            blast.run.log_file_path = log_file_path
            blast.num_threads = self.num_threads
            blast.makedb(os.path.join(output_db_path, 'COG.fa'))
        else:
            self.run.warning("BLAST tools do not seem to be installed on this system, so anvi'o is not going to\
                              generate a search database for them to be used. Keep this in mind for later.")

        os.remove(temp_fasta_path)


    def get_raw_data(self):
        if not os.path.exists(self.raw_NCBI_files_dir):
            os.mkdir(self.raw_NCBI_files_dir)
            open(self.COG_data_dir_version, 'w').write(COG_DATA_VERSION)

        for file_name in self.files:
            if not 'url' in self.files[file_name]:
                continue

            file_path = J(self.raw_NCBI_files_dir, file_name)
            if not os.path.exists(file_path):
                utils.download_file(self.files[file_name]['url'], file_path, progress=progress, run=run)


    def setup_raw_data(self):
        for file_name in self.files:
            file_path = J(self.raw_NCBI_files_dir, file_name)

            if not 'func' in self.files[file_name]:
                continue

            if not os.path.exists(file_path):
                raise ConfigError("Something is wrong :/ Raw files are not in place...")

            self.files[file_name]['func'](file_path, J(self.COG_data_dir, self.files[file_name]['formatted_file_name']))


    def wait_for_the_user(self):
        if self.just_do_it:
            return

        try:
            input("Press ENTER to continue, or press CTRL + C to cancel...\n")
        except:
            sys.exit()


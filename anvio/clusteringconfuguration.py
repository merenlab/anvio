# -*- coding: utf-8 -*-
# pylint: disable=line-too-long
"""To make sense of config files for mixed clustering"""

import os
import argparse
import configparser

import anvio
import anvio.db as db
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths

from anvio.utils import check_sample_id
from anvio.utils import store_array_as_TAB_delimited_file as store_array
from anvio.utils import store_dict_as_TAB_delimited_file
from anvio.utils import is_all_columns_present_in_TAB_delim_file as cols_present
from anvio.utils import get_vectors_from_TAB_delim_matrix as get_vectors
from anvio.errors import ConfigError
from anvio.tables.miscdata import TableForItemAdditionalData

run = terminal.Run()
progress = terminal.Progress()

__author__ = "Developers of anvi'o (see AUTHORS.txt)"
__copyright__ = "Copyleft 2015-2018, the Meren Lab (http://merenlab.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "A. Murat Eren"
__email__ = "a.murat.eren@gmail.com"
__status__ = "Development"


config_template = {
    'general': {
                'name': {'mandatory': False, 'test': lambda x: len(x.split()) == 1, 'required': 'a single-word string'},
                'output_file': {'mandatory': False, 'test': lambda x: filesnpaths.is_output_file_writable(x)},
                'num_components': {'mandatory': False, 'test': lambda x: RepresentsInt(x) and int(x) > 0 and int(x) <= 256,
                                   'required': "an integer value between 1 and 256"},
                'distance': {'mandatory': False, 'test': lambda x: len(x.split()) == 1, 'required': 'a single-word string'},
                'linkage': {'mandatory': False, 'test': lambda x: len(x.split()) == 1, 'required': 'a single-word string'},
                'seed': {'mandatory': False, 'test': lambda x: RepresentsInt(x), 'required': 'an integer'}
    },
    'matrix': {
                'table_form': {'mandatory': False, 'test': lambda x: x in ['dataframe', 'matrix'], 'required': 'matrix or dataframe'},
                'columns_to_use': {'mandatory': False, 'test': lambda x: len(x.strip().replace(' ', '').split(',')) > 0,
                            'required': 'one or more comma-separated column names'},
                'ratio': {'mandatory': False, 'test': lambda x: RepresentsInt(x) and int(x) > 0 and int(x) <= 256,
                          'required': "an integer value between 1 and 256."},
                'normalize': {'mandatory': False, 'test': lambda x: x in ['True', 'False'], 'required': 'True or False'},
                'log': {'mandatory': False, 'test': lambda x: x in ['True', 'False'], 'required': 'True or False'},
               },
}


def RepresentsInt(s):
    try:
        int(s)
        return True
    except ValueError:
        return False

def NameIsOK(n):
    try:
        check_sample_id(n)
    except ConfigError:
        return False
    return True


class ClusteringConfiguration:
    def __init__(self, config_file_path, input_directory=None, db_paths={}, row_ids_of_interest=[], r=run, p=progress):
        self.run = r
        self.progress = p

        self.input_directory = input_directory or os.path.abspath(os.getcwd())
        self.config_file_path = config_file_path

        # `row_ids_of_interest` gives opportunity to filter out irrelevant entries quickly
        # while vectors are being obtained from each matrix described in the config file.
        # to see why it is important in the context of anvi'o, see
        # https://github.com/meren/anvio/issues/100
        self.row_ids_of_interest = set(row_ids_of_interest)

        # these are the database files that may be referenced from within the config files
        # with !DATABASE.db::table notation. If a database entry has an exclamation mark,
        # it will be searched for in the db_paths dict to associate it with the relative
        # path that is only known to the client
        self.db_paths = db_paths

        # read the config
        filesnpaths.is_file_exists(self.config_file_path)
        config = configparser.ConfigParser()
        config.read(self.config_file_path)

        # this will keep the actual paths for each matrix:
        self.matrix_paths = {}
        self.set_default_paths(config)

        self.check_for_db_requests(config)

        # and sanity check.
        self.sanity_check(config)

        if self.get_option(config, 'general', 'output_file', str):
            self.output_file_name = self.get_option(config, 'general', 'output_file', str)
            self.output_file_path = os.path.join(self.input_directory, self.output_file_name)
        else:
            self.output_file_name = None
            self.output_file_path = None

        self.name = self.get_option(config, 'general', 'name', str) or filesnpaths.get_name_from_file_path(self.config_file_path)
        self.distance = self.get_option(config, 'general', 'distance', str)
        self.linkage = self.get_option(config, 'general', 'linkage', str)

        self.num_components = self.get_option(config, 'general', 'num_components', int)
        self.seed = self.get_option(config, 'general', 'seed', int)
        self.master = None

        self.matrices_dict = {}
        self.matrices = []
        for section in self.get_other_sections(config):
            alias, matrix = section.split()

            self.matrices.append(alias)

            m = {}
            columns_to_use = self.get_option(config, section, 'columns_to_use', str)
            table_form = self.get_option(config, section, 'table_form', str)
            m['alias'] = alias
            m['matrix'] = matrix
            m['table_form'] = table_form
            m['columns_to_use'] = [c.strip() for c in columns_to_use.split(',')] if columns_to_use else None
            m['ratio'] = self.get_option(config, section, 'ratio', int)
            m['path'] = self.matrix_paths[alias]
            m['normalize'] = False if self.get_option(config, section, 'normalize', str) == 'False' else True
            m['log'] = True if self.get_option(config, section, 'log', str) == 'True' else False
            # next two variables are necessary to follow the order of vectors
            m['id_to_sample'], m['sample_to_id'], m['cols'], m['vectors'] = get_vectors(m['path'], m['columns_to_use'], self.row_ids_of_interest)
            self.matrices_dict[alias] = m

        # make sure all matrices have identical rows:
        if len(set([list(m['id_to_sample'].values()).__str__() for m in list(self.matrices_dict.values())])) > 1:
            master_rows, master_matrix = sorted([(len(self.matrices_dict[m]['id_to_sample']), list(self.matrices_dict[m]['id_to_sample'].values()), m)\
                                                            for m in self.matrices_dict])[0][1:]
            self.master = master_matrix
            self.master_rows = master_rows
            # the smallest matrix is 'master_matrix', and the rows it has is master_rows. so every other matrix
            # must match that, or we will throw a tantrum.
            for matrix in [m for m in self.matrices if m != master_matrix]:
                m = self.matrices_dict[matrix]

                # get reduced set of vectors from rows that match `master_rows`:
                m['id_to_sample'], m['sample_to_id'], m['cols'], m['vectors'] = get_vectors(m['path'], m['columns_to_use'], master_rows)

                if len(m['vectors']) != len(master_rows):
                    raise ConfigError('The content of rows differed between input matrices. So I tried to\
                                        match all other matrices to the matrix with the smallest number of\
                                        rows (which was "%s"). However, not all other matrices contained\
                                        the small set of rows.' % (master_matrix))
        else:
            self.master_rows = sorted(self.matrices_dict[self.matrices[0]]['sample_to_id'].keys())

        self.num_matrices = len(self.matrices)
        self.multiple_matrices = self.num_matrices > 1

        # following section will be irrelevant for a while (this is tied to working on
        # the clustering.order_contigs_experimental() function:
        #if not self.multiple_matrices:
        #    # there is only one matrix, we don't expect to see a ratio.
        #    if self.matrices_dict[self.matrices[0]]['ratio']:
        #        raise ConfigError, 'There is only one matrix declared in the config file. Which renders the\
        #                            "ratio" variables irrelevant. Please make sure it is not set in the\
        #                             config file.'
        #    if self.num_components:
        #        raise ConfigError, 'There is only one matrix declared in the config file. In this there\
        #                            will be no scaling step. Therefore the "num_components" variable will\
        #                            not be used. Please make sure it is not set under the general section.'
        #else:
        #    # if there are multiple matrices, it means this config is going to be used to
        #    # scale and mix the data, therefore it is mandatory to have num_components
        #    # defined.
        #    if self.multiple_matrices and not self.get_option(config, 'general', 'num_components', int):
        #        raise ConfigError, 'When multiple matrices are defined, it is mandatory to define the number of\
        #                            components under the general section ("num_components").'


    def print_summary(self, r=run):
        r.info_single('Summary of the config file:', mc='green', nl_before=2)
        r.warning('', header='General')
        r.info('Input directory', self.input_directory)
        r.info('Number of components', self.num_components)
        r.info('Seed', self.seed)
        r.info('Output file', self.output_file_name)
        for alias in self.matrices:
            m = self.matrices_dict[alias]
            r.warning('', header='%s (%s) %s' % (alias, m['alias'], '[MASTER]' if alias == self.master else ''))
            r.info('Ratio', m['ratio'])
            r.info('Normalize', m['normalize'])
            r.info('Log', m['log'])
            r.info('Columns to use', m['columns_to_use'] or 'ALL')
            r.info('Num Columns', len(m['cols']))
            r.info('Num Rows', len(m['vectors']))


    def get_option(self, config, section, option, cast):
        try:
            return cast(config.get(section, option).strip())
        except configparser.NoOptionError:
            return None


    def get_other_sections(self, config):
        return [s for s in config.sections() if s != 'general']


    def check_section(self, config, section, template_class):
        """`section` is the actual section name in the config file, `template_class`
            corresponds to what type of section it is..."""
        for option, value in config.items(section):
            if option not in list(config_template[template_class].keys()):
                raise ConfigError('Unknown option, "%s", under section "%s".' % (option, section))
            if 'test' in config_template[template_class][option] and not config_template[template_class][option]['test'](value):
                if 'required' in config_template[template_class][option]:
                    r = config_template[template_class][option]['required']
                    raise ConfigError('Unexpected value ("%s") for option "%s", under section "%s".\
                                        What is expected is %s.' % (value, option, section, r))
                else:
                    raise ConfigError('Unexpected value ("%s") for option "%s", under section "%s".' % (value, option, section))

        for option in config_template[template_class]:
            if 'mandatory' in config_template[template_class][option] and config_template[template_class][option]['mandatory'] and not config.has_option(section, option):
                raise ConfigError('Missing mandatory option for section "%s": %s' % (section, option))


    def set_default_paths(self, config):
        sections = self.get_other_sections(config)

        # set default paths:
        for section in sections:
            try:
                alias, matrix = section.split()
            except:
                raise ConfigError('Each section must have "alias" and "matrix" fields separated by\
                                    a white space.')
            self.matrix_paths[alias] = os.path.join(self.input_directory, matrix)


    def check_for_db_requests(self, config):
        sections = self.get_other_sections(config)
        # look for requests from the database, create temporary tab delimited files:
        for section in sections:
            alias, matrix = section.split()
            if matrix.find('::') > -1:
                if matrix.startswith('!'):
                    database, table = matrix.split('::')
                    database = database[1:]

                    if database not in self.db_paths:
                        raise ConfigError('anvio could not recover the actual path of the database\
                                            (!%s) referenced in the config file, because the database\
                                            paths variable sent from the client does not have an entry\
                                            for it :( There are two options. One is to get a db_paths\
                                            dictionary sent to this class that contains a key for %s\
                                            with the full path to the dataase as a value. Or the table\
                                            "%s" can be exported to a TAB-delimited matrix and declared in\
                                            the config file. If you are experimenting and stuck here, please\
                                            see the documentation or send an e-mail to the developers.'\
                                                                                % (database, database, table))
                    database_path = self.db_paths[database]
                else:
                    database, table = matrix.split('::')
                    database_path = os.path.abspath(self.db_paths[database]) if database in self.db_paths else os.path.abspath(database)

                    # if its not there, let's try one more thing
                    if not os.path.exists(database_path):
                        database_path = os.path.abspath(os.path.join(self.input_directory, database))

                if not os.path.exists(database_path):
                    raise ConfigError("The database you requested (%s) is not where it was supposed to be ('%s') :/" % (database, database_path))

                dbc = db.DB(database_path, None, ignore_version=True)

                if not table in dbc.get_table_names():
                    raise ConfigError('The table you requested (%s) does not seem to be in %s :/' % (table, database))

                # here we know we are working with a database table that we have access to. however, in anvi'o database
                # tables in two forms: dataframe form, and matrix form. in dataframe form, we have key/value pairs rather
                # than MxN matrices where each N is a column for an attribute. while the latter is easier to export as a
                # matrix the clustering module can work with, the former requires extra attention. so here we need to first
                # figure out whether which form the table is in. why this even became necessary? taking a look at this issue
                # may help: https://github.com/merenlab/anvio/issues/662
                table_form = None
                if config.has_option(section, 'table_form'):
                    table_form = config.get(section, 'table_form')

                table_rows = dbc.get_all_rows_from_table(table)

                if self.row_ids_of_interest:
                    if table_form == 'dataframe':
                        raise ConfigError("Oops .. anvi'o does not know how to deal with specific row ids of interest when a table\
                                           refernced from a clustering recipe is in dataframe form :(")
                    table_rows = [r for r in table_rows if r[0] in self.row_ids_of_interest]

                if not len(table_rows):
                    raise ConfigError("It seems the table '%s' in the database it was requested from is empty. This\
                                        is not good. Here is the section that is not working for you: '%s' :/" \
                                                                % (table, section))

                tmp_file_path = filesnpaths.get_temp_file_path()

                # time to differentially store table contents.
                if table_form == 'dataframe':
                    args = argparse.Namespace(pan_or_profile_db=database_path, table_name=table)
                    table = TableForItemAdditionalData(args)
                    table_keys_list, table_data_dict = table.get()
                    store_dict_as_TAB_delimited_file(table_data_dict, tmp_file_path)
                else:
                    table_structure = dbc.get_table_structure(table)
                    columns_to_exclude = [c for c in ['entry_id', 'sample_id'] if c in table_structure]
                    store_array(table_rows, tmp_file_path, table_structure, exclude_columns=columns_to_exclude)

                self.matrix_paths[alias] = tmp_file_path


    def sanity_check(self, config):
        filesnpaths.is_file_exists(self.input_directory)

        if 'general' not in config.sections():
            raise ConfigError("[general] section is mandatory.")

        if len(config.sections()) < 2:
            raise ConfigError("Config file must contain at least one matrix section.")

        self.check_section(config, 'general', 'general')

        sections = self.get_other_sections(config)

        for section in sections:
            alias, matrix = section.split()

            if not os.path.exists(self.matrix_paths[alias]):
                raise ConfigError('The matrix file "%s" you mentioned in %s has not been found in the\
                                    input directory :/' % (matrix, os.path.basename(self.config_file_path)))

            self.check_section(config, section, 'matrix')

            # it is not very elegant to do this here, but carrying this test in the template was going to
            # cause a lot of uncalled for complexity (or reporting was going to be very vague):
            columns_to_use_str = self.get_option(config, section, 'columns_to_use', str)
            columns_to_use = [c.strip() for c in columns_to_use_str.split(',')] if columns_to_use_str else None
            if columns_to_use and not cols_present(columns_to_use, self.matrix_paths[alias]):
                raise ConfigError('One or more of the columns declared for "%s" in the config file\
                                    seem(s) to be missing in the matrix :/' % (matrix))


        # 'ratio' must be defined either for all, or for none of the matrices
        with_ratio = len([True for section in sections if self.get_option(config, section, 'ratio', int)])
        if with_ratio and with_ratio != len(sections):
            raise ConfigError('Ratio value must be defined either for all, or none of the matrices. In your\
                                configuration only %d of %d matrices have ratio values defined. Either remove\
                                all, or complete the remaining one%s.' % (with_ratio, len(sections),
                                                                          's' if (len(sections) - with_ratio) > 1 else ''))



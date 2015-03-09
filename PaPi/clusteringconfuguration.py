# -*- coding: utf-8 -*-

# Copyright (C) 2014, A. Murat Eren
#
# This program is free software; you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free
# Software Foundation; either version 2 of the License, or (at your option)
# any later version.
#
# Please read the COPYING file.


import os
import sys
import ConfigParser

import PaPi.db as db
import PaPi.terminal as terminal
import PaPi.filesnpaths as filesnpaths

from PaPi.utils import ConfigError
from PaPi.utils import check_sample_id
from PaPi.constants import allowed_chars
from PaPi.utils import store_array_as_TAB_delimited_file as store_array
from PaPi.utils import is_all_columns_present_in_TAB_delim_file as cols_present
from PaPi.utils import get_vectors_from_TAB_delim_matrix as get_vectors


config_template = {
    'general': {
                'output_file'    : {'mandatory': False, 'test': lambda x: filesnpaths.is_output_file_writable(x)},
                'num_components': {'mandatory': False, 'test': lambda x: RepresentsInt(x) and int(x) > 0 and int(x) <= 256,
                                   'required': "an integer value between 1 and 256"},
                'seed': {'mandatory': False, 'test': lambda x: RepresentsInt(x), 'required': 'an integer'}
    },
    'matrix': {
                'columns_to_use': {'mandatory': False, 'test': lambda x: len(x.strip().replace(' ','').split(',')) > 1,
                            'required': 'more than one, comma-separated column names'},
                'ratio': {'mandatory': False, 'test': lambda x: RepresentsInt(x) and int(x) > 0 and int(x) <= 256,
                          'required': "an integer value between 1 and 256."},
                'alias': {'mandatory': True, 'test': lambda x: NameIsOK(x),
                         'required': 'a single word alias composed of these characters alone: "%s"' % allowed_chars},
                'normalize': {'mandatory': False, 'test': lambda x: x in ['True', 'False'], 'required': 'True or False'},
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
    def __init__(self, config_file_path, input_directory = None, version = None):
        self.input_directory = input_directory or os.getcwd()
        self.config_file_path = config_file_path
        self.version = version

        # read the config
        filesnpaths.is_file_exists(self.config_file_path)
        config = ConfigParser.ConfigParser()
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

        self.num_components = self.get_option(config, 'general', 'num_components', int)
        self.seed = self.get_option(config, 'general', 'seed', int)
        self.master = None

        self.matrices_dict = {}
        self.matrices = []
        for matrix in self.get_other_sections(config):
            self.matrices.append(matrix)
            m = {}
            columns_to_use = self.get_option(config, matrix, 'columns_to_use', str)
            m['name'] = matrix
            m['columns_to_use'] = [c.strip() for c in columns_to_use.split(',')] if columns_to_use else None
            m['ratio'] = self.get_option(config, matrix, 'ratio', int)
            m['alias'] = self.get_option(config, matrix, 'alias', str)
            m['path'] = self.matrix_paths[matrix]
            m['normalize'] = False if self.get_option(config, matrix, 'normalize', str) == 'False' else True 
            # next two variables are necessary to follow the order of vectors
            m['id_to_sample'], m['cols'], m['vectors'] = get_vectors(m['path'], m['columns_to_use'])
            m['sample_to_id'] = dict([(v, k) for k, v in m['id_to_sample'].iteritems()])
            self.matrices_dict[matrix] = m

        # make sure all matrices have equal number of entries:
        if len(set([m['id_to_sample'].values().__str__() for m in self.matrices_dict.values()])) > 1:
            master_rows, master_matrix = sorted([(len(self.matrices_dict[m]['id_to_sample']), self.matrices_dict[m]['id_to_sample'].values(), m)\
                                                            for m in self.matrices_dict])[0][1:]
            self.master = master_matrix
            self.master_rows = master_rows
            # the smallest matrix is 'master_matrix', and the rows it has is master_rows. so every other matrix
            # must match that, or we will throw a tantrum.
            for matrix in [m for m in self.matrices if m != master_matrix]:
                m = self.matrices_dict[matrix]
                m['id_to_sample'], m['cols'], m['vectors'] = get_vectors(m['path'], m['columns_to_use'], master_rows)
                if len(m['vectors']) != len(master_rows):
                    raise ConfigError, 'The content of rows differed between input matrices. So I tried to\
                                        match all other matrices to the matrix with the smallest number of\
                                        rows (which was "%s"). However, not all other matrices contained\
                                        the small set of rows.' % (master_matrix)
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


    def print_summary(self, r):
        r.warning('', header = 'General')
        r.info('Input directory', self.input_directory)
        r.info('Number of components', self.num_components)
        r.info('Seed', self.seed)
        r.info('Output file', self.output_file_name)
        for matrix in self.matrices:
            m = self.matrices_dict[matrix]
            r.warning('', header = '%s (%s) %s' % (matrix, m['alias'], '[MASTER]' if matrix == self.master else ''))
            r.info('Ratio', m['ratio'])
            r.info('Normalize', m['normalize'])
            r.info('Columns to use', m['columns_to_use'] or 'ALL')
            r.info('Num Columns', len(m['cols']))
            r.info('Num Rows', len(m['vectors']))


    def get_option(self, config, section, option, cast):
        try:
            return cast(config.get(section, option).strip())
        except ConfigParser.NoOptionError:
            return None


    def get_other_sections(self, config):
        return [s for s in config.sections() if s != 'general']


    def check_section(self, config, section, template_class):
        """`section` is the actual section name in the config file, `template_class`
            corresponds to what type of section it is..."""
        for option, value in config.items(section):
            if option not in config_template[template_class].keys():
                raise ConfigError, 'Unknown option, "%s", under section "%s".' % (option, section)
            if config_template[template_class][option].has_key('test') and not config_template[template_class][option]['test'](value):
                if config_template[template_class][option].has_key('required'):
                    r = config_template[template_class][option]['required']
                    raise ConfigError, 'Unexpected value ("%s") for option "%s", under section "%s".\
                                        What is expected is %s.' % (value, option, section, r)
                else:
                    raise ConfigError, 'Unexpected value ("%s") for option "%s", under section "%s".' % (value, option, section)

        for option in config_template[template_class]:
            if config_template[template_class][option].has_key('mandatory') and config_template[template_class][option]['mandatory'] and not config.has_option(section, option):
                raise ConfigError, 'Missing mandatory option for section "%s": %s' % (section, option)


    def set_default_paths(self, config):
        matrices = self.get_other_sections(config)

        # set default paths:
        for matrix in matrices:
            self.matrix_paths[matrix] = os.path.join(self.input_directory, matrix)


    def check_for_db_requests(self, config):
        matrices = self.get_other_sections(config)
        # look for requests from the database, create temporary tab delimited files:
        for matrix in matrices:
            if matrix.find('::') > -1:
                database, table = matrix.split('::')
                database_path = os.path.join(self.input_directory, database)
                if not os.path.exists(database_path):
                    raise ConfigError, 'The database you requested (%s) is not in the input directory :/' % database

                dbc = db.DB(database_path, self.version)

                if not table in dbc.get_table_names():
                    raise ConfigError, 'The table you requested (%s) does not seem to be in %s :/' % (table, database)

                table_rows = dbc.get_all_rows_from_table(table)
                tmp_file_path = filesnpaths.get_temp_file_path()
                table_structure = dbc.get_table_structure(table)
                columns_to_exclude = [c for c in ['entry_id', 'sample_id'] if c in table_structure]
                store_array(table_rows, tmp_file_path, table_structure, exclude_columns = columns_to_exclude)
                self.matrix_paths[matrix] = tmp_file_path


    def sanity_check(self, config):
        filesnpaths.is_file_exists(self.input_directory)

        if 'general' not in config.sections():
            raise ConfigError, "[general] section is mandatory."

        if len(config.sections()) < 2:
            raise ConfigError, "Config file must contain at least one matrix section."

        self.check_section(config, 'general', 'general')

        matrices = self.get_other_sections(config)


        for matrix in matrices:
            if not os.path.exists(self.matrix_paths[matrix]):
                raise ConfigError, 'The matrix file "%s" you mentioned in %s has not been found in the\
                                    input directory :/' % (matrix, os.path.basename(self.config_file_path))


            self.check_section(config, matrix, 'matrix')

            # it is not very elegant to do this here, but carrying this test in the template was going to
            # cause a lot of uncalled for complexity (or reporting was going to be very vague):
            columns_to_use_str = self.get_option(config, matrix, 'columns_to_use', str)
            columns_to_use = [c.strip() for c in columns_to_use_str.split(',')] if columns_to_use_str else None
            if columns_to_use and not cols_present(columns_to_use, self.matrix_paths[matrix]):
                raise ConfigError, 'One or more of the columns declared for "%s" in the config file\
                                    seem(s) to be missing in the matrix :/' % (matrix)

        # 'ratio' must be defined either for all, or for none of the matrices
        with_ratio = len([True for matrix in matrices if self.get_option(config, matrix, 'ratio', int)])
        if with_ratio and with_ratio != len(matrices):
            raise ConfigError, 'Ratio value must be defined either for all, or none of the matrices. In your\
                                configuration only %d of %d matrices have ratio values defined. Either remove\
                                all, or complete the remaining one%s.' % (with_ratio, len(matrices),
                                                                          's' if (len(matrices) - with_ratio) > 1 else '')



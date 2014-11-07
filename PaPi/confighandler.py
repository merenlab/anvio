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

import PaPi.terminal as terminal
import PaPi.filesnpaths as filesnpaths

from PaPi.utils import ConfigError
from PaPi.utils import check_sample_id
from PaPi.constants import allowed_chars
from PaPi.utils import is_all_columns_present_in_TAB_delim_file as cols_present


config_template = {
    'general': {
                'output_file'    : {'mandatory': True, 'test': lambda x: filesnpaths.is_output_file_writable(x)},
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
                         'required': 'a single word alias composed of these characters alone: "%s"' % allowed_chars}
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
    def __init__(self, config_file_path, input_directory = None):
        self.input_directory = input_directory or os.getcwd()

        # read the config
        filesnpaths.is_file_exists(config_file_path)
        config = ConfigParser.ConfigParser()
        config.read(config_file_path)

        # and sanity check.
        self.sanity_check(config)

        self.output_file_name = self.get_option(config, 'general', 'output_file', str)
        self.output_file_path = os.path.join(self.input_directory, self.output_file_name)
        self.num_components = self.get_option(config, 'general', 'num_components', int)
        self.seed = self.get_option(config, 'general', 'seed', int)

        self.matrices_dict = {}
        self.matrices = []
        for matrix in self.get_other_sections(config):
            self.matrices.append(matrix)
            m = {}
            columns_to_use = self.get_option(config, matrix, 'columns_to_use', str)
            m['columns_to_use'] = [c.strip() for c in columns_to_use.split(',')] if columns_to_use else None
            m['ratio'] = self.get_option(config, matrix, 'ratio', int)
            m['alias'] = self.get_option(config, matrix, 'alias', str)
            m['path'] = os.path.join(self.input_directory, matrix)
            self.matrices_dict[matrix] = m

        self.num_matrices = len(self.matrices)
        self.multiple_matrices = self.num_matrices > 1

        if not self.multiple_matrices:
            # there is only one matrix, we don't expect to see a ratio.
            if self.matrices_dict[self.matrices[0]]['ratio']:
                raise ConfigError, 'There is only one matrix declared in the config file. Which renders the\
                                    "ratio" variables irrelevant. Please make sure it is not set in the\
                                     config file.'
            if self.num_components:
                raise ConfigError, 'There is only one matrix declared in the config file. In this there\
                                    will be no scaling step. Therefore the "num_components" variable will\
                                    not be used. Please make sure it is not set under the general section.'
        else:
            # if there are multiple matrices, it means this config is going to be used to
            # scale and mix the data, therefore it is mandatory to have num_components
            # defined.
            if self.multiple_matrices and not self.get_option(config, 'general', 'num_components', int):
                raise ConfigError, 'When multiple matrices are defined, it is mandatory to define the number of\
                                    components under the general section ("num_components").'


    def print_summary(self):
        r = terminal.Run(width=35)
        r.info('General', '', header=True)
        r.info('Input directory', self.input_directory)
        r.info('Number of components', self.num_components)
        r.info('Seed', self.seed)
        r.info('Output file', self.output_file_name)
        for matrix in self.matrices:
            r.info('%s (%s)' % (matrix, self.matrices_dict[matrix]['alias']), '', header=True)
            r.info('Columns to use', self.matrices_dict[matrix]['columns_to_use'])
            r.info('Ratio', self.matrices_dict[matrix]['ratio'])


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


    def sanity_check(self, config):
        filesnpaths.is_file_exists(self.input_directory)

        if 'general' not in config.sections():
            raise ConfigError, "[general] section is mandatory."

        if len(config.sections()) < 2:
            raise ConfigError, "Config file must contain at least one matrix sextion."

        self.check_section(config, 'general', 'general')

        matrices = self.get_other_sections(config)

        for matrix in matrices:
            if not os.path.exists(os.path.join(self.input_directory, matrix)):
                raise ConfigError, 'The matrix file "%s" you mentioned in the config file is not in the\
                                    input directory (if you have not specify an input directory, it is\
                                    assumed to be the "current working directory". If you have\
                                    not specified one, please specify the correct input directory' % (matrix)

            self.check_section(config, matrix, 'matrix')

            # it is not very elegant to do this here, but carrying this test in the template was going to
            # cause a lot of uncalled for complexity (or reporting was going to be very vague):
            columns_to_use_str = self.get_option(config, matrix, 'columns_to_use', str)
            columns_to_use = [c.strip() for c in columns_to_use_str.split(',')] if columns_to_use_str else None
            if columns_to_use and not cols_present(columns_to_use, os.path.join(self.input_directory, matrix)):
                raise ConfigError, 'One or more of the columns declared for "%s" in the config file\
                                    seem(s) to be missing in the matrix :/' % (matrix)

        # 'ratio' must be defined either for all, or for none of the matrices
        with_ratio = len([True for matrix in matrices if self.get_option(config, matrix, 'ratio', int)])
        if with_ratio and with_ratio != len(matrices):
            raise ConfigError, 'Ratio value must be defined either for all, or none of the matrices. In your\
                                configuration only %d of %d matrices have ratio values defined. Either remove\
                                all, or complete the remaining one%s.' % (with_ratio, len(matrices),
                                                                          's' if (len(matrices) - with_ratio) > 1 else '')



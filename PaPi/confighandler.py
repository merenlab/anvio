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

import PaPi.filesnpaths as filesnpaths
import PaPi.terminal as terminal
from PaPi.utils import ConfigError as ConfigError

config_template = {
    'general': {
                'output_file'    : {'mandatory': True, 'test': lambda x: filesnpaths.is_output_file_writable(x)},
                'num_components': {'mandatory': False, 'test': lambda x: RepresentsInt(x) and int(x) > 0 and int(x) <= 256,
                                   'required': "an integer value between 1 and 256"},
                'seed': {'mandatory': False, 'test': lambda x: RepresentsInt(x), 'required': 'an integer'}
    },
    'matrix': {
                'columns': {'mandatory': False, 'test': lambda x: len(x.strip().replace(' ','').split(',')) > 1,
                            'required': 'more than one, comma-separated sample names'},
                'ratio': {'mandatory': False, 'test': lambda x: RepresentsInt(x) and int(x) > 0 and int(x) <= 256,
                          'required': "an integer value between 1 and 256."},
               },
}


def RepresentsInt(s):
    try: 
        int(s)
        return True
    except ValueError:
        return False


class RunConfiguration:
    def __init__(self, config_file_path, input_directory = None):
        self.input_directory = input_directory or os.getcwd()

        # read the config
        filesnpaths.is_file_exists(config_file_path)
        config = ConfigParser.ConfigParser()
        config.read(config_file_path)

        # and sanity check.
        self.sanity_check(config)

        self.output_file = self.get_option(config, 'general', 'output_file', str)
        self.num_components = self.get_option(config, 'general', 'num_components', int)
        self.seed = self.get_option(config, 'general', 'seed', int)

        self.matrices = {}
        for matrix in self.get_other_sections(config):
            m = {}
            columns = self.get_option(config, matrix, 'columns', str)
            m['columns'] = [c.strip() for c in columns.split(',')] if columns else None
            m['ratio'] = self.get_option(config, matrix, 'ratio', int)
            self.matrices[matrix] = m


    def print_summary(self):
        r = terminal.Run()
        r.info('General', '', header=True)
        r.info('Input directory', self.input_directory)
        r.info('Number of components', self.num_components)
        r.info('Seed', self.seed)
        r.info('Output file', self.output_file)
        for matrix in self.matrices:
            r.info(matrix, '', header=True)
            r.info('Columns', self.matrices[matrix]['columns'])
            r.info('Ratio', self.matrices[matrix]['ratio'])


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
                raise ConfigError, 'Unknown option under "%s" section: "%s"' % (section, option)
            if config_template[template_class][option].has_key('test') and not config_template[template_class][option]['test'](value):
                if config_template[template_class][option].has_key('required'):
                    r = config_template[template_class][option]['required']
                    raise ConfigError, 'Unexpected value for "%s" section "%s": "%s".\
                                        What is expected is %s.' % (option, section, value, r)
                else:
                    raise ConfigError, 'Unexpected value for "%s" section "%s": "%s".' % (option, section, value)

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
        for section in self.get_other_sections(config):
            if not os.path.exists(os.path.join(self.input_directory, section)):
                raise ConfigError, 'The matrix file "%s" you mentioned in the config file is not in the\
                                    input directory (if you have not specify an input directory, it is\
                                    assumed to be the "current working directory". If you have\
                                    not specified one, please specify the correct input directory' % (section)
            self.check_section(config, section, 'matrix')


# coding: utf-8
# pylint: disable=line-too-long
"""This is the main script for anvi-run-batch program."""

import os
import re
import sys
import json
import yaml
import anvio
import subprocess
import anvio.utils as utils
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths
from anvio.argparse import ArgumentParser
from anvio.errors import ConfigError, FilesNPathsError

__author__ = "Developers of anvi'o (see AUTHORS.txt)"
__copyright__ = "Copyleft 2015-2018, the Meren Lab (http://merenlab.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['metehan']
__requires__ = ['contigs-db', 'profile-db']
__description__ = ("This is helper script for anvi-run-batch program. It is used to run the commands in the yaml file.")

progress = terminal.Progress(verbose=False)
run = terminal.Run(verbose=False)
pp = terminal.pretty_print

class AnvioBatchWork():
    """ This class is used to do some magic after getting yaml file. """

    def __init__(self, args, skip_sanity_check=False, run=terminal.Run(), progress=terminal.Progress()):
        self.args = args
        self.run = run
        
        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        self.yaml_file_path = A('yaml')
        self.setup_arg = A('setup')
        self.given_work_dir = A('work_directory')
        self.contigs_db_path = A('contigs_db')
        self.profile_db_path = A('profile_db')
        self.pan_db_path = A('pan_db')

        if not skip_sanity_check:
            self.sanity_check()

        # Open and parse the file
        self.yaml_file = utils.get_yaml_as_dict(self.yaml_file_path)

        # to be filled later if necessary
        self.contigs_basic_info = {}


    def sanity_check(self):
        if not self.yaml_file_path:
            raise ConfigError("You must provide a YAML file path.")
        
        if not self.yaml_file_path[-4:] in ['yaml', '.yml']:
            raise ConfigError("You must provide a file with the file extension 'yaml' or 'yml.")
        filesnpaths.is_anvio_batch_yaml(self.yaml_file_path)

        if self.contigs_db_path or self.profile_db_path or self.pan_db_path:
            raise ConfigError("You should remove all the db_paths (CONTIGS.db, PROFILE.db...) to run this command. :/")


    def execute(self):
        """This function is used to run the commands in the yaml file."""

        self.run.info('Project Title', self.yaml_file.get('project').get('title'), mc='green', nl_before= 1)
        self.run.info('Project Description', self.yaml_file.get('project').get('description'), mc='green')
        self.run.info('Author Name', self.yaml_file.get('author').get('name'), mc='green')
        self.run.info('Author Email', self.yaml_file.get('author').get('email'), mc='green')
        self.run.info('Author Affiliation', self.yaml_file.get('author').get('affiliation'), mc='green')
        self.run.info('Author Web', self.yaml_file.get('author').get('web'), mc='green')

        #Trimming working directory.
        work_dir = re.sub('^[a-z0-9](?!.*?[^\na-z0-9]{2}).*?[a-z0-9]$', '', self.yaml_file.get('work_directory'))
        
        setup = self.yaml_file.get('setup')
        setup_command_counter = 0
        cwd = os.getcwd()
        unix_commands = ['cd','rm', 'ls', 'll', 'pwd', 'cat', 'mkdir', 'cp', 'mv', 'rmdir']

        if not work_dir:
            raise ConfigError('You must give your Working Directory!!')

        if self.setup_arg:
            while setup_command_counter < len(setup):
                if 'cd' in setup[setup_command_counter]:
                    self.run.warning('You should remove "cd" commands in your yaml file')
                    sys.exit(-1)
                else:
                    subprocess.run(str(setup[setup_command_counter]), shell=True)
                    setup_command_counter += 1
                    
        running_command = self.yaml_file.get('run')
        run_command_counter = 0

        #Change directory anyway!
        cwd_file = cwd + '/' + work_dir
        os.chdir(cwd_file)

        # WE ALWAYS RUN MIGRATION EVEN USER GIVE IN SETUP FILE. 2 is better than 1
        subprocess.call('anvi-migrate --migrate-dbs-safely --migrate-safely *.db', shell=True)

        while run_command_counter < len(running_command):
            subprocess.call(running_command[run_command_counter].get('command'), shell=True)
            run_command_counter += 1            




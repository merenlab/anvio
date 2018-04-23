# coding: utf-8
"""Interface to PyANI."""

import os

import anvio
import anvio.utils as utils
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths

from anvio.errors import ConfigError


__author__ = "Developers of anvi'o (see AUTHORS.txt)"
__copyright__ = "Copyleft 2015-2018, the Meren Lab (http://merenlab.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "Özcan Esen"
__email__ = "ozcanesen@gmail.com"


class PyANI:
    def __init__(self, args={}, run=terminal.Run(), progress=terminal.Progress()):
        self.run = run
        self.progress = progress
        self.program_name = 'average_nucleotide_identity.py'
        utils.is_program_exists(self.program_name)

        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        self.num_threads = A('num_threads') or 1
        self.method = A('method') or 'ANIb'
        self.log_file_path = os.path.abspath(A('log_file') or filesnpaths.get_temp_file_path())

        self.run.warning("Anvi'o will use 'PyANI' by Pritchard et al. (DOI: 10.1039/C5AY02550H) to compute ANI. If you publish your findings, \
                            please do not forget to properly credit their work.", lc='green', header="CITATION")

        self.run.info('[PyANI] Num threads to use', self.num_threads)
        self.run.info('[PyANI] Alignment method', self.method)
        self.run.info('[PyANI] Log file path', self.log_file_path)


    def run_command(self, input_path):
        # backup the old working directory before changing the directory
        old_wd = os.getcwd()
        os.chdir(input_path)


        full_command = [self.program_name,
                        '--outdir', 'output',
                        '--indir', input_path,
                        '--method', self.method,
                        '--workers', self.num_threads]

        self.progress.new('PyANI')
        self.progress.update('Running ...')
        exit_code = utils.run_command(full_command, self.log_file_path)
        self.progress.end()

        if int(exit_code):
            raise ConfigError("PyANI returned with non-zero exit code, there may be some errors. \
                              please check the log file for details.")

        with open(os.path.join(input_path, 'output', self.method + '_percentage_identity.tab'), 'r') as f:
            percent_identity = f.read()

        # restore old working directory
        os.chdir(old_wd)

        return percent_identity


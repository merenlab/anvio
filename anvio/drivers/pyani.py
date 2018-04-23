# coding: utf-8
"""Interface to PyANI."""

import os
import io
import sys

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


run = terminal.Run()
progress = terminal.Progress()
pp = terminal.pretty_print


class PyANI:
    def __init__(self, run=run):
        self.run = run
        self.progress = progress
        self.program_name = 'average_nucleotide_identity.py'
        utils.is_program_exists(self.program_name)

        self.run.warning("Anvi'o will use 'PyANI' by Pritchard et al. (DOI: 10.1039/C5AY02550H) to compute ANI. If you publish your findings, \
                            please do not forget to properly credit their work.", lc='green', header="CITATION")

    def run_command(self, input_path, method='ANIb'):
        # backup the old working directory before changing the directory
        old_wd = os.getcwd()
        os.chdir(input_path)

        log_file_path = filesnpaths.get_temp_file_path()
        self.run.info('Log file path', log_file_path)

        full_command = [self.program_name, '--outdir', 'output', '--indir', input_path, '-g', '-m', method]
        exit_code = utils.run_command(full_command, log_file_path)

        if exit_code != 0:
            self.run.warning("PyANI returned with non-zero exit code, there may be some errors. \
                              please check the log file for details.")

        J = lambda name: os.path.join(input_path, 'output', method + name)

        percent_identity = utils.get_TAB_delimited_file_as_dictionary(J('_percentage_identity.tab'))

        # restore old working directory
        os.chdir(old_wd)

        return percent_identity


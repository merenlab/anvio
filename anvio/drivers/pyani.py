# coding: utf-8
"""Interface to PyANI."""

import os
import io
from subprocess import Popen, PIPE

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

    def run_command(self, input_path, method='ANIb'):
        output_path = filesnpaths.get_temp_directory_path()

        full_command = [self.program_name, '-i', input_path, '-o', output_path, '-g', '-m', method]
        program = Popen(full_command, stdout=PIPE, stdin=PIPE, stderr=PIPE)
        sStdout, sStdErr = program.communicate()

        with open(os.path.join(output_path, method + '_percentage_identity.tab'), 'r') as f:
            percent_identity = f.read()

        return percent_identity


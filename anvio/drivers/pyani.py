# coding: utf-8
"""Interface to PyANI."""

import os
import io
import sys
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

        self.run.warning("Anvi'o will use 'PyANI' by Pritchard et al. (DOI: 10.1039/C5AY02550H) to compute ANI. If you publish your findings, \
                            please do not forget to properly credit their work.", lc='green', header="CITATION")

    def run_command(self, input_path, method='ANIb'):
        old_wd = os.getcwd()
        os.chdir(input_path)

        full_command = [self.program_name, '--outdir', 'output', '--indir', input_path, '-g', '-m', method]

        program = Popen(full_command, stdin=PIPE, stderr=PIPE)
        sStdout, sStdErr = program.communicate()

        if len(sStdErr) > 0 and 'ERROR' in str(sStdErr):
            print(sStdErr.decode('utf-8'))
            sys.exit(1)

        with open(os.path.join(input_path, 'output', method + '_percentage_identity.tab'), 'r') as f:
            percent_identity = f.read()

        os.chdir(old_wd)
        return percent_identity


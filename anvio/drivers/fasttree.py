# coding: utf-8
"""Interface to FastTree."""

import os
import io
from subprocess import Popen, PIPE

import anvio
import anvio.utils as utils
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths

from anvio.errors import ConfigError


__author__ = "Özcan Esen"
__copyright__ = "Copyright 2017, The anvio Project"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "Özcan Esen"
__email__ = "ozcanesen@gmail.com"


run = terminal.Run()
progress = terminal.Progress()
pp = terminal.pretty_print


class FastTree:
    def __init__(self):
        self.run = run
        self.progress = progress
        self.command = ['FastTree']

        utils.is_program_exists('FastTree')

    def create_tree_from_aligment_file(self, file_path):
        filesnpaths.is_file_fasta_formatted(file_path)
        return self.run_command(open(file_path, 'rb'))

    def create_tree_from_aligment_text(self, aligment_text):
        return self.run_command(io.BytesIO(aligment_text))

    def run_command(self, input_buffer):
        fasttree = Popen(self.command, stdout=PIPE, stdin=PIPE, stderr=PIPE)
        newick_tree = fasttree.communicate(input=input_buffer.read())[0].decode('utf-8')
        return newick_tree

# coding: utf-8
"""Interface to muscle."""

import os
import shutil

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
__maintainer__ = "A. Murat Eren"
__email__ = "a.murat.eren@gmail.com"


run = terminal.Run()
progress = terminal.Progress()
pp = terminal.pretty_print


class Muscle:
    def __init__(self, progress=progress, run=run, program_name = 'muscle'):
        """A class to take care of muscle alignments."""
        self.progress = progress
        self.run = run

        self.program_name = program_name

        utils.is_program_exists(self.program_name)

        self.citation = "Edgar, doi:10.1093/nar/gkh340"
        self.web = "http://www.drive5.com/muscle"


    def run_stdin(self, sequences_list, debug=False):
        """Takes a list of tuples for sequences, performs MSA using muscle, returns a dict.

            >>> from anvio.drivers.muscle import Muscle
            >>> m = Muscle()
            >>> m.run_stdin([('seq1', 'ATCATCATCGA'), ('seq2', 'ATCGAGTCGAT')])
            {u'seq1': u'ATCATCATCGA-', u'seq2': u'ATCG-AGTCGAT'}

        """

        tmp_dir = filesnpaths.get_temp_directory_path()
        log_file_path = os.path.join(tmp_dir, '00_log.txt')

        self.run.info('Running %s' % self.program_name, '%d seqeunces will be aligned' % len(sequences_list))
        self.run.info('Log file path', log_file_path)

        sequences_data = ''.join(['>%s\n%s\n' % (t[0], t[1]) for t in sequences_list])
        cmd_line = [self.program_name, '-quiet']

        output = utils.run_command_STDIN(cmd_line, log_file_path, sequences_data)

        if not (len(output) and output[0] == '>'):
            with open(log_file_path, "a") as log_file: log_file.write('# THIS IS THE OUTPUT YOU ARE LOOKING FOR:\n\n%s\n' % (output))
            raise ConfigError("Drivers::Muscle: Something went wrong with this alignment that was working on %d\
                               sequences :/ You can find the output in this log file: %s" % (len(sequences_list), log_file_path))

        alignments = {}

        # parse the output, and fill alignments
        defline, seq = None, None
        for line in [o for o in output.split('\n') if len(o)] + ['>']:
            if line.startswith('>'):
                if defline:
                    alignments[defline[1:]] = seq
                defline, seq = line, None
            else:
                if not seq:
                    seq = line
                else:
                    seq += line

        if not debug:
            shutil.rmtree(tmp_dir)

        return alignments

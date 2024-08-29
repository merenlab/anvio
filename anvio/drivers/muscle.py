# coding: utf-8
"""Interface to muscle."""

import os
import shutil

import anvio
import anvio.fastalib as f
import anvio.utils as utils
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths

from anvio.errors import ConfigError


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
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


    def run_stdin(self, sequences_list, debug=False, clustalw_format=False):
        """Takes a list of tuples for sequences, performs MSA using muscle, returns a dict.

            >>> from anvio.drivers.muscle import Muscle
            >>> m = Muscle()
            >>> m.run_stdin([('seq1', 'ATCATCATCGA'), ('seq2', 'ATCGAGTCGAT')])
            {u'seq1': u'ATCATCATCGA-', u'seq2': u'ATCG-AGTCGAT'}

        PARAMETERS
        ==========
        sequences_list : List of tuples
            each tuple should contain the sequence name as the first element, and the sequence itself 
            as the second element
        debug : Boolean
            controls whether or not the temporary directory is removed after execution. Probably a good
            idea to set this to anvio.DEBUG
        clustalw_format : Boolean
            if True, will return alignment in CLUSTALW format (obtained using Muscle's -clw flag) rather than
            a dict of FASTA alignments
        """

        tmp_dir = filesnpaths.get_temp_directory_path()
        log_file_path = os.path.join(tmp_dir, '00_log.txt')

        self.run.info('Running %s' % self.program_name, '%d sequences will be aligned' % len(sequences_list))
        self.run.info('Log file path', log_file_path)

        sequences_data = ''.join(['>%s\n%s\n' % (t[0], t[1]) for t in sequences_list])
        cmd_line = [self.program_name, '-quiet']

        if clustalw_format:
            cmd_line += ['-clw']

        additional_params = self.get_additional_params_from_shell()
        if additional_params:
            cmd_line += additional_params

        output = utils.run_command_STDIN(cmd_line, log_file_path, sequences_data)

        if not clustalw_format and not (len(output) and output[0] == '>'):
            with open(log_file_path, "a") as log_file: log_file.write('# THIS IS THE OUTPUT YOU ARE LOOKING FOR:\n\n%s\n' % (output))
            raise ConfigError("Drivers::Muscle: Something went wrong with this alignment that was working on %d "
                              "sequences :/ You can find the output in this log file: %s" % (len(sequences_list), log_file_path))

        if clustalw_format:
            return output

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


    def run_default(self, sequences_list, debug=False):
        """Takes a list of tuples for sequences, performs MSA using muscle, returns a dict.

        Unlike `run_stdin`, runs everything through files rather than passing data via STDIN.

            >>> from anvio.drivers.muscle import Muscle
            >>> m = Muscle()
            >>> m.run_default([('seq1', 'ATCATCATCGA'), ('seq2', 'ATCGAGTCGAT')])
            {u'seq1': u'ATCATCATCGA-', u'seq2': u'ATCG-AGTCGAT'}

        """

        tmp_dir = filesnpaths.get_temp_directory_path()
        log_file_path = os.path.join(tmp_dir, '00_log.txt')
        input_file_path = os.path.join(tmp_dir, 'input.fa')
        output_file_path = os.path.join(tmp_dir, 'output.fa')

        self.run.info('Running %s' % self.program_name, '%d sequences will be aligned' % len(sequences_list))
        self.run.info('Log file path', log_file_path)
        self.run.info('Input file path', input_file_path)
        self.run.info('Output file path', output_file_path)

        sequences_data = ''.join(['>%s\n%s\n' % (t[0], t[1]) for t in sequences_list])

        with open(input_file_path, 'w') as input_file:
            input_file.write(sequences_data)

        cmd_line = [self.program_name, '-in', input_file_path, '-out', output_file_path]

        additional_params = self.get_additional_params_from_shell()
        if additional_params:
            cmd_line += additional_params

        output = utils.run_command(cmd_line, log_file_path)

        if not os.path.exists(output_file_path) or os.path.getsize(output_file_path) == 0:
            with open(log_file_path, "a") as log_file: log_file.write('# THIS IS THE OUTPUT YOU ARE LOOKING FOR:\n\n%s\n' % (output))
            raise ConfigError("Drivers::Muscle: Something went wrong with this alignment that was working on %d "
                              "sequences :/ You can find the output in this log file: %s" % (len(sequences_list), log_file_path))

        alignments = {}

        # parse the output, and fill alignments
        output = f.SequenceSource(output_file_path)

        while next(output):
            alignments[output.id] = output.seq

        if not debug:
            shutil.rmtree(tmp_dir)

        return alignments


    def get_additional_params_from_shell(self):
        """Get additional user-defined params from environmental variables"""

        if 'MUSCLE_PARAMS' in os.environ:
            return os.environ['MUSCLE_PARAMS'].split()
        else:
            return None

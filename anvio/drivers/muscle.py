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

        self.major_version = self.get_major_version()


    def run_default(self, sequences_list, debug=False):
        """Takes a list of tuples for sequences, performs MSA using muscle, returns a dict.

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

        cmd_line = [self.program_name, '-align', input_file_path, '-output', output_file_path]

        additional_params = self.get_additional_params_from_shell()
        if additional_params:
            if '-super5' in additional_params:
                cmd_line = [self.program_name, '-super5', input_file_path, '-output', output_file_path]
                additional_params.remove('-super5')

            cmd_line += additional_params

        ret_val = utils.run_command(cmd_line, log_file_path)

        if ret_val:
            raise ConfigError("Drivers::Muscle: Something went wrong with this alignment that was working on %d "
                              "sequences :/ You can find the output in this log file: %s" % (len(sequences_list), log_file_path))

        if not os.path.exists(output_file_path) or os.path.getsize(output_file_path) == 0:
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


    def get_major_version(self):
        """Get the MUSCLE major version."""

        output, ret_code = utils.get_command_output_from_shell('%s -version' % self.program_name)
        output = output.decode('utf-8', errors='replace') if isinstance(output, bytes) else output
        output = output.lower()

        if output.startswith('muscle 5') or output.startswith('muscle v5'):
            return 5

        if output.startswith('muscle v3'):
            raise ConfigError("Anvi'o recently started using a newer version of MUSCLE (you know, the "
                              "sequence alignment software), but your installed version in this environment "
                              "appears have the old version still :/ You can solve this issue by simply "
                              "installing MUSCLE v5 or newer. If you are in a conda environment, you can "
                              "try running the following: `conda install -c conda-forge -c bioconda ""\"muscle>=5\"`.")

        raise ConfigError("The anvi'o MUSCLE driver requires MUSCLE 5, but failed to recognize the installed "
                          "MUSCLE version from `muscle -version` output: %s" % output)


    def get_additional_params_from_shell(self):
        """Get additional user-defined params from environmental variables"""

        if 'MUSCLE_PARAMS' in os.environ:
            return os.environ['MUSCLE_PARAMS'].split()
        else:
            return None

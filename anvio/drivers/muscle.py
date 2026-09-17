"""Interface to muscle."""

import os
import re
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

# asked about muscle version only once per process rather than once per alignment
major_version_cache = {}


class Muscle:
    def __init__(self, progress=progress, run=run, program_name = 'muscle', num_threads=1):
        """A class to take care of muscle alignments."""
        self.progress = progress
        self.run = run

        self.program_name = program_name

        # MUSCLE 5 helps itself to every core on the machine unless it is told otherwise, and
        # anvi'o already runs many of these alignments in parallel. one thread per alignment
        # leaves the parallelism to the caller, who is the one that knows how much of the
        # machine it is using.
        self.num_threads = num_threads

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

        # a user who asks for a specific number of threads through MUSCLE_PARAMS gets it
        if '-threads' not in cmd_line:
            cmd_line += ['-threads', str(self.num_threads)]

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

        if self.program_name in major_version_cache:
            return major_version_cache[self.program_name]

        output, ret_code = utils.get_command_output_from_shell('%s -version' % self.program_name)
        output = output.decode('utf-8', errors='replace') if isinstance(output, bytes) else output

        # MUSCLE reports itself as `muscle 5.3.osx64 []` or `MUSCLE v3.8.1551 by Robert C. Edgar`.
        # anything else the shell has to say (an environment activation notice, a warning about a
        # library) is merged into this output too, so the version is looked for on every line
        # rather than only at the very beginning of it.
        major_version = None
        for line in output.lower().splitlines():
            version_match = re.search(r'^muscle\s+v?(\d+)', line.strip())
            if version_match:
                major_version = int(version_match.group(1))
                break

        if major_version == 5:
            major_version_cache[self.program_name] = major_version
            return major_version

        if major_version and major_version < 5:
            raise ConfigError("Anvi'o recently started using a newer version of MUSCLE (you know, the "
                              "sequence alignment software), but the one in this environment is still "
                              "MUSCLE v%d :/ You can solve this issue by simply installing MUSCLE v5. "
                              "If you are in a conda environment, you can try running the following: "
                              "`conda install -c conda-forge -c bioconda \"muscle>=5,<6\"`." % major_version)

        if major_version:
            raise ConfigError("The anvi'o MUSCLE driver knows how to talk to MUSCLE v5, and the one in this "
                              "environment is MUSCLE v%d :/ Every major version of MUSCLE so far has come "
                              "with its own command line, so anvi'o would rather say this out loud than "
                              "assume it knows how to run this one and risk making a mess of your "
                              "alignments. Installing MUSCLE v5 will get you going. We would also love to "
                              "hear about this at https://github.com/merenlab/anvio/issues so the driver "
                              "can catch up with MUSCLE." % major_version)

        raise ConfigError("The anvi'o MUSCLE driver requires MUSCLE v5, but it could not tell which "
                          "version of MUSCLE is installed in this environment. This is what `%s -version` "
                          "had to say for itself (with an exit code of %d): \"%s\"." \
                                % (self.program_name, ret_code, output.strip() or '(nothing at all)'))


    def get_additional_params_from_shell(self):
        """Get additional user-defined params from environmental variables"""

        if 'MUSCLE_PARAMS' in os.environ:
            return os.environ['MUSCLE_PARAMS'].split()
        else:
            return None

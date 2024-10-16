# coding: utf-8
"""An interface for tRNAScan-SE"""

import anvio
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


class tRNAScanSE:
    """Parent class for all tRNAScan-SE usage

    Parameters
    ==========
    args: Namespace object
        All the arguments from up above.
    program_name: str, tRNAscan-SE
        The program name that gives access to tRNAscan-SE functionality
    """

    def __init__(self, args, program_name='tRNAscan-SE', run=None, progress=None, skip_sanity_check=False):
        self.program_name = program_name

        self.tested_versions = ['2.0.5', '2.0.7', '2.0.9', '2.0.12']

        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        self.num_threads = A('num_threads') or 1
        self.fasta_file_path = A('fasta_file')
        self.trna_hits_file_path = A('trna_hits_file')
        self.log_file_path = A('log_file')
        self.cutoff_score = A('trna_cutoff_score') or 20
        self.quiet = A('quiet')

        self.run = run or terminal.Run(verbose=(not self.quiet))
        self.progress = progress or terminal.Progress(not self.quiet)

        if not skip_sanity_check:
            self.sanity_check()


    def sanity_check(self):
        """Executes rudimentary checks

        Parameters
        ==========
        N\A

        Returns
        =======
        N\A
        """

        self.check_programs()

        for name, variable in [('a log file', self.log_file_path), ('a FASTA file path', self.fasta_file_path)]:
            if not variable:
                raise ConfigError("A proper instance of this driver must have %s variable set." % name)

        filesnpaths.is_output_file_writable(self.log_file_path)
        filesnpaths.is_file_exists(self.fasta_file_path)

        try:
            self.cutoff_score = float(self.cutoff_score)
        except:
            raise ConfigError("Cutoff score must be a float.")

        if self.cutoff_score < 20 or self.cutoff_score > 100:
            raise ConfigError("The cutoff score must be between 20 and 100.")


    def check_programs(self, quiet=False):
        utils.is_program_exists(self.program_name)

        output, ret_code = utils.get_command_output_from_shell('%s -h' % self.program_name)
        try:
            version_found = output.split(b'\n')[1].split()[1].split(b':')[0].lower().decode("utf-8")
            if not quiet:
                self.run.info('%s version found' % self.program_name, version_found, mc="green", nl_after=1)
        except:
            version_found = 'Unknown'
            self.run.warning("Anvi'o failed to learn the version of %s installed on this system :/" % self.program_name)

        if version_found not in self.tested_versions:
            self.run.warning("The version of %s installed on your system ('%s') is not one of those that we tested its anvi'o driver "
                             "with. Anvi'o will continue to try to run everything as if this didn't happen. If you see this warning "
                             "but everything works fine, let us know so we can include this version number into the list of 'tested' "
                             "version numbers. If you see an unexpected error, please consider installing one of these versions "
                             "of tRNAScan-SE (and again please let us know anyway so we can address it for later): '%s'" % \
                                           (self.program_name, version_found, ', '.join(list(self.tested_versions))))

        self.installed_version = version_found


    def parse_output(self):
        """Parses the tRNAScan-SE output file

        Parameters
        ==========
        N\A

        Returns
        =======
        results: `dict`
            A Dictionary of hits
        """

        num_lines = filesnpaths.get_num_lines_in_file(self.trna_hits_file_path)

        if not num_lines:
            self.run.warning("No tRNA genes found in tRNAScan-SE output.")
            return {}

        d = {}

        self.progress.new("Parsing the output ...")
        with open(self.trna_hits_file_path) as output:
            # first three lines are garbage
            for i in range(0, 3):
                output.readline()

            entry_no = 0
            num_introns = 0
            while 1:
                self.progress.update(entry_no)
                line = output.readline().strip('\n')

                if not line:
                    break

                entry_no += 1

                fields = [f.strip() for f in line.split('\t')]

                if not len(fields) == 10:
                    raise ConfigError("The expected output of tRNAScan-SE includes exactly 10 columns. However, the output "
                                      "anvi'o is working contains at least one line with %d columns :/ This doesn't look "
                                      "good. Here is the list of columns data of that line for your reference: '%s'." \
                                                            % (len(fields), fields))

                d[entry_no] = {'contig_name': fields[0],
                               'trna_no': fields[1],
                               'start': int(fields[2]),
                               'stop': int(fields[3]),
                               'amino_acid': fields[4],
                               'anticodon': fields[5],
                               'intron_start': int(fields[6]),
                               'intron_end': int(fields[7]),
                               'score': float(fields[8])}

                if d[entry_no]['intron_start']:
                    num_introns += 1
                    if d[entry_no]['start'] < d[entry_no]['stop']:
                        d[entry_no]['intron_start'] = int(fields[6]) - d[entry_no]['start']
                        d[entry_no]['intron_end'] = int(fields[7]) - int(fields[6])
                    else:
                        d[entry_no]['intron_start'] = d[entry_no]['start'] - int(fields[6])
                        d[entry_no]['intron_end'] = int(fields[6]) - int(fields[7]) 

        self.progress.end()

        self.run.info("Num tRNA genes parsed", entry_no)
        self.run.info("Num tRNA genes with introns", num_introns, nl_after=1)

        return d


    def process(self):
        """Processes a given FASTA file with tRNAScan-SE and reports all hits

        Parameters
        ==========
        N\A

        Returns
        =======
        results: `dict`
            A Dictionary of hits
        """

        filesnpaths.is_output_file_writable(self.trna_hits_file_path, ok_if_exists=False)

        command = [self.program_name,
                   self.fasta_file_path,
                   '--score', self.cutoff_score,
                   '-G',
                   '-o', self.trna_hits_file_path,
                   '--thread', self.num_threads]

        self.run.warning("Anvi'o will use 'tRNAScan-SE' by Chan and Lowe (doi:10.1007/978-1-4939-9173-0_1) to identify tRNA "
                         "sequences in your data. When you publish your findings, please do not forget to properly credit "
                         "their work.", lc='green', header="CITATION")

        self.progress.new('Processing')
        self.progress.update('Running tRNAScan-SE scan in %d threads...' % self.num_threads)

        exit_code = utils.run_command(command, self.log_file_path)

        self.progress.end()

        if exit_code:
            raise ConfigError("tRNAScan-SE finished with a non-zero exit code, which indicates that something went "
                              "wrong with it :/ Please check the log file to see learn more :/")

        d = self.parse_output()

        return d



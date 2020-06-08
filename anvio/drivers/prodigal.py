# coding: utf-8
"""Interface to Prodigal."""

# Note: The Prodigal command will be run in a "mutlithreaded" way even if the user selects only one thread.  In other
#   words, it will still "split" the input fasta into a single split, and run Prodigal in a new thread.  It seems
#   wasteful, but with the Infant Gut test dataset, it only added a couple of seconds of overhead vs. the original
#   single threaded code.  So it seems to not be a major performance issue as of now.

import os

import anvio
import anvio.utils as utils
import anvio.fastalib as fastalib
import anvio.terminal as terminal
import anvio.constants as constants
import anvio.filesnpaths as filesnpaths

from anvio.anvi_thread import AnviThread
from anvio.errors import CommandError, ConfigError

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

# How many times to retry the gene calling.
NUM_RETRIES = 3


def run_command_wrapper(command, log_file_path, run, retries=NUM_RETRIES):
    """A wrapper around `utils.run_command`.

    It will retry the command `retries` time.  Sometimes commands can fail for weird reasons so retrying often fixes
    things.

    If the return value of the command is zero, then it is returned.  If not, an instance of CommandError is returned.
    """
    for attempt in range(retries):
        try:
            ret_val = utils.run_command(command, log_file_path)
        except ConfigError as e:
            # utils.run_command will raise ConfigError on certain types of failures.
            run.warning(f'Caught error: {e}, while running command {command}, attempt {attempt}, '
                        'but we are going to keep trying!')
            continue

        # Return code of zero means success in the unix world.
        if ret_val == 0:
            return ret_val

    # Failed to run the command a couple of times, so return an error.
    return CommandError(f'Failed to run the command after {retries} tries. Final exit code was: {ret_val}')


def split_fasta_wrapper(fasta_file_path, num_splits):
    """A wrapper around `utils.split_fasta`.  It checks that the split files actually exist.

    Returns output path.
    """
    paths = utils.split_fasta(fasta_file_path, parts=num_splits)

    # Make sure all the splits actually exist.
    for path in paths:
        assert os.path.exists(path), f'I expected {path} to exist, but it did not.  Something weird happened.'

    return paths


def add_suffix(strings, suffix, separator='.'):
    """Adds a suffix on to a list of `strings`.  Joins with `separator`."""
    return [separator.join([string, suffix]) for string in strings]


# Todo: Consider using shutil.copyfileobj instead
def concatenate_files(intermediate_files, final_file):
    """Concatenate `intermediate_files` into `final_file`.

    Raises `CommandError` if any of the `intermediate_files` don't exist.
    """
    with open(final_file, 'a') as outfile:
        for file_path in intermediate_files:
            if os.path.exists(file_path):
                with open(file_path, 'r') as infile:
                    # Read line by line in case the logfile is large.
                    for line in infile:
                        outfile.write(line)
            else:
                raise CommandError('Something went wrong in the command, and we could not find the expected '
                                   f'logfile: {file_path}')


def run_commands_in_threads(commands, target, final_logfile, num_retries):
    """Runs the given `commands` each in its own thread.

    `target` is the function/method to run in the thread (e.g., run_command_wrapper).

    `log_file_path` will contain logfiles from each thread, concatenated.

    Commands are tried `num_retries` times each, as often bash commands can fail for weird reasons and sometimes rerunning them solves the issue.
    """
    threads = []
    intermediate_logfiles = []
    for index, command in enumerate(commands):
        # Each thread will get its own logfile so we don't have to worry about locking on the main logfile.
        this_logfile_path = '.'.join([final_logfile, f'thread-{index}'])
        intermediate_logfiles.append(this_logfile_path)

        # Spawn a new thread and start it.
        thread = AnviThread(target=target, args=(command, this_logfile_path, num_retries))
        thread.start()
        threads.append(thread)

    # Wait for all jobs to finish.
    for thread in threads:
        thread.join()

    concatenate_files(intermediate_logfiles, final_logfile)

    for file_path in intermediate_logfiles:
        os.remove(file_path)

    return threads


class Prodigal:
    def __init__(self, args=None, progress=progress, run=run):
        self.progress = progress
        self.run = run
        self.args = args

        A = lambda x: (args.__dict__[x] if x in args.__dict__ else None) if args else None
        self.prodigal_translation_table = A('prodigal_translation_table')

        self.num_threads = A('num_threads')

        self.run.info('No. threads for gene calling', self.num_threads)

        self.parser = None
        self.installed_version = None
        self.available_parsers = {'v2.6.3': self.__parser_1,
                                  'v2.6.2': self.__parser_1,
                                  'v2.60': self.__parser_1}

        self.check_version()

    def __parser_1(self, defline):
        """parses this:

            204_10M_MERGED.PERFECT.gz.keep_contig_1720_7 # 7086 # 7262 # 1 # ID=3_7;partial=00;start_type=ATG;rbs_motif=GGAG/GAGG;rbs_spacer=5-10bp;gc_cont=0.294

        """

        # FIXME: This bullshit is here due to the ridiculous design of GFF3 format, for which I don't
        #        want to write a generic parser. What is available to parse GFF3 files are not any impressive
        #        than the file format itself either. I tried one library through PyPI, which failed to parse
        #        prodigal GFF3 output. Seriously. This is very upsetting. Is it even possible to build something
        #        long-lasting for proper bioinformatics? ... This class will need to be replaced with something more
        #        reasonable at some point. The best solution would be to have a generic GFF3 parser that maps that
        #        BS into a lightweight data structure for anvi'o internals.

        hit = {}

        fields = defline.split()

        if not len(fields) != 8 or fields[6] not in ['1', '-1']:
            raise ConfigError('Somethings is wrong with this prodigal output :( The parser for the '
                              'version %s is failing to make sense of it.' % (self.installed_version))

        hit['direction'] = 'f' if fields[6] == '1' else 'r'
        hit['start'] = int(fields[2]) - 1
        hit['stop'] = int(fields[4])
        hit['contig'] = '_'.join(fields[0].split('_')[:-1])

        additional_attributes = fields[8].split(';')

        hit['partial'] = False if additional_attributes[1] == 'partial=00' else True

        hit['source'] = 'prodigal'
        hit['call_type'] = constants.gene_call_types['CODING']
        hit['version'] = self.installed_version

        return hit

    def __make_prodigal_commands(self, input_paths, gff_paths, peptide_paths):
        prodigal_commands = []
        for i in range(len(input_paths)):
            command = ['prodigal',
                       '-i', input_paths[i],
                       '-o', gff_paths[i],
                       '-a', peptide_paths[i]]

            # Use either custom translation table or 'meta' mode.
            if self.prodigal_translation_table:
                command.extend(['-g', self.prodigal_translation_table])
                self.run.warning("Prodigal translation table is set to '%s' (whatever you did has worked so far, but "
                                 "keep an eye for errors from prodigal in case it doesn't like your translation table "
                                 "parameter). This means we will not use prodigal in the metagenomics mode, due to this"
                                 "issue: https://github.com/hyattpd/Prodigal/issues/19. If that issue is closed, and "
                                 "you are reading this message, then please contact an anvi'o developer."
                                 % str(self.prodigal_translation_table))
            else:
                command.extend(['-p', 'meta'])

            prodigal_commands.append(command)

        return prodigal_commands

    def __check_for_prodigal_errors(self, threads):
        """The `threads` have information regarding errors during gene calling.  This function checks for those
        errors.
        """
        for thread in threads:
            if isinstance(thread.target_return_value, CommandError):
                self.progress.end()
                raise thread.target_return_value
            elif thread.target_return_value != 0:
                # Technically should never get to this point as `target_return_value` will either be a CommandError
                # or zero.
                #
                # Todo replace BaseException
                self.progress.end()
                raise BaseException('Something went horribly wrong.')

    # Check the gff files and append them into the single outfile.
    def __process_gff_files(self, gff_paths):
        """Check the gff files and combine them.

        Note:  is genes_in_contigs (the gff file) ever actually used...each of the threads will start renumbering
        sequences from 1 in each subfile.  Will this break something downstream?
        """
        for path in gff_paths:
            if not os.path.exists(path):
                self.progress.end()
                raise ConfigError(f'Something bad happened and {path} does not exist!')

            # Some splits may not actually have gene calls.  Skip them.
            if filesnpaths.is_file_empty(path):
                continue

            with open(self.genes_in_contigs, 'a') as outfile:
                with open(path) as infile:
                    for line in infile:
                        outfile.write(line)

            os.remove(path)

    def __process_peptide_files(self, peptide_paths, log_file_path):
        """Checks that `peptide_paths` files actually exist, then combines them."""
        # For renaming fasta headers
        hit_id = 0

        # Set up data storage.
        gene_calls_dict = {}  # each entry must contain {'contig', 'start', stop, 'direction', 'partial'} items.
        amino_acid_sequences_dict = {}

        for peptide_path in peptide_paths:
            if not os.path.exists(peptide_path):
                self.progress.end()
                raise ConfigError("Something went wrong with prodigal, and it failed to generate the "
                                  "expected output ('%s') :/ Fortunately, this log file should tell you what "
                                  "might be the problem: '%s'. Please do not forget to include this "
                                  "file if you were to ask for help." % (peptide_path, log_file_path))

            # Some splits may not actually have gene calls.  Skip them.
            if filesnpaths.is_file_empty(peptide_path):
                continue

            # If we get here, the fasta file will not be empty.
            fasta = fastalib.SequenceSource(peptide_path)

            while next(fasta):
                gene_calls_dict[hit_id] = self.parser(fasta.id)
                amino_acid_sequences_dict[hit_id] = fasta.seq.replace('*', '')
                hit_id += 1

            fasta.close()

            # Remove the split peptide file.
            os.remove(peptide_path)

        # If no genes were predicted across all output files, warn the user.
        if len(amino_acid_sequences_dict) == 0:
            self.run.info('Result',
                          f'Prodigal ({self.installed_version}) has identified no genes :/',
                          nl_after=1,
                          mc="red")
        else:  # Write out the final gene file
            with open(self.amino_acid_sequences_in_contigs, 'w') as f:
                for hit_id, sequence in amino_acid_sequences_dict.items():
                    f.write(f">{hit_id}\n{sequence}\n")

        return gene_calls_dict, amino_acid_sequences_dict

    def check_version(self):
        """checks the installed version of prodigal, sets the parser"""

        utils.is_program_exists('prodigal')
        output, ret_code = utils.get_command_output_from_shell('prodigal -v')

        version_found = output.split(b'\n')[1].split()[1].split(b':')[0].lower().decode("utf-8")

        if version_found not in self.available_parsers:
            raise ConfigError("The prodigal version installed on your system is not compatible "
                              "with any of the versions anvi'o can work with. Please install "
                              "any of the following versions: %s" % (', '.join(list(self.available_parsers.keys()))))

        self.installed_version = version_found
        self.parser = self.available_parsers[version_found]

    def process(self, fasta_file_path, output_dir):
        """Take the fasta file, run prodigal on it, and make sense of the output

           Returns a gene calls dict, and amino acid sequences dict.
        """

        # Generate the needed input/output file paths.
        input_paths = split_fasta_wrapper(fasta_file_path, self.num_threads)
        assert len(input_paths) == self.num_threads

        gff_paths = add_suffix(input_paths, 'prodigal-gff')
        peptide_paths = add_suffix(input_paths, 'prodigal-peptides')
        assert len(input_paths) == len(gff_paths) and len(input_paths) == len(peptide_paths)

        # Commands to run prodigal on fastas in `input_paths`.  Output data will be in `gff_paths` and `peptide_paths`.
        prodigal_commands = self.__make_prodigal_commands(input_paths, gff_paths, peptide_paths)

        self.genes_in_contigs = os.path.join(output_dir, 'contigs.genes')
        self.amino_acid_sequences_in_contigs = os.path.join(output_dir, 'contigs.amino_acid_sequences')

        log_file_path = os.path.join(output_dir, '00_log.txt')

        self.run.warning('', header='Finding ORFs in contigs', lc='green')
        self.run.info('Genes', self.genes_in_contigs)
        self.run.info('Amino acid sequences', self.amino_acid_sequences_in_contigs)
        self.run.info('Log file', log_file_path)

        self.run.warning("Anvi'o will use 'prodigal' by Hyatt et al (doi:10.1186/1471-2105-11-119) to identify open "
                         "reading frames in your data. When you publish your findings, please do not forget to "
                         "properly credit their work.", lc='green', header="CITATION")

        self.progress.new('Processing')
        self.progress.update('Identifying ORFs in contigs ...')

        threads = run_commands_in_threads(prodigal_commands, run_command_wrapper, log_file_path, NUM_RETRIES)

        self.progress.update('Checking for Prodigal errors ...')
        self.__check_for_prodigal_errors(threads)

        self.progress.update('Processing GFF files ...')
        self.__process_gff_files(gff_paths)

        self.progress.update('Processing peptide files ...')
        gene_calls_dict, amino_acid_sequences_dict = self.__process_peptide_files(peptide_paths, log_file_path)

        self.progress.update('Removing fasta splits ...')
        for path in input_paths:
            os.remove(path)

        self.progress.end()

        self.run.info('Result',
                      'Prodigal (%s) has identified %d genes.' % (self.installed_version, len(gene_calls_dict)),
                      nl_after=1)

        return gene_calls_dict, amino_acid_sequences_dict

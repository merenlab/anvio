# coding: utf-8
"""Interface to Prodigal."""

# Note: The Prodigal command will be run in a "mutlithreaded" way even if the user selects only one thread.  In other
#   words, it will still "split" the input fasta into a single split, and run Prodigal in a new thread.  It seems
#   wasteful, but with the Infant Gut test dataset, it only added a couple of seconds of overhead vs. the original
#   single threaded code.  So it seems to not be a major performance issue as of now.
import argparse
import os

import anvio
import anvio.utils as utils
import anvio.terminal as terminal
import anvio.constants as constants

from anvio.errors import ConfigError

from anvio.threadingops import ThreadedProdigalRunner

__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "A. Murat Eren"
__email__ = "a.murat.eren@gmail.com"

run = terminal.Run()
progress = terminal.Progress()
pp = terminal.pretty_print


class Prodigal:
    def __init__(self, args=None, progress=progress, run=run):
        self.progress = progress
        self.run = run
        self.args = args
        A = lambda x: (args.__dict__[x] if x in args.__dict__ else None) if args else None
        self.prodigal_translation_table = A('prodigal_translation_table')
        self.num_threads = A('num_threads')
        self.prodigal_single_mode = A('prodigal_single_mode')

        self.run.info('Num threads for gene calling', self.num_threads)

        self.parser = None
        self.installed_version = None
        self.available_parsers = {'v2.6.3': self.__parser_1,
                                  'v2.6.2': self.__parser_1,
                                  'v2.6.0': self.__parser_1}

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

        # Set up the output files.
        log_file_path = os.path.join(output_dir, '00_log.txt')
        self.genes_in_contigs = os.path.join(output_dir, 'contigs.genes')
        self.amino_acid_sequences_in_contigs = os.path.join(output_dir, 'contigs.amino_acid_sequences')

        self.run.warning("Anvi'o will use 'prodigal' by Hyatt et al (doi:10.1186/1471-2105-11-119) to identify open "
                         "reading frames in your data. When you publish your findings, please do not forget to "
                         "properly credit their work.", lc='green', header="CITATION")

        # Put some nice logging info.
        self.run.warning('', header='Finding ORFs in contigs using prodigal', lc='green')
        self.run.info('Procedure', 'single' if self.prodigal_single_mode else 'meta')
        self.run.info('Genes', self.genes_in_contigs)
        self.run.info('Amino acid sequences', self.amino_acid_sequences_in_contigs)
        self.run.info('Log file', log_file_path)

        self.progress.new('Processing')
        self.progress.update(f"Identifying ORFs using {terminal.pluralize('thread', self.num_threads)}.")

        # if more threads assigned to gene calling than the number of sequences in the FASTA file,
        # it can cause issues downstream, and we should set the number of threads accordingly.
        num_sequences_in_fasta_file = utils.get_num_sequences_in_fasta(fasta_file_path)
        if num_sequences_in_fasta_file < self.num_threads:
            self.progress.reset()
            self.run.warning(f"Even though you set the number of threads to {self.num_threads}, your FASTA file contains only "
                             f"{terminal.pluralize('sequence', num_sequences_in_fasta_file)}. "
                             f"To avoid any hiccups later, anvi'o will set the number of threads to match the number of "
                             f"sequences in your FASTA (who would've thought a perfect assembly can have a downside?).")
            self.num_threads = num_sequences_in_fasta_file

            self.progress.update(f"Identifying ORFs using {terminal.pluralize('thread', self.num_threads)}")

        # Set up the prodigal runner.
        collated_output_file_paths = {'gff_path': self.genes_in_contigs,
                                      'peptide_path': self.amino_acid_sequences_in_contigs}

        args = argparse.Namespace(input_file_path=fasta_file_path,
                                  collated_output_file_paths=collated_output_file_paths,
                                  number_of_splits=self.num_threads,
                                  log_file_path=log_file_path,
                                  logger=terminal.Logger(progress=self.progress, run=self.run),
                                  installed_version=self.installed_version,
                                  parser=self.parser,
                                  translation_table=self.prodigal_translation_table,
                                  prodigal_single_mode=self.prodigal_single_mode)

        prodigal_runner = ThreadedProdigalRunner(args)

        # Run the pipeline!
        state = prodigal_runner.run()

        self.progress.end()
        self.run.info('Result',
                      'Prodigal (%s) has identified %d genes.' % (self.installed_version,
                                                                  len(state['gene_calls_dict'])),
                      nl_after=1)

        return state['gene_calls_dict'], state['amino_acid_sequences_dict']

# coding: utf-8
"""Interface to Prodigal."""

import os

import anvio
import anvio.utils as utils
import anvio.fastalib as fastalib
import anvio.terminal as terminal

from anvio.errors import ConfigError


__author__ = "A. Murat Eren"
__copyright__ = "Copyright 2016, The anvio Project"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "A. Murat Eren"
__email__ = "a.murat.eren@gmail.com"


run = terminal.Run()
progress = terminal.Progress()
pp = terminal.pretty_print


class Prodigal:
    def __init__(self, progress=progress, run=run):
        self.progress = progress
        self.run = run

        self.parser = None
        self.installed_prodigal_version = None
        self.ok_prodigal_versions = {'v2.6.2': self.__parser_1,
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
            raise ConfigError, 'Somethings is wrong with this prodigal output :( The parser for the\
                                version %s is failing to make sense of it.' % (self.installed_prodigal_version)

        hit['direction'] = 'f' if fields[6] == '1' else 'r'
        hit['start'] = int(fields[2]) - 1
        hit['stop'] = int(fields[4])
        hit['contig'] = '_'.join(fields[0].split('_')[:-1])

        additional_attributes = fields[8].split(';')

        hit['partial'] = False if additional_attributes[1] == 'partial=00' else True

        hit['source'] = 'prodigal'
        hit['version'] = self.installed_prodigal_version

        return hit


    def check_version(self):
        """checks the installed version of prodigal, sets the parser"""

        utils.is_program_exists('prodigal')
        output, ret_code = utils.get_command_output_from_shell('prodigal -v')

        prodigal_version_found = output.split('\n')[1].split()[1].split(':')[0].lower()

        if prodigal_version_found not in self.ok_prodigal_versions:
            raise ConfigError, "The prodigal version installed on your system is not compatible\
                                with any of the versions anvi'o can work with. Please install\
                                any of the following versions: %s" % (', '.join(self.ok_prodigal_versions.keys()))

        self.installed_prodigal_version = prodigal_version_found
        self.parser = self.ok_prodigal_versions[prodigal_version_found]


    def process(self, fasta_file_path, output_dir):
        """Take the fasta file, run prodigal on it, and make sense of the output

           Returns a gene calls dict, and protein sequences dict.
        """
        gene_calls_dict = {} # each entry must contain {'contig', 'start', stop, 'direction', 'partial'} items.
        protein_sequences_dict = {}

        self.genes_in_contigs = os.path.join(output_dir, 'contigs.genes')
        self.proteins_in_contigs = os.path.join(output_dir, 'contigs.proteins')

        log_file_path = os.path.join(output_dir, '00_log.txt')

        self.run.warning('', header='Finding ORFs in contigs', lc='green')
        self.run.info('Genes', self.genes_in_contigs)
        self.run.info('Proteins', self.proteins_in_contigs)
        self.run.info('Log file', log_file_path)

        self.progress.new('Processing')
        self.progress.update('Identifying ORFs in contigs ...')
        cmd_line = ('prodigal -i "%s" -o "%s" -a "%s" -p meta >> "%s" 2>&1' % (fasta_file_path,
                                                                               self.genes_in_contigs,
                                                                               self.proteins_in_contigs,
                                                                               log_file_path))
        with open(log_file_path, "a") as myfile: myfile.write('CMD: ' + cmd_line + '\n')
        utils.run_command(cmd_line)

        if not os.path.exists(self.proteins_in_contigs):
            self.progress.end()
            raise ConfigError, "Something went wrong with prodigal, and it failed to generate the\
                                expected output :/ Fortunately, this log file should tell you what\
                                might be the problem: '%s'. Please do not forget to include this\
                                file if you were to ask for help." % log_file_path

        self.progress.update('Processing gene calls ...')

        fasta = fastalib.SequenceSource(self.proteins_in_contigs)

        hit_id = 0
        while fasta.next():
            gene_calls_dict[hit_id] = self.parser(fasta.id)
            protein_sequences_dict[hit_id] = fasta.seq.replace('*', '')
            hit_id += 1

        fasta.close()

        self.progress.end()

        self.run.info('Result', 'Prodigal (%s) has identified %d genes.' % (self.installed_prodigal_version, len(gene_calls_dict)), nl_after=1)

        return gene_calls_dict, protein_sequences_dict

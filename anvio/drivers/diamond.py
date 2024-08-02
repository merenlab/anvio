# coding: utf-8
"""Interface to Diamond."""

import os
import pandas as pd
import tempfile

import anvio
import anvio.utils as utils
import anvio.terminal as terminal

from anvio.errors import ConfigError


__copyright__ = "Copyleft 2015-2018, the Meren Lab (http://merenlab.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "A. Murat Eren"
__email__ = "a.murat.eren@gmail.com"


run = terminal.Run()
progress = terminal.Progress()
pp = terminal.pretty_print


class Diamond:
    def __init__(self, query_fasta=None, target_fasta=None, outfmt=None, run=run, progress=progress, num_threads=1, overwrite_output_destinations=False):
        self.run = run
        self.progress = progress

        self.num_threads = num_threads
        self.overwrite_output_destinations = overwrite_output_destinations

        utils.is_program_exists('diamond')

        self.tmp_dir = tempfile.gettempdir()
        self.evalue = 1e-05
        self.min_pct_id = None
        self.max_target_seqs = 100000
        self.outfmt = outfmt or '6'

        self.query_fasta = query_fasta
        self.target_fasta = target_fasta

        if not self.target_fasta:
            self.target_fasta = self.query_fasta

        self.search_output_path = 'diamond-search-results'
        self.tabular_output_path = 'diamond-search-results.txt'

        if not self.run.log_file_path:
            self.run.log_file_path = 'diamond-log-file.txt'

        self.additional_params_for_blastp = ""

        # if names_dict is None, all fine. if not, the query_fasta is assumed to be uniqued, and names_dict is
        # the dictionary that connects the ids in the fasta file, to ids that were identical to it.
        self.names_dict = None


    def get_blast_results(self):
        force_makedb, force_blastp, force_view = False, False, False

        if self.overwrite_output_destinations:
            force_makedb = True

        if os.path.exists(self.target_fasta + '.dmnd') and not force_makedb:
            self.run.warning("Notice: A diamond database is found in the output directory, and will be used!")
        else:
            self.makedb()
            force_blastp, force_view = True, True

        if os.path.exists(self.tabular_output_path) and not force_blastp:
            self.run.warning("Notice: A DIAMOND search result is found in the output directory: skipping BLASTP!")
        else:
            self.blastp()
            force_view = True

        if os.path.exists(self.tabular_output_path) and not force_view:
            self.run.warning("Notice: A DIAMOND tabular output is found in the output directory. Anvi'o will not generate another one!")
        else:
            self.view()

        return self.tabular_output_path


    def check_output(self, expected_output, process='diamond'):
        if not os.path.exists(expected_output):
            self.progress.end()
            raise ConfigError("Pfft. Something probably went wrong with Diamond's '%s' since one of the expected output files are missing. "
                               "Please check the log file here: '%s'. IT IS VERY LIKELY to get these kinds of errors if the version of "
                               "DIAMOND installed on your system differs from the one you had used to first setup your databases. Some "
                               "errors may disappear if you were to setup your search databases from scratch." % (process, self.run.log_file_path))


    def makedb(self, output_file_path=None):
        self.run.warning(None, header="DIAMOND MAKEDB", lc="green")
        self.progress.new('DIAMOND')
        self.progress.update('creating the search database (using %d thread(s)) ...' % self.num_threads)

        # NOTE Question from Evan. Why is the query_fasta the input to the database?
        cmd_line = ['diamond',
                    'makedb',
                    '--in', self.query_fasta,
                    '-d', output_file_path or self.target_fasta,
                    '-p', self.num_threads]

        utils.run_command(cmd_line, self.run.log_file_path)

        self.progress.end()

        expected_output = (output_file_path or self.target_fasta) + '.dmnd'

        self.run.info('Command line', ' '.join([str(x) for x in cmd_line]), quiet=True)
        self.run.info('Diamond search DB', expected_output)


    def blastp(self):
        self.run.warning(None, header="DIAMOND BLASTP", lc="green")
        self.run.info("Additional params for blastp", self.additional_params_for_blastp, mc='green')

        cmd_line = ['diamond',
                    'blastp',
                    '-q', self.query_fasta,
                    '-d', self.target_fasta,
                    '-o', self.tabular_output_path,
                    '-t', self.tmp_dir,
                    '-p', self.num_threads,
                    '--outfmt', *self.outfmt.split()]

        if self.additional_params_for_blastp:
            cmd_line.extend(self.additional_params_for_blastp.split())

        if self.max_target_seqs:
            cmd_line.extend(['--max-target-seqs', self.max_target_seqs])

        if self.min_pct_id:
            cmd_line.extend(['--id', self.min_pct_id])

        if self.evalue:
            cmd_line.extend(['--evalue', self.evalue])

        self.run.info('DIAMOND blastp cmd', ' '.join([str(p) for p in cmd_line]), quiet=(not anvio.DEBUG))

        self.progress.new('DIAMOND')
        self.progress.update('Running blastp (using %d thread(s)) ...' % self.num_threads)

        utils.run_command(cmd_line, self.run.log_file_path)

        self.progress.end()

        self.run.info('Search results', self.tabular_output_path)


    def blastp_stdout(self):
        self.run.warning(None, header="DIAMOND BLASTP STDOUT", lc="green")
        self.run.info("Additional params for blastp", self.additional_params_for_blastp, mc='green')

        cmd_line = ['diamond',
                    'blastp',
                    '-q', self.query_fasta,
                    '-d', self.target_fasta,
                    '-p', self.num_threads,
                    '--outfmt', *self.outfmt.split()]

        if self.max_target_seqs:
            cmd_line.extend(['--max-target-seqs', self.max_target_seqs])

        if self.min_pct_id:
            cmd_line.extend(['--id', self.min_pct_id])

        if self.evalue:
            cmd_line.extend(['--evalue', self.evalue])

        if self.additional_params_for_blastp:
            cmd_line.extend(self.additional_params_for_blastp.split())

        self.run.info('Command line', ' '.join([str(p) for p in cmd_line]), quiet=(not anvio.DEBUG))


        shell_cmd_line=' '.join(str(x) for x in cmd_line)
        out_bytes, ret_code = utils.get_command_output_from_shell(shell_cmd_line)

        try:
            decode_out=out_bytes.decode("utf-8")
        except:
            decode_out=out_bytes

        return(decode_out)


    def blastp_stdin(self, sequence):
        self.run.warning(None, header="DIAMOND BLASTP STDIN", lc="green")
        self.run.info("Additional params for blastp", self.additional_params_for_blastp, mc='green')

        cmd_line = ['diamond',
                    'blastp',
                    '-d', self.target_fasta,
                    '-p', self.num_threads,
                    '--outfmt', *self.outfmt.split()]

        if self.max_target_seqs:
            cmd_line.extend(['--max-target-seqs', self.max_target_seqs])

        if self.min_pct_id:
            cmd_line.extend(['--id', self.min_pct_id])

        if self.evalue:
            cmd_line.extend(['--evalue', self.evalue])

        if self.additional_params_for_blastp:
            cmd_line.extend(self.additional_params_for_blastp.split())

        self.run.info('DIAMOND blastp stdin cmd', ' '.join([str(p) for p in cmd_line]), quiet=(not anvio.DEBUG))

        output = utils.run_command_STDIN(cmd_line, self.run.log_file_path, '>seq\n%s' % sequence,remove_log_file_if_exists=False)

        self.progress.end()

        self.run.info('Diamond blastp results', '%d lines were returned from STDIN call' % len(output))

        return(output)


    def blastp_stdin_multi(self, multisequence):
        self.run.warning(None, header="DIAMOND BLASTP STDIN MULTI", lc="green")
        self.run.info("Additional params for blastp", self.additional_params_for_blastp, mc='green')

        cmd_line = ['diamond',
                    'blastp',
                    '-d', self.target_fasta,
                    '-p', self.num_threads,
                    '--outfmt', *self.outfmt.split()]

        if self.max_target_seqs:
            cmd_line.extend(['--max-target-seqs', self.max_target_seqs])

        if self.min_pct_id:
            cmd_line.extend(['--id', self.min_pct_id])

        if self.evalue:
            cmd_line.extend(['--evalue', self.evalue])

        if self.additional_params_for_blastp:
            cmd_line.extend(self.additional_params_for_blastp.split())

        self.run.info('Command line', ' '.join([str(p) for p in cmd_line]), quiet=(not anvio.DEBUG))

        self.progress.new('DIAMOND')
        self.progress.update('running blastp (using %d thread(s)) ...' % self.num_threads)

        output = utils.run_command_STDIN(cmd_line, self.run.log_file_path, multisequence,remove_log_file_if_exists=False)

        self.progress.end()

        self.run.info('Diamond blastp results', '%d lines were returned from STDIN call' % len(output))
        return(output)


    def makedb_stdin(self, sequence, output_file_path=None):
        self.run.warning(None, header="DIAMOND MAKEDB STDIN", lc="green")
        self.progress.new('DIAMOND')
        self.progress.update('creating the search database (using %d thread(s)) ...' % self.num_threads)

        cmd_line = ['diamond'
                    'makedb',
                    '-d', output_file_path or self.target_fasta,
                    '-p', self.num_threads]

        utils.run_command_STDIN(cmd_line, self.run.log_file_path,sequence)

        self.progress.end()

        expected_output = utils.run_command_STDIN(cmd_line, self.run.log_file_path,sequence)

        expected_output = (output_file_path or self.target_fasta) + '.dmnd'
        self.check_output(expected_output, 'makedb')

        self.run.info('diamond makedb cmd', ' '.join([str(x) for x in cmd_line]), quiet=True)


    def view(self):
        self.run.warning(None, header="DIAMOND VIEW", lc="green")
        self.progress.new('DIAMOND')
        self.progress.update('generating tabular output (using %d thread(s)) ...' % self.num_threads)

        cmd_line = ['diamond',
                    'view',
                    '-a', self.search_output_path + '.daa',
                    '-o', self.tabular_output_path,
                    '-p', self.num_threads,
                    '--outfmt', *self.outfmt.split()]

        self.run.info('Command line', ' '.join([str(x) for x in cmd_line]), quiet=True)

        utils.run_command(cmd_line, self.run.log_file_path)

        self.check_output(self.tabular_output_path, 'view')

        if self.names_dict:
            # if we are here, this means the self.tabular_output_path contains information only about unique
            # entries. We will expand it here so downstream analyses do not need to pay attention to this
            # detail.
            self.run.info('self.names_dict is found', 'Un-uniqueing the tabular output', quiet=True)
            self.progress.update('Un-uniqueing the tabular output ...')

            try:
                int(self.outfmt)
            except:
                if not self.outfmt.startswith("6 qseqid sseqid"):
                    self.progress.end()
                    raise ConfigError("drivers.diamond :: You can't supply a names_dict when running "
                                      "view(...) unless your outfmt starts with '6 qseqid sseqid'. Update "
                                      "utils.ununique_BLAST_tabular_output to fix this problem. If you're a "
                                      "user, please report this on github.")

            utils.ununique_BLAST_tabular_output(self.tabular_output_path, self.names_dict)

        self.progress.end()

        self.run.info('Diamond %s tabular output file' % ('un-uniqued' if self.names_dict else ''), self.tabular_output_path)


    def view_as_dataframe(self, results_path):
        """Returns a dataframe of the results path"""

        if self.outfmt == '6':
            cols = (
                'qseqid',
                'sseqid',
                'pident',
                'length',
                'mismatch',
                'gapopen',
                'qstart',
                'qend',
                'sstart',
                'send',
                'evalue',
                'bitscore',
            )
        else:
            cols = tuple(self.outfmt.split(" ")[1:])

        return pd.read_csv(results_path, sep='\t', comment='#', names=cols, index_col=False)



# coding: utf-8
"""Interface to Diamond."""

import os
import tempfile

import anvio
import anvio.utils as utils
import anvio.terminal as terminal

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


class Diamond:
    def __init__(self, query_fasta=None, target_fasta=None, run=run, progress=progress, num_threads=1, overwrite_output_destinations=False):
        self.run = run
        self.progress = progress

        self.num_threads = num_threads
        self.overwrite_output_destinations = overwrite_output_destinations

        utils.is_program_exists('diamond')

        self.tmp_dir = tempfile.gettempdir()
        self.evalue = 1e-05
        self.min_pct_id = None
        self.max_target_seqs = 100000

        self.query_fasta = query_fasta
        self.target_fasta = target_fasta

        if not self.target_fasta:
            self.target_fasta = self.query_fasta

        self.search_output_path = 'diamond-search-results'
        self.tabular_output_path = 'diamond-search-results.txt'

        if not self.run.log_file_path:
            self.run.log_file_path = 'diamond-log-file.txt'

        self.sensitive = False

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

        if os.path.exists(self.search_output_path + '.daa') and not force_blastp:
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
            raise ConfigError("Pfft. Something probably went wrong with Diamond's '%s' since one of the expected output files are missing.\
                                Please check the log file here: '%s." % (process, self.run.log_file_path))


    def makedb(self, output_file_path=None):
        self.progress.new('DIAMOND')
        self.progress.update('creating the search database (using %d thread(s)) ...' % self.num_threads)

        cmd_line = ['diamond',
                    'makedb',
                    '--in', self.query_fasta,
                    '-d', output_file_path or self.target_fasta,
                    '-p', self.num_threads]

        utils.run_command(cmd_line, self.run.log_file_path)

        self.progress.end()

        expected_output = (output_file_path or self.target_fasta) + '.dmnd'

        self.run.info('diamond makedb cmd', ' '.join([str(x) for x in cmd_line]), quiet=True)
        self.run.info('Diamond search db', expected_output)


    def blastp(self):
        self.run.info('DIAMOND is set to be', 'Sensitive' if self.sensitive else 'Fast')

        cmd_line = ['diamond',
                    'blastp',
                    '-q', self.query_fasta,
                    '-d', self.target_fasta,
                    '-o', self.tabular_output_path,
                    '-t', self.tmp_dir,
                    '-p', self.num_threads]

        cmd_line.append('--sensitive') if self.sensitive else None

        if self.max_target_seqs:
            cmd_line.extend(['--max-target-seqs', self.max_target_seqs])

        if self.min_pct_id:
            cmd_line.extend(['--id', self.min_pct_id])

        if self.evalue:
            cmd_line.extend(['--evalue', self.evalue])

        self.run.info('DIAMOND blastp cmd', ' '.join([str(p) for p in cmd_line]), quiet=(not anvio.DEBUG))

        self.progress.new('DIAMOND')
        self.progress.update('running blastp (using %d thread(s)) ...' % self.num_threads)

        utils.run_command(cmd_line, self.run.log_file_path)

        self.progress.end()

        self.expected_output = self.search_output_path + '.daa'

        self.run.info('Diamond blastp results', self.expected_output)

    def blastp_stdout(self):

        cmd_line = ['diamond',
                    'blastp',
                    '-q', self.query_fasta,
                    '-d', self.target_fasta,
                    '-p', self.num_threads]

        cmd_line.append('--sensitive') if self.sensitive else None

        if self.max_target_seqs:
            cmd_line.extend(['--max-target-seqs', self.max_target_seqs])

        if self.min_pct_id:
            cmd_line.extend(['--id', self.min_pct_id])

        if self.evalue:
            cmd_line.extend(['--evalue', self.evalue])

        self.run.info('DIAMOND blastp cmd', ' '.join([str(p) for p in cmd_line]), quiet=(not anvio.DEBUG))


        shell_cmd_line=' '.join(str(x) for x in cmd_line)
        out_bytes, ret_code = utils.get_command_output_from_shell(shell_cmd_line)

        try:
            decode_out=out_bytes.decode("utf-8")
        except:
            decode_out=out_bytes

        return(decode_out)




    def blastp_stdin(self, sequence):
        self.run.info('DIAMOND is set to be', 'Sensitive' if self.sensitive else 'Fast')

        cmd_line = ['diamond',
                    'blastp',
                    '-d', self.target_fasta,
                    '-p', self.num_threads]

        cmd_line.append('--sensitive') if self.sensitive else None

        if self.max_target_seqs:
            cmd_line.extend(['--max-target-seqs', self.max_target_seqs])

        if self.min_pct_id:
            cmd_line.extend(['--id', self.min_pct_id])

        if self.evalue:
            cmd_line.extend(['--evalue', self.evalue])



        self.run.info('DIAMOND blastp stdin cmd', ' '.join([str(p) for p in cmd_line]), quiet=(not anvio.DEBUG))

        output = utils.run_command_STDIN(cmd_line, self.run.log_file_path, '>seq\n%s' % sequence,remove_log_file_if_exists=False)

        self.progress.end()

        self.run.info('Diamond blastp results', '%d lines were returned from STDIN call' % len(output))

        return(output)

    def blastp_stdin_multi(self, multisequence):
        self.run.info('DIAMOND is set to be', 'Sensitive' if self.sensitive else 'Fast')

        cmd_line = ['diamond',
                    'blastp',
                    '-d', self.target_fasta,
                    '-p', self.num_threads]

        cmd_line.append('--sensitive') if self.sensitive else None

        if self.max_target_seqs:
            cmd_line.extend(['--max-target-seqs', self.max_target_seqs])

        if self.min_pct_id:
            cmd_line.extend(['--id', self.min_pct_id])

        if self.evalue:
            cmd_line.extend(['--evalue', self.evalue])

        self.run.info('DIAMOND blastp stdin cmd', ' '.join([str(p) for p in cmd_line]), quiet=(not anvio.DEBUG))

        self.progress.new('DIAMOND')
        self.progress.update('running blastp (using %d thread(s)) ...' % self.num_threads)

        output = utils.run_command_STDIN(cmd_line, self.run.log_file_path, multisequence,remove_log_file_if_exists=False)

        self.progress.end()

        self.run.info('Diamond blastp results', '%d lines were returned from STDIN call' % len(output))
        return(output)


    def makedb_stdin(self, sequence, output_file_path=None):
        self.progress.new('DIAMOND')
        self.progress.update('creating the search database (using %d thread(s)) ...' % self.num_threads)

        cmd_line = ['diamond',
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
        self.progress.new('DIAMOND')
        self.progress.update('generating tabular output (using %d thread(s)) ...' % self.num_threads)

        cmd_line = ['diamond',
                    'view',
                    '-a', self.search_output_path + '.daa',
                    '-o', self.tabular_output_path,
                    '-p', self.num_threads]

        self.run.info('diamond view cmd', ' '.join([str(x) for x in cmd_line]), quiet=True)

        utils.run_command(cmd_line, self.run.log_file_path)

        self.check_output(self.tabular_output_path, 'view')

        if self.names_dict:
            self.run.info('self.names_dict is found', 'Un-uniqueing the tabular output', quiet=True)
            self.progress.update('Un-uniqueing the tabular output ...')
            # if we are here, this means the self.tabular_output_path contains information only about unique
            # entries. We will expand it here so downstream analyses do not need to pay attention to this
            # detail.
            utils.ununique_BLAST_tabular_output(self.tabular_output_path, self.names_dict)

        self.progress.end()

        self.run.info('Diamond %stabular output file' % ('un-uniqued' if self.names_dict else ''), self.tabular_output_path)

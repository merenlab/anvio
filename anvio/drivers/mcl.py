# coding: utf-8
"""Interface to MCL."""

import os

import anvio
import anvio.utils as utils
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


class MCL:
    def __init__(self, mcl_input_file_path, run=run, progress=progress, num_threads=1):
        self.run = run
        self.progress = progress

        self.inflation = 2.0
        self.mcl_input_file_path = mcl_input_file_path
        self.num_threads = num_threads

        utils.is_program_exists('mcl')

        self.clusters_file_path = 'mcl-clusters.txt'

        if not self.run.log_file_path:
            self.run.log_file_path = 'mcl-log-file.txt'


    def check_output(self, expected_output, process='diamond'):
        if not os.path.exists(expected_output):
            self.progress.end()
            raise ConfigError, "Pfft. Something probably went wrong with MCL's '%s' since one of the expected output files are missing.\
                                Please check the log file here: '%s." % (process, self.run.log_file_path)


    def get_clusters_dict(self):
        self.cluster()

        clusters_dict = {}

        line_no = 1
        for line in open(self.clusters_file_path).readlines():
            clusters_dict['PC_%08d' % line_no] = line.strip().split('\t')

            line_no += 1

        self.run.info('Number of protein clusters', '%s' % pp(len(clusters_dict)))

        return clusters_dict


    def cluster(self):
        self.run.info('MCL inflation', self.inflation)

        self.progress.new('MCL')
        self.progress.update('clustering (using %d thread(s)) ...' % self.num_threads)

        cmd_line = ['mcl', 
                    self.mcl_input_file_path,
                    '--abc',
                    '-I', self.inflation,
                    '-o', self.clusters_file_path,
                    '-te', self.num_threads]

        self.run.info('mcl cmd', ' '.join(map(lambda x: str(x), cmd_line)), quiet=True)

        utils.run_command(cmd_line, self.run.log_file_path)

        self.progress.end()

        self.check_output(self.clusters_file_path, 'makedb')

        self.run.info('MCL output', self.clusters_file_path)

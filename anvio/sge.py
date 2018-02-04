# -*- coding: utf-8
# pylint: disable=line-too-long
"""Module to submit/track jobs for SUN Grid Engine"""

import os
import time
import glob
import random
import string
import subprocess

import anvio
import anvio.fastalib as u
import anvio.utils as utils
import anvio.filesnpaths as filesnpaths

from anvio.errors import ConfigError
from anvio.terminal import pretty_print as pp


__author__ = "Developers of anvi'o (see AUTHORS.txt)"
__copyright__ = "Copyleft 2015-2018, the Meren Lab (http://merenlab.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "A. Murat Eren"
__email__ = "a.murat.eren@gmail.com"
__status__ = "Development"


QSUB_SCRIPT = """#!/bin/sh
#$ -j y
#$ -o %(log)s
#$ -e %(log)s
#$ -N %(identifier)s
#$ -V

%(command)s"""


class Progress:
    def update(self, str):
        print(str)


class Run:
    def info(self, str_1, str_2):
        print("%s: %s" % (str(str_1), str(str_2)))


class SGE:
    def __init__(self):
        """
        This is a simple class to send jobs to Sun Grid Engine and to merge partial
        results. This runs well for me, but hasn't been really tested well with
        different versions of SGE. An example usage is follows:

            ----8<-----8<-----8<-----8<-----8<-----8<-----8<-----8<-----8<-----8<-----
                import os
                import utils
                from anvio.sge import SGE
                sge = SGE()
                sge.check_sge_binaries()
                sge.input_file_path = ...
                sge.tmp_dir = ...
                sge.input_is_fasta = True
                sge.merged_results_file_path = ...
                sge.binary = ... (full path to the binary file)
                sge.command = 'perl %(binary)s %(part)s'
                sge.wild_card_for_partial_results = "results.01.phymm*part-*.txt"

                try:
                    sge._run()
                except ConfigError, e:
                    print e
                    sys.exit(-1)

                os.deltree(sge.tmp_dir)

                # at this point the merged results file must be found here:
                print '%s' % sge.merged_results_file_path
            ----8<-----8<-----8<-----8<-----8<-----8<-----8<-----8<-----8<-----8<-----

        So that's that. Here are some critical variables before calling sge._run:

            * `self.tmp_dir` is the directory to store all parts.
            * `sge.command` is the command template that has to have %(binary) and %(part)s;
               here is an example:

                    sge.command = 'perl %(binary)s %(part)s'

            * `sge.wild_card_for_partial_results`. so everything is split into parts, sent to
               the cluster, and followed until all processes are done. output files expected to
               be found in the self.tmp_dir directory, and it is the user's responsibility to
               format sge.command in a proper way to make sure that is happening. this property
               is a wildcard to look for in self.tmp_dir to merge partial results into one
               results file. here is an example:

                    sge.wild_card_for_partial_results = "results.01.phymm*part-*.txt"

            * `self.merged_results_file_path` is the file all partial results will be merged into one
               file for downstream analyses.
        """

        self.input_file_path = None
        self.merged_results_file_path = None
        self.num_entries_per_file = 10
        self.tmp_dir = None
        self.wild_card_for_partial_results = None

        self.progress = Progress()
        self.run = Run()

        self.input_is_fasta = True

        self.binary = None
        self.command = None


    def _run(self):
        self.check_sge_binaries()

        if not self.binary:
            raise ConfigError('A binary has to be declared.')
        if not self.command:
            raise ConfigError('SGE module cannot run without a command.')
        if not self.tmp_dir:
            raise ConfigError('SGE module needs a tmp dir.')

        filesnpaths.is_file_exists(self.input_file_path)
        filesnpaths.is_output_file_writable(self.merged_results_file_path)

        self.run.info('temp_directory', self.tmp_dir)

        parts = self.split_input_file()

        old_workdir = os.getcwd()
        os.chdir(os.path.dirname(self.tmp_dir))
        self.clusterize(parts)

        if self.wild_card_for_partial_results:
            self.merge_partial_results()
        os.chdir(old_workdir)


    def merge_partial_results(self):
        self.progress.update('Partial results file are being concatenated ...')
        files_to_concat = glob.glob(os.path.join(self.tmp_dir, self.wild_card_for_partial_results))
        if not files_to_concat:
            raise ConfigError("Wild card '%s' didn't return any files to concatenate." % self.wild_card_for_partial_results)

        utils.concatenate_files(self.merged_results_file_path, files_to_concat)


    def check_sge_binaries(self):
        filesnpaths.is_program_exists('qsub')
        filesnpaths.is_program_exists('qstat')


    def split_input_file(self):
        parts = []
        next_part = 1
        part_obj = None

        if self.input_is_fasta:
            fasta = u.SequenceSource(self.input_file_path)

            while next(fasta):
                if (fasta.pos - 1) % self.num_entries_per_file == 0:
                    self.progress.update('Creating part: ~ %s' % (pp(next_part)))

                    if part_obj:
                        part_obj.close()

                    file_path = os.path.join(self.tmp_dir, 'part-%08d.fa' % next_part)
                    parts.append(file_path)
                    next_part += 1
                    part_obj = open(file_path, 'w')

                part_obj.write('>%s\n' % fasta.id)
                part_obj.write('%s\n' % fasta.seq)

            if part_obj:
                part_obj.close()

        return parts


    def clusterize(self, parts):
        # create a 8 digits random identifier for cluster jobs:
        identifier = ''.join(random.choice(string.ascii_uppercase) for x in range(10))

        for part in parts:
            command = self.command % {'binary': self.binary, 'part': part}

            # create sh file
            shell_script = part + '.sh'
            open(shell_script, 'w').write(QSUB_SCRIPT % {'log': part + '.log',
                                                         'identifier': identifier,
                                                         'command': command})

            # submit script to cluster
            utils.run_command('qsub %s' % shell_script)


        while True:
            qstat_info = self.get_qstat_info(identifier)
            total_processes = sum(qstat_info.values())
            if total_processes == 0:
                break

            self.progress.update('Qstat Info :: Total Jobs: %s, %s' % (pp(total_processes),
                       ', '.join(['%s: %s' % (x, pp(qstat_info[x])) for x in qstat_info])))

            time.sleep(5)

        return True


    def get_qstat_info(self, job_identifier):
        try:
            proc = subprocess.Popen(['qstat'], stdout=subprocess.PIPE)
        except OSError as e:
            raise ConfigError("qstat command was failed for the following reason: '%s'" % (e))

        qstat_state_codes = {'Pending': ['qw', 'hqw', 'hRwq'],
                             'Running': ['r', 't', 'Rr', 'Rt'],
                             'Suspended': ['s', 'ts', 'S', 'tS', 'T', 'tT'],
                             'Error': ['Eqw', 'Ehqw', 'EhRqw'],
                             'Deleted': ['dr', 'dt', 'dRr', 'dRt', 'ds', 'dS', 'dT', 'dRs', 'dRS', 'dRT']}

        info_dict = {'Pending': 0, 'Running': 0, 'Suspended': 0, 'Error': 0, 'Deleted': 0}
        line_no = 0

        while True:
            line = proc.stdout.readline()

            # skip the first two lines
            if line_no < 2:
                line_no += 1
                continue

            if line != '':
                id, priority, name, user, state = line.strip().split()[0:5]
                if name == job_identifier:
                    found = False
                    for s in qstat_state_codes:
                        if state in qstat_state_codes[s]:
                            found = True
                            info_dict[s] += 1
                    if not found:
                        raise ConfigError("Unknown state for qstat: '%s' (known states: '%s')"\
                                 % (state, ', '.join(list(info_dict.keys()))))

                line_no += 1
            else:
                break

        return info_dict


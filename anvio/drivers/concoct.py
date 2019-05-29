# -*- coding: utf-8
# pylint: disable=line-too-long
"""Anvi'o - CONCOCT interface for unsupervised clustering

   There are two classes in this file, `CONCOCT` and `CONCOCT_INTERFACE`.
   `CONCOCT` is pretty straightforward to use using anvi'o resources. Here
   is an example:

   >>> import anvio.concoct as concoct
   >>> class Args:
          pass
   >>> args = Args()
   >>> args.profile_db = '/path/to/PROFILE.db'
   >>> args.contigs_db = '/path/to/CONTIGS.db'
   >>> c = concoct.CONCOCT(args)
   >>> c.cluster()
   >>> print c.clusters

   The other class `CONCOCT_INTERFACE`, handles more low-level access
   to the vbgmm module.
"""
import os

import anvio
import anvio.utils as utils
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths

from anvio.errors import ConfigError


__author__ = "Developers of anvi'o (see AUTHORS.txt)"
__copyright__ = "Copyleft 2015-2018, the Meren Lab (http://merenlab.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "A. Murat Eren"
__email__ = "a.murat.eren@gmail.com"
__status__ = "Development"


run = terminal.Run()
progress = terminal.Progress()


class CONCOCT:
    arguments = {
        'seed': (
                ['--seed'],
                {'metavar': "INT",
                 'required': False,
                 'help': "Seed for random numbers"}
                    ),
        'threads': (
                ['-T', '--threads'],
                {'metavar': "INT",
                 'required': False,
                 'help': "Number of threads"}
                    ),
    }

    def __init__(self, run=run, progress=progress):
        self.run = run
        self.progress = progress

        self.program_name = 'concoct'

        utils.is_program_exists(self.program_name)
        self.temp_path = filesnpaths.get_temp_directory_path()



    def cluster(self, input_files, args, threads=1):
        cwd_backup = os.getcwd()
        os.chdir(self.temp_path)
        log_path = os.path.join(self.temp_path, 'logs.txt')
        
        if anvio.DEBUG:
            self.run.info('Working directory', self.temp_path)


        cmd_line = [self.program_name,
            '--coverage_file', input_files.coverage, 
            '--composition_file', input_files.fasta,
            '--threads', threads]

        self.progress.new(self.program_name)
        self.progress.update('Running using %d threads...' % threads)
        utils.run_command(cmd_line, log_path)
        self.progress.end()

        clusters = {}
        with open(os.path.join(self.temp_path, 'clustering_gt1000.csv'), 'r') as f:
            lines = f.readlines()[1:]

            for entry in lines:
                contig, bin_name = map(str.strip, entry.split(','))

                pretty_bin_name = 'Bin_' + bin_name

                if pretty_bin_name not in clusters:
                    clusters[pretty_bin_name] = []

                clusters[pretty_bin_name].append(contig)

        # restore cwd
        os.chdir(cwd_backup)

        return clusters




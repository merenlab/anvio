# coding: utf-8
"""Interface to MetaBAT."""

from subprocess import Popen, PIPE

import anvio
import anvio.utils as utils
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths


__author__ = "Developers of anvi'o (see AUTHORS.txt)"
__copyright__ = "Copyleft 2015-2019, the Meren Lab (http://merenlab.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "Özcan Esen"
__email__ = "ozcanesen@gmail.com"


run = terminal.Run()
progress = terminal.Progress()
pp = terminal.pretty_print


class MetaBAT:
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

    def __init__(self, run=run):
        self.run = run
        self.progress = progress

        utils.is_program_exists('metabat2')

        self.temp_path = filesnpaths.get_temp_directory_path()


    def get_parser(self):
        pass


    def run_command(self, args):
        utils.is_profile_db_and_contigs_db_compatible(args.profile_db, args.contigs_db)

        merged_profile_db = dbops.ProfileDatabase(args.profile_db)

        if(merged_profile_db.meta['merged'] != True):
            raise ConfigError("'%s' does not seem to be a merged profile database :/" % args.profile_db)

        
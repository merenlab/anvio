# -*- coding: utf-8
# pylint: disable=line-too-long
"""
"""
import os

import anvio
import anvio.utils as utils
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths
import anvio.ccollections as ccollections

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


class DAS_Tool:
    arguments = {
        'source-collections': (
                ['-S', '--source-collections'],
                {'metavar': "COLLECTION_LIST",
                 'required': True,
                 'help': "Comma-separated list of collections, case sensitive."}
                    ),
    }

    def __init__(self, run=run, progress=progress):
        self.run = run
        self.progress = progress

        self.program_name = 'DAS_Tool'

        utils.is_program_exists(self.program_name)
        self.temp_path = filesnpaths.get_temp_directory_path()



    def cluster(self, input_files, args, threads=1):
        P = lambda x: os.path.join(self.temp_path, x)

        cwd_backup = os.getcwd()
        os.chdir(self.temp_path)
        log_path = P('logs.txt')
        
        if anvio.DEBUG:
            self.run.info('Working directory', self.temp_path)


        c = ccollections.Collections(r = run, p = progress)
        c.populate_collections_dict(input_files.profile_db)

        source_collections = set(map(str.strip, args.source_collections.split(',')))

        missing_collections = source_collections - set(c.collections_dict.keys())

        if len(missing_collections):
            raise ConfigError("Some of the collections you wanted are missing in the database.\
                              Here is the list of missing collections: %s" % (", ".join(missing_collections)))

        c_names = []
        c_files = []

        for collection_name in source_collections:
            prefix = P(collection_name)

            c_names.append(collection_name)
            c_files.append(prefix + '.txt')

            c.export_collection(collection_name, output_file_prefix=prefix, include_unbinned=False)

        cmd_line = [self.program_name,
            '-c', input_files.fasta,
            '-i', ','.join(c_files),
            '-l', ','.join(c_names),
            '-o', P('OUTPUT'),
            '--threads', str(threads)]

        self.progress.new(self.program_name)
        self.progress.update('Running using %d threads...' % threads)
        utils.run_command(cmd_line, log_path)
        self.progress.end()


        clusters = {}
        with open(P('OUTPUT_DASTool_scaffolds2bin.txt'), 'r') as f:
            lines = f.readlines()

            for entry in lines:
                contig, bin_name = map(str.strip, entry.split())

                pretty_bin_name = 'Bin_' + bin_name

                if pretty_bin_name not in clusters:
                    clusters[pretty_bin_name] = []

                clusters[pretty_bin_name].append(contig)

        # restore cwd
        os.chdir(cwd_backup)

        return clusters

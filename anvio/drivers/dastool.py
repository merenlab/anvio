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
        'search-engine': (
                ['--search-engine'],
                {'metavar': "PROGRAM",
                 'required': False,
                 'default': 'usearch',
                 'help': "Engine used for single copy gene identification [blast/diamond/usearch].\
                              (default: usearch)"}
                    ),
        'score-threshold': (
                ['--score-threshold'],
                {'metavar': "FLOAT",
                 'required': False,
                 'default': 0.5,
                 'help': "Score threshold until selection algorithm will keep selecting bins [0..1].\
                              (default: 0.5)"}
                    ),
        'duplicate-penalty': (
                ['--duplicate-penalty'],
                {'metavar': "FLOAT",
                 'required': False,
                 'default': 0.6,
                 'help': "Penalty for duplicate single copy genes per bin (weight b).\
                              Only change if you know what you're doing. [0..3]\
                              (default: 0.6)"}
                    ),
        'megabin-penalty': (
                ['--megabin-penalty'],
                {'metavar': "FLOAT",
                 'required': False,
                 'default': 0.5,
                 'help': "Penalty for megabins (weight c). Only change if you know what you're doing. [0..3]\
                              (default: 0.5)"}
                    ),
        'db-directory': (
                ['--db-directory'],
                {'metavar': "PATH",
                 'required': False,
                 'default': None,
                 'help': "Directory of single copy gene database. (default: install_dir/db)"}
                    ),
    }
    citation = "Christian M. K. Sieber, Alexander J. Probst, Allison Sharrar, Brian C. Thomas, \
                Matthias Hess, Susannah G. Tringe & Jillian F. Banfield (2018). \
                Recovery of genomes from metagenomes via a dereplication, aggregation \
                and scoring strategy. Nature Microbiology. https://doi.org/10.1038/s41564-018-0171-1."

    def __init__(self, run=run, progress=progress):
        self.run = run
        self.progress = progress

        self.program_name = 'DAS_Tool'

        utils.is_program_exists(self.program_name)
        self.run.info_single("If you publish results from this workflow, \
                               please do not forget to cite \n%s" % DAS_Tool.citation,
                               nl_before=1, nl_after=1, mc='green')
        self.temp_path = filesnpaths.get_temp_directory_path()



    def cluster(self, input_files, args, threads=1, splits_mode=False):
        if not splits_mode:
            raise ConfigError("DAS_Tool can only be run in splits mode. See --help for details.")

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
            '--threads', str(threads),
            *utils.serialize_args(args, 
                use_underscore=True, 
                skip_keys=['source_collections'])]

        self.progress.new(self.program_name)
        self.progress.update('Running using %d threads...' % threads)
        utils.run_command(cmd_line, log_path)
        self.progress.end()


        clusters = {}
        with open(P('OUTPUT_DASTool_scaffolds2bin.txt'), 'r') as f:
            lines = f.readlines()

            for entry in lines:
                contig, bin_name = map(str.strip, entry.split())

                pretty_bin_name = 'Bin_' + bin_name.replace('.', '_')

                if pretty_bin_name not in clusters:
                    clusters[pretty_bin_name] = []

                clusters[pretty_bin_name].append(contig)

        # restore cwd
        os.chdir(cwd_backup)

        return clusters

# -*- coding: utf-8
"""
    Classes for pan operations.

    anvi-pan-genome is the default client of this module
"""

import os

import anvio
import anvio.tables as t
import anvio.utils as utils
import anvio.dbops as dbops
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths

from anvio.errors import ConfigError, FilesNPathsError


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


class Pangenome:
    def __init__(self, args = None, run = run, progress = progress):
        self.run = run
        self.progress = progress

        A = lambda x: args.__dict__[x] if args.__dict__.has_key(x) else None
        input_file_for_contig_dbs = A('input_contig_dbs')
        self.num_CPUs = A('num_CPUs')
        self.output_dir = A('output_dir')
        self.overwrite_output_destinations = A('overwrite_output_destinations')
        self.debug = A('debug')

        self.temp_files_to_remove_later = []

        self.contig_dbs = utils.get_TAB_delimited_file_as_dictionary(input_file_for_contig_dbs, expected_fields = ['name', 'path']) if input_file_for_contig_dbs else {}

        # convert relative paths to absolute paths
        for contigs_db in self.contig_dbs:
            path = self.contig_dbs[contigs_db]['path']
            if not path.startswith('/'):
                self.contig_dbs[contigs_db]['path'] = os.path.abspath(os.path.join(os.path.dirname(input_file_for_contig_dbs), path))


    def get_output_file_path(self, file_name, temp_file = False):
        output_file_path = os.path.join(self.output_dir, file_name)

        if temp_file:
            self.temp_files_to_remove_later.append(output_file_path)

        return output_file_path


    def remove_temp_files(self):
        if self.debug:
            return

        for file_path in self.temp_files_to_remove_later:
            os.remove(file_path)


    def check_programs(self):
        utils.is_program_exists('diamond')


    def check_params(self):
        # deal with the output directory:
        try:
            filesnpaths.is_file_exists(self.output_dir)
        except FilesNPathsError:
            filesnpaths.gen_output_directory(self.output_dir, delete_if_exists = self.overwrite_output_destinations)

        filesnpaths.is_output_dir_writable(self.output_dir)
        self.output_dir = os.path.abspath(self.output_dir)


    def init_contig_dbs(self):
        if type(self.contig_dbs) != type({}):
            raise ConfigError, "self.contig_dbs must be of type dict. Anvi'o needs an adult :("

        if not len(self.contig_dbs):
            raise ConfigError, "There is no contig databases to process..."

        if len(self.contig_dbs) < 2:
            raise ConfigError, "There must be at least two contigs databases for this to work :/"

        if len([c for c in self.contig_dbs.values() if 'path' not in c]):
            raise ConfigError, "self.contig_dbs does not seem to be a properly formatted dictionary for\
                                the anvi'o class Pangenome. You did something very wrong."

        missing_dbs = [c['path'] for c in self.contig_dbs.values() if not os.path.exists(c['path'])]
        if len(missing_dbs):
            raise ConfigError, "%d of %d of your contigs databases are not found where they were supposed to be \
                                based on the description you provided :( Here is one that is missing: '%s'" \
                                                % (len(missing_dbs), len(self.contig_dbs), missing_dbs[0])

        # just go over the contig dbs to make sure they all are OK.
        for contigs_db_name in self.contig_dbs:
            c = self.contig_dbs[contigs_db_name]
            c['name'] = contigs_db_name
            contigs_db = dbops.ContigsDatabase(c['path'])
            for key in contigs_db.meta:
                c[key] = contigs_db.meta[key]
            contigs_db.disconnect()

        # if two contigs db has the same hash, we are kinda f'd:
        if len(set([c['contigs_db_hash'] for c in self.contig_dbs.values()])) != len(self.contig_dbs):
            raise ConfigError, 'Not all hash values are unique across all contig databases you provided. Something\
                                very fishy is going on :/'

        # make sure genes are called in every contigs db:
        if len([c['genes_are_called'] for c in self.contig_dbs.values()]) != len(self.contig_dbs):
            raise ConfigError, 'Genes are not called in every contigs db in the collection :/'

        self.run.info('Contig DBs', '%d contig databases have been found.' % len(self.contig_dbs))


    def gen_combined_proteins_fasta(self):
        self.progress.new('Storing combined protein sequences')
        output_file_path = self.get_output_file_path('combined_proteins.fa', temp_file = True)
        output_file = open(output_file_path, 'w')

        for c in self.contig_dbs.values():
            self.progress.update('Working on %s ...' % c['name'])
            num_genes = 0
            contigs_db = dbops.ContigsDatabase(c['path'])
            protein_sequences = contigs_db.db.get_table_as_dict(t.gene_protein_sequences_table_name)
            for protein_id in protein_sequences:
                num_genes += 1
                output_file.write('>%s_%d\n' % (c['contigs_db_hash'], protein_id))
                output_file.write('%s\n' % protein_sequences[protein_id]['sequence'])
            contigs_db.disconnect()
            c['num_genes'] = num_genes

        output_file.close()
        self.progress.end()

        self.run.info('ORFs', '%s protein sequences are stored for analysis.' % pp(sum([c['num_genes'] for c in self.contig_dbs.values()])))

        return output_file_path


    def run_diamond(self, combined_proteins_fasta_path):
        diamond = Diamond(run = self.run, progress = self.progress)

        log_file_path = self.get_output_file_path('log.txt')
        db_path = self.get_output_file_path('.'.join(combined_proteins_fasta_path.split('.')[:-1]))

        force_makedb, force_blastp, force_view = False, False, False

        if self.overwrite_output_destinations:
            force_makedb = True

        if os.path.exists(db_path + '.dmnd') and not force_makedb:
            run.info_single("Notice: A diamond database is found in the output directory, and will be used!", mc = 'red', nl_before = 1)
        else:
            diamond.makedb(combined_proteins_fasta_path, db_path, log_file_path)
            force_blastp, forrce_view = True, True

        search_output_path = self.get_output_file_path('search_results')
        if os.path.exists(search_output_path + '.daa') and not force_blastp:
            run.info_single("Notice: A DIAMOND search result is found in the output directory: skipping BLASTP!", mc = 'red', nl_before = 1)
        else:
            diamond.blastp(combined_proteins_fasta_path, db_path, search_output_path, self.output_dir, log_file_path)
            force_view = True

        tabular_output_path = self.get_output_file_path('search_results.txt')
        if os.path.exists(tabular_output_path) and not force_view:
            run.info_single("Notice: A DIAMOND tabular output is found in the output directory. Anvi'o will not generate another one!", mc = 'red', nl_before = 1)
        else:
            diamond.view(search_output_path, tabular_output_path, log_file_path)


    def sanity_check(self):
        self.check_programs()
        self.check_params()
        self.init_contig_dbs()


    def process(self):
        self.sanity_check()

        # first we will export all proteins 
        combined_proteins_fasta_path = self.gen_combined_proteins_fasta()

        # run diamond
        self.run_diamond(combined_proteins_fasta_path)


class Diamond:
    def __init__(self, run = run, progress = progress):
        self.run = run
        self.progress = progress

    def makedb(self, combined_proteins_fasta_path, db_path, log_file_path):
        self.progress.new('DIAMOND')
        self.progress.update('creating the search database ...')
        cmd_line = ('diamond makedb --in %s -d %s >> "%s" 2>&1' % (combined_proteins_fasta_path,
                                                                   db_path,
                                                                   log_file_path))

        with open(log_file_path, "a") as log: log.write('CMD: ' + cmd_line + '\n')

        utils.run_command(cmd_line)

        self.progress.end()
        self.run.info('Diamond temp search db', db_path + '.dmnd')


    def blastp(self, combined_proteins_fasta_path, db_path, search_output_path, output_dir, log_file_path):
        self.progress.new('DIAMOND')
        self.progress.update('running blastp ...')
        cmd_line = ('diamond blastp -q %s -d %s -a %s -t %s >> "%s" 2>&1' % (combined_proteins_fasta_path,
                                                                             db_path,
                                                                             search_output_path,
                                                                             output_dir,
                                                                             log_file_path))
        with open(log_file_path, "a") as log: log.write('CMD: ' + cmd_line + '\n')

        utils.run_command(cmd_line)

        self.progress.end()
        self.run.info('Diamond blastp results', search_output_path + '.daa')


    def view(self, search_output_path, tabular_output_path, log_file_path):
        self.progress.new('DIAMOND')
        self.progress.update('creating a tabular output file ...')
        cmd_line = ('diamond view -a %s -o %s >> "%s" 2>&1' % (search_output_path + '.daa',
                                                               tabular_output_path,
                                                               log_file_path))
        with open(log_file_path, "a") as log: log.write('CMD: ' + cmd_line + '\n')

        utils.run_command(cmd_line)

        self.progress.end()
        self.run.info('Diamond tabular output file', tabular_output_path)

# -*- coding: utf-8
# pylint: disable=line-too-long
"""Module to deal with HDF5 files"""

import h5py
import numpy as np

import anvio
import anvio.tables as t
import anvio.utils as utils
import anvio.terminal as terminal
import anvio.constants as constants
import anvio.filesnpaths as filesnpaths

from anvio.errors import HDF5Error


__author__ = "A. Murat Eren"
__copyright__ = "Copyright 2015, The anvio Project"
__credits__ = []
__license__ = "GPL 3.0"
__maintainer__ = "A. Murat Eren"
__email__ = "a.murat.eren@gmail.com"


run = terminal.Run()
progress = terminal.Progress()
pp = terminal.pretty_print


class HDF5_IO(object):
    def __init__(self, file_path, unique_hash, create_new=False, ignore_hash=False):
        self.file_path = file_path

        if create_new:
            if ignore_hash:
                raise HDF5Error, "When creating a new database, you can't use the 'ignore_hash'\
                                  parameter."

            if not unique_hash:
                raise HDF5Error, "When creating a new database, the 'unique_hash' cannot be None"

            self.fp = h5py.File(self.file_path, 'w')
            self.fp.attrs['hash'] = unique_hash
            self.fp.attrs['version'] = self.version
        else:
            filesnpaths.is_file_exists(self.file_path)
            self.fp = h5py.File(self.file_path, 'r')

            if self.fp.attrs['version'] != self.version:
                raise HDF5Error, "The data file for %s ('%s') is at version '%s', however, your client is at\
                                  version '%s'. This is bad news, because your version of anvi'o can't work with\
                                  this file. You can regenerate the data file using the current version of anvi'o,\
                                  or look around to see whether there is an upgrade script is available (a good start\
                                  would be to type 'anvi-script-upgrade-' and then click TAB key twice). Otherwise you\
                                  may want to consider sending an e-mail to the anvi'o developers to find out what's up.\
                                  We heard that they love them some e-mails." % (self.db_type, self.file_path, self.fp.attrs['version'], self.version)

            if not ignore_hash and self.fp.attrs['hash'] != unique_hash:
                raise HDF5Error, "The database at '%s' does not seem to be compatible with the client :/\
                                  (i.e., the hash values do not match)." % self.file_path

            self.unique_hash = self.fp.attrs['hash']


    def add_integer_list(self, path, l, data_type='uint16'):
        """Add an array into the the HDF5 file.
        
            >>> h = HDF5_IO('test.h5')
            >>> l = [1, 2, 3, 4, 5]
            >>> h.add_integer_list('/split_1/sample_x', l)
            >>> h.close()
        """

        new_data_obj = self.fp.create_dataset(path, (len(l),), dtype=np.dtype(data_type))
        new_data_obj[...] = np.array(l)


    def get_integer_list(self, path):
        l = self.fp[path]
        return l.value


    def path_exists(self, path):
        return path in self.fp


    def close(self):
        self.fp.close()


class AuxiliaryDataForSplitCoverages(HDF5_IO):
    """A class to handle HDF5 operations to store and access split coverages"""
    def __init__(self, file_path, db_hash, create_new = False, ignore_hash = False, run=run, progress=progress, quiet = False):
        self.db_type = 'auxiliary data for coverages'
        self.version = anvio.__hdf5__version__

        HDF5_IO.__init__(self, file_path, db_hash, create_new = create_new, ignore_hash = ignore_hash)

        self.quiet = quiet


    def is_known_split(self, split_name):
        if not self.path_exists('/data/coverages/%s' % split_name):
            raise HDF5Error, 'The database at "%s" does not know anything about "%s" :(' % (self.file_path, split_name)


    def append(self, split_name, sample_id, coverage_list):
        self.add_integer_list('/data/coverages/%s/%s' % (split_name, sample_id), coverage_list)


    def get(self, split_name):
        self.is_known_split(split_name)

        d = {}

        sample_names = self.fp['/data/coverages/%s' % split_name].keys()

        for sample_name in sample_names:
            d[sample_name] = self.get_integer_list('/data/coverages/%s/%s' % (split_name, sample_name))

        return d


class AuxiliaryDataForNtPositions(HDF5_IO):
    """A class to handle HDF5 operations to store and access split coverages"""
    def __init__(self, file_path, db_hash, create_new = False, run=run, progress=progress, quiet = False):
        self.db_type = 'auxiliary data for nt positions'
        self.version = anvio.__hdf5__version__

        HDF5_IO.__init__(self, file_path, db_hash, create_new = create_new)

        self.quiet = quiet


    def is_known_contig(self, contig_name):
        path = '/data/nt_position_info/%s' % contig_name
        return self.path_exists(path)


    def append(self, contig_name, position_info_list):
        self.add_integer_list('/data/nt_position_info/%s' % contig_name, position_info_list, data_type = 'uint8')


    def get(self, contig_name):
        if not self.is_known_contig(contig_name):
            return []

        return self.get_integer_list('/data/nt_position_info/%s' % contig_name)


class GenomesDataStorage(HDF5_IO):
    """A class to handle HDF5 operations to store and access protein sequnces in pan genome analyses.
    
       An example:

           >>> x = a.GenomesDataStorage('test.h5', 'unique_hash', create_new=True)
           >>> x.add_gene_call_data('genome_name', int_gene_caller_id,
                                    sequence = 'IMLQWIVIIYFLVINLVLFSMMGYDKKQAKRGNWRIPERRLLTIGLVGGGLGGLMGQKKFHHKTQKPVFALCYSIGVIAMISCIYLTFK',
                                    partial = 0,
                                    functions = [('PFAM', 'PFAM FUNC_1'), ('TIGRFAM', 'sik')],
                                    taxonomy_dict = {'t_phylum': 'phy', 't_class': 'clas', 't_order': 'ord', 't_family': None, 't_genus': 'genus', 't_species': 'sp'})
           >>> x.close()
           >>> x = a.GenomesDataStorage('test.h5', 'unique_hash')
           >>> x.get_gene_sequence('genome_name', int_gene_caller_id)
           IMLQWIVIIYFLVINLVLFSMMGYDKKQAKRGNWRIPERRLLTIGLVGGGLGGLMGQKKFHHKTQKPVFALCYSIGVIAMISCIYLTFK
    """

    def __init__(self, file_path, db_hash, create_new=False, ignore_hash=False, run=run, progress=progress, quiet=False):
        self.version = anvio.__genomes_storage_version__
        self.db_type = 'genomes data storage'

        HDF5_IO.__init__(self, file_path, db_hash, create_new = create_new, ignore_hash = ignore_hash)

        self.run = run
        self.progress = progress
        self.quiet = quiet

        self.essential_genome_info = constants.essential_genome_info + ['genome_hash', 'external_genome']

        if not create_new:
            self.num_genomes = len(self.fp['/info/genomes'])

        self.D = lambda genome_name: self.fp['/data/genomes/%s' % genome_name]
        self.G = lambda gene_callers_id, genome_data: genome_data['%d' % gene_callers_id]


    def is_known_genome(self, genome_name, throw_exception=True):
        if not self.path_exists('/info/genomes/%s' % genome_name):
            if throw_exception:
                raise HDF5Error, 'The database at "%s" does not know anything about "%s" :(' % (self.file_path, genome_name)
            else:
                return False


    def is_known_gene_call(self, genome_name, gene_caller_id):
        if not self.path_exists('/data/genomes/%s/%d' % (genome_name, gene_caller_id)):
            raise HDF5Error, 'The genome "%s" does not know anything about the gene caller id "%d" :(' % (genome_name, gene_caller_id)


    def add_gene_call_data(self, genome_name, gene_caller_id, sequence, partial=0, functions = [], taxonomy_dict = None):
        """Add a gene call in a genome into the database"""
        self.fp['/data/genomes/%s/%d/sequence' % (genome_name, gene_caller_id)] = sequence
        self.fp['/data/genomes/%s/%d/length' % (genome_name, gene_caller_id)] = len(sequence)
        self.fp['/data/genomes/%s/%d/partial' % (genome_name, gene_caller_id)] = partial

        if taxonomy_dict:
            for t_level in t.taxon_names_table_structure[1:]:
                self.fp['/data/genomes/%s/%d/taxonomy/%s' % (genome_name, int(gene_caller_id), t_level)] = taxonomy_dict[t_level] or ''

        for source, function in functions:
            self.fp['/data/genomes/%s/%d/functions/%s' % (genome_name, gene_caller_id, source)] = function


    def is_partial_gene_call(self, genome_name, gene_caller_id):
        self.is_known_genome(genome_name)
        self.is_known_gene_call(genome_name, gene_caller_id)

        d = self.fp['/data/genomes/%s/%d/partial' % (genome_name, gene_caller_id)]

        return d.value


    def get_gene_sequence(self, genome_name, gene_caller_id):
        self.is_known_genome(genome_name)
        self.is_known_gene_call(genome_name, gene_caller_id)

        d = self.fp['/data/genomes/%s/%d/sequence' % (genome_name, gene_caller_id)]

        return d.value


    def add_genome(self, genome_name, info_dict):
        if self.is_known_genome(genome_name, throw_exception=False):
            raise "Genome '%s' is already in this data storage :/" % genome_name

        for key in self.essential_genome_info:
            self.fp['/info/genomes/%s/%s' % (genome_name, key)] = info_dict[key]


    def get_storage_hash(self):
        return self.fp.attrs['hash']


    def get_genomes_dict(self):
        genomes_dict = {}

        for d in self.fp['/info/genomes']:
            genomes_dict[d] = {}

        for genome_name in genomes_dict:
            # add every key-value pair we know of in to the dict:
            for key in self.fp['/info/genomes/%s' % genome_name]:
                genomes_dict[genome_name][key] = self.fp['/info/genomes/%s/%s' % (genome_name, key)].value

            # add in AA sequence lengths dict for each gene caller uisng '/data/genomes':
            genomes_dict[genome_name]['gene_lengths'] = {}
            genome_data = self.D(genome_name)
            for gene_caller_id in genome_data:
                genomes_dict[genome_name]['gene_lengths'][int(gene_caller_id)] = genome_data[gene_caller_id]['length'].value

        return genomes_dict


    def gen_combined_protein_sequences_FASTA(self, output_file_path, exclude_partial_gene_calls=False):
        self.run.info('Exclude partial gene calls', exclude_partial_gene_calls, nl_after=1)

        genomes = self.get_genomes_dict()

        total_num_protein_sequences = 0
        total_num_excluded_protein_sequences = 0

        output_file = open(output_file_path, 'w')

        for genome_name in genomes:
            self.progress.new('Storing protein sequences')
            self.progress.update('%s ...' % genome_name)

            genome_data = self.D(genome_name)
            gene_caller_ids = sorted([int(i[0]) for i in genome_data.items()])

            for gene_caller_id in gene_caller_ids:
                partial = self.G(gene_caller_id, genome_data)['partial'].value

                if exclude_partial_gene_calls and partial:
                    total_num_excluded_protein_sequences += 1
                    continue

                sequence = self.G(gene_caller_id, genome_data)['sequence'].value

                output_file.write('>%s_%d\n' % (genomes[genome_name]['genome_hash'], int(gene_caller_id)))
                output_file.write('%s\n' % sequence)

                total_num_protein_sequences += 1

            self.progress.end()

        self.progress.new('Uniquing the output FASTA file')
        self.progress.update('...')
        unique_proteins_FASTA_path, unique_proteins_names_file_path, unique_proteins_names_dict = utils.unique_FASTA_file(output_file_path, store_frequencies_in_deflines=False)
        self.progress.end()

        self.run.info('Unique protein sequences FASTA', output_file_path)
        self.run.info('Num protein sequences reported', '%s' % pp(total_num_protein_sequences), nl_before=1)
        self.run.info('Num excluded gene calls', '%s' % pp(total_num_excluded_protein_sequences))

        return unique_proteins_FASTA_path, unique_proteins_names_dict

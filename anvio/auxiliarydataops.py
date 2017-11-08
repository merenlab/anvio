# -*- coding: utf-8
# pylint: disable=line-too-long
"""Module to deal with HDF5 files"""

import time
import numpy as np

import anvio
import anvio.db as db
import anvio.tables as t
import anvio.utils as utils
import anvio.terminal as terminal
import anvio.constants as constants
import anvio.filesnpaths as filesnpaths

from anvio.errors import HDF5Error, AuxiliaryDataError


__author__ = "A. Murat Eren"
__copyright__ = "Copyright 2015, The anvio Project"
__credits__ = []
__license__ = "GPL 3.0"
__maintainer__ = "A. Murat Eren"
__email__ = "a.murat.eren@gmail.com"


run = terminal.Run()
progress = terminal.Progress()
pp = terminal.pretty_print

class AuxiliaryDataForSplitCoverages(object):
    def __init__(self, file_path, db_hash, create_new=False, ignore_hash=False, run=run, progress=progress, quiet=False):
        self.db_type = 'auxiliary data for coverages'
        self.db_hash = db_hash
        self.version = anvio.__auxiliary_data_version__
        self.file_path = file_path
        self.quiet = quiet
        self.run = run
        self.progress = progress
        self.numpy_data_type = 'uint16'

        self.db = db.DB(self.file_path, self.version, new_database=create_new)

        if create_new:
            self.create_tables()

        if not ignore_hash:
            self.check_hash()


    def create_tables(self):
        self.db.set_meta_value('db_type', self.db_type)
        self.db.set_meta_value('contigs_db_hash', self.db_hash)
        self.db.set_meta_value('creation_date', time.time())

        self.db.create_table(t.split_coverages_table_name, t.split_coverages_table_structure, t.split_coverages_table_types)


    def check_hash(self):
        actual_db_hash = self.db.get_meta_value('contigs_db_hash')
        if self.db_hash != actual_db_hash:
            raise AuxiliaryDataError('The hash value inside Auxiliary Database "%s" does not match with Contigs Database hash "%s",\
                                      this files probaby belong to different projects.' % (actual_db_hash, self.db_hash))


    def append(self, split_name, sample_name, coverage_list):
        coverage_list_blob = utils.convert_numpy_array_to_binary_blob(np.array(coverage_list, dtype=self.numpy_data_type))
        self.db.insert(t.split_coverages_table_name, values=(split_name, sample_name, coverage_list_blob, ))


    def get_all_known_split_names(self):
        return set(self.db.get_single_column_from_table(t.split_coverages_table_name, 'split_name'))


    def get_all(self):
        return self.get_coverage_for_multiple_splits(self.get_all_known_split_names())


    def get_coverage_for_multiple_splits(self, split_names):
        self.progress.new('Recovering split coverages')
        self.progress.update('...')

        split_coverages = {}
        all_known_splits = self.get_all_known_split_names()

        for split_name in split_names:
            self.progress.update('Processing split "%s"' % split_name)

            split_coverages[split_name] = self.get(split_name)

        self.progress.end()
        return split_coverages        


    def get(self, split_name):
        cursor = self.db._exec('''SELECT sample_name, coverages FROM %s WHERE split_name = "%s"''' % 
                                                 (t.split_coverages_table_name, split_name))

        rows = cursor.fetchall()

        if len(rows) == 0:
            raise AuxiliaryDataError('Database does not know anything about split "%s"' % split_name)

        split_coverage = {}
        for row in rows:
            sample_name, coverage_blob = row # unpack sqlite row tuple

            split_coverage[sample_name] = utils.convert_binary_blob_to_numpy_array(coverage_blob, dtype=self.numpy_data_type).tolist()
        
        return split_coverage


    def close(self):
        self.db.disconnect()


class AuxiliaryDataForNtPositions(object):
    def __init__(self, file_path, db_hash, create_new=False, ignore_hash=False, progress=progress, run=run):
        self.file_path = file_path
        self.db_type = 'auxiliary data for nt positions'
        self.version = anvio.__auxiliary_data_version__
        self.db_hash = db_hash
        self.nt_position_info = {}
        self.numpy_data_type = 'uint8'

        self.db = db.DB(self.file_path, self.version, new_database=create_new)

        if create_new:
            self.create_tables()

        if not ignore_hash:
            self.check_hash()

        self.load_all_position_info()


    def create_tables(self):
        self.db.set_meta_value('db_type', self.db_type)
        self.db.set_meta_value('contigs_db_hash', self.db_hash)
        self.db.set_meta_value('creation_date', time.time())

        self.db.create_table(t.nt_position_info_table_name, t.nt_position_info_table_structure, t.nt_position_info_table_types)


    def check_hash(self):
        actual_db_hash = self.db.get_meta_value('contigs_db_hash')
        if self.db_hash != actual_db_hash:
            raise AuxiliaryDataError('The hash value inside Auxiliary Database "%s" does not match with Contigs Database hash "%s",\
                                      this files probaby belong to different projects.' % (actual_db_hash, self.db_hash))


    def is_known_contig(self, contig_name):
        return contig_name in self.nt_position_info


    def load_all_position_info(self):
        cursor = self.db._exec('''SELECT * FROM "%s"''' % (t.nt_position_info_table_name))

        for row in cursor.fetchall():
            contig_name, position_info_blob = row

            self.nt_position_info[contig_name] = utils.convert_binary_blob_to_numpy_array(position_info_blob, dtype=self.numpy_data_type)


    def get_nt_position_info(self, contig_name, pos_in_contig):
        """This function returns a tuple with three items for each nucleotide position.

            (in_partial_gene_call, in_complete_gene_call, base_pos_in_codon)

        See `init_nt_position_info_dict` for more info."""

        if not self.nt_position_info:
            raise AuxiliaryDataError("get_nt_position_info: I am asked to return stuff, but self.nt_position_info is None!\
                                This may happen if you don't have the '.h5' file for your contigs database in the same\
                                directory with your contigs database. But if you do have it there, then anvi'o really\
                                needs an adult :(")

        if not self.nt_positions_info.is_known_contig(contig_name):
            return (0, 0, 0)

        position_info = self.nt_positions_info[contig_name][pos_in_contig]

        if not position_info:
            return (0, 0, 0)
        if position_info == 8:
            return (1, 0, 0)
        if position_info == 4:
            return (0, 1, 1)
        if position_info == 2:
            return (0, 1, 2)
        if position_info == 1:
            return (0, 1, 3)


    def compress_nt_position_info(self, contig_length, genes_in_contig, genes_in_contigs_dict):
        """This function compresses information regarding each nucleotide position in a given contig
           into a small int. Every nucleotide position is represented by four bits depending on whether
           they occur in a complete opoen reading frame, and which base they correspond to in a codon.

                0000
                ||||
                ||| \
                |||   Third codon?
                || \
                ||   Second codon?
                | \
                |   First codon?
                 \
                   Whether the position in a partial gene call

           8: int('1000', 2); nt position is in a partial gene call
           4: int('0100', 2); nt position is in a complete gene call, and is at the 1st position in the codon
           2: int('0010', 2); nt position is in a complete gene call, and is at the 2nd position in the codon
           1: int('0001', 2); nt position is in a complete gene call, and is at the 3rd position in the codon
        """

        # first we create a list of zeros for each position of the contig
        nt_position_info_list = [0] * contig_length

        for gene_unique_id, start, stop in genes_in_contig:
            gene_call = genes_in_contigs_dict[gene_unique_id]

            # if the gene call is a partial one, meaning the gene was cut at the beginning or
            # at the end of the contig, we are not going to try to make sense of synonymous /
            # non-synonmous bases in that. the clients who wish to use these variables must also
            # be careful about the difference
            if gene_call['partial']:
                for nt_position in range(start, stop):
                    nt_position_info_list[nt_position] = 8
                continue

            if gene_call['direction'] == 'f':
                for nt_position in range(start, stop, 3):
                    nt_position_info_list[nt_position] = 4
                    nt_position_info_list[nt_position + 1] = 2
                    nt_position_info_list[nt_position + 2] = 1
            elif gene_call['direction'] == 'r':
                for nt_position in range(stop - 1, start - 1, -3):
                    nt_position_info_list[nt_position] = 4
                    nt_position_info_list[nt_position - 1] = 2
                    nt_position_info_list[nt_position - 2] = 1

        return nt_position_info_list


    def append(self, contig_name, contig_length, genes_in_contig, genes_in_contigs_dict):
        position_info_list = self.compress_nt_position_info(contig_length, genes_in_contig, genes_in_contigs_dict)
        position_info_blob = utils.convert_numpy_array_to_binary_blob(np.array(position_info_list, dtype=self.numpy_data_type))

        self.db.insert(t.nt_position_info_table_name, values=(contig_name, position_info_blob, ))


    def close(self):
        self.db.disconnect()

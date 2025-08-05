# -*- coding: utf-8
# pylint: disable=line-too-long
"""Module to deal with HDF5 files"""

import time
import numpy as np

import anvio
import anvio.db as db
import anvio.tables as t
import anvio.terminal as terminal

from anvio.errors import AuxiliaryDataError
from anvio.utils.algorithms import convert_binary_blob_to_numpy_array, convert_numpy_array_to_binary_blob


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__maintainer__ = "A. Murat Eren"
__email__ = "a.murat.eren@gmail.com"

run = terminal.Run()
progress = terminal.Progress()
pp = terminal.pretty_print

DEFAULT_COVERAGE_DTYPE = 'uint16'
DEFAULT_COVERAGE_MAX_VALUE = np.iinfo(DEFAULT_COVERAGE_DTYPE).max
TRNASEQ_COVERAGE_DTYPE = 'uint32'
TRNASEQ_COVERAGE_MAX_VALUE = np.iinfo(TRNASEQ_COVERAGE_DTYPE).max

class AuxiliaryDataForSplitCoverages(object):
    def __init__(self, db_path, db_hash, db_variant='unknown', create_new=False, ignore_hash=False, run=run, progress=progress, quiet=False):
        self.db_type = 'auxiliary data for coverages'
        self.db_hash = str(db_hash)
        self.db_variant = str(db_variant)
        self.version = anvio.__auxiliary_data_version__
        self.db_path = db_path
        self.quiet = quiet
        self.run = run
        self.progress = progress
        self.coverage_entries = []
        self.create_new = create_new

        self.coverage_dtype = TRNASEQ_COVERAGE_DTYPE if db_variant == 'trnaseq' else DEFAULT_COVERAGE_DTYPE
        self.coverage_max_value = TRNASEQ_COVERAGE_MAX_VALUE if db_variant == 'trnaseq' else DEFAULT_COVERAGE_MAX_VALUE
        self.db = db.DB(self.db_path, self.version, new_database=self.create_new)

        if self.create_new:
            self.create_tables()

        if not ignore_hash:
            self.check_hash()


    def create_tables(self):
        self.db.set_meta_value('db_type', self.db_type)
        self.db.set_meta_value('contigs_db_hash', self.db_hash)
        self.db.set_meta_value('creation_date', time.time())

        self.db.create_table(t.split_coverages_table_name, t.split_coverages_table_structure, t.split_coverages_table_types)


    def check_hash(self):
        actual_db_hash = str(self.db.get_meta_value('contigs_db_hash'))
        if self.db_hash != actual_db_hash:
            raise AuxiliaryDataError('The hash value inside Auxiliary Database "%s" does not match with Contigs Database hash "%s",\
                                      these files probaby belong to different projects.' % (actual_db_hash, self.db_hash))


    def append(self, split_name, sample_name, coverage_array):
        if (coverage_array == 0).all():
            # Oh my. This split has yielded 0 reads, resulting in a coverage array of strictly zeros.
            # Instead of storing a gzipped byte representation of `coverage_array`, we save the
            # space and instead store a simple integer value equal to the length of the split.
            # Later, when downstream applications query the coverage of this split, the dedicated
            # method self.get will know how to interpret this integer value and return the
            # appropriate array of zeros. This saves some storage space, and also makes querying the
            # coverage of this split faster
            stored_coverage = len(coverage_array)
        else:
            stored_coverage = convert_numpy_array_to_binary_blob(coverage_array.astype(self.coverage_dtype))

        self.coverage_entries.append((split_name, sample_name, stored_coverage, ))


    def store(self):
        self.db.insert_many(t.split_coverages_table_name, entries=self.coverage_entries)
        self.coverage_entries = []


    def get_all_known_split_names(self):
        return set(self.db.get_single_column_from_table(t.split_coverages_table_name, 'split_name'))


    def get_all(self, split_names=None):
        if split_names and not isinstance(split_names, set):
            raise AuxiliaryDataError("Split names for auxiliarydataops::get_all must be of type `set`.`")
        else:
            split_names = self.get_all_known_split_names()

        return self.get_coverage_for_multiple_splits(split_names)


    def get_coverage_for_multiple_splits(self, split_names):
        self.progress.new('Recovering split coverages')
        self.progress.update('...')

        split_coverages = {}

        num_split_names = len(split_names)
        counter = 0
        for split_name in split_names:
            if counter % 10 == 0:
                self.progress.update('Processing split %d of %d (%s) ...' % (counter + 1, num_split_names, split_name))

            split_coverages[split_name] = self.get(split_name)
            counter += 1

        self.progress.end()
        return split_coverages


    def get(self, split_name):
        cursor = self.db._exec('''SELECT sample_name, coverages FROM %s WHERE split_name = "%s"''' %
                                                 (t.split_coverages_table_name, split_name))

        rows = cursor.fetchall()

        if len(rows) == 0:
            raise AuxiliaryDataError('The auxiliary database at "%s" does not know anything about the split "%s"' % (self.db_path, split_name))

        split_coverage = {}
        for row in rows:
            sample_name, blob = row # unpack sqlite row tuple

            if isinstance(blob, int):
                # Look what we have here! Most typically, we store split coverage data as a gzipped
                # byte representation of a numpy array, where each element in the numpy array
                # denotes the coverage of a given nucleotide. However, an exception is made whenever
                # a split has a numpy array of strictly zeros. In this case, we instead store `blob`
                # as an integer representing the split length, which is used in the line below to
                # construct an array of zeros on the fly. This saves us the requirement of storing
                # these zeros, but it also makes querying faster, since no array has to be
                # decompressed. So if you've ended up here, you're enjoying a 2X speed gain.
                coverage_array = np.zeros(blob, dtype=self.coverage_dtype)
            else:
                coverage_array = convert_binary_blob_to_numpy_array(blob, dtype=self.coverage_dtype)

            split_coverage[sample_name] = coverage_array

        return split_coverage


    def close(self):
        """Carries out teardown operations for the table

        An index is created on (split_name, sample_name) if self.create_new == True and if one
        doesn't already exist.

        Notes
        =====
        - For fast accession, it is important to have an index for the table. However, INSERT
          operations are very slow for tables with indices. Hence, the index is created in this
          teardown method, rather than when the table is first created.
        """

        if self.create_new:
            self.db._exec("""CREATE INDEX IF NOT EXISTS covering_index ON %s(split_name, sample_name)""" % (t.split_coverages_table_name))

        self.db.disconnect()

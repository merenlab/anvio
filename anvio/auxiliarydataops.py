"""Module to deal with HDF5 files"""

import json
import time
import numpy as np

import anvio
import anvio.db as db
import anvio.tables as t
import anvio.utils as utils
import anvio.terminal as terminal

from anvio.errors import AuxiliaryDataError


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
IDENTITY_HISTOGRAM_DTYPE = 'uint32'

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
        self.identity_histogram_entries = []
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
        self.db.create_table(t.split_identity_histograms_table_name, t.split_identity_histograms_table_structure, t.split_identity_histograms_table_types)

        # becomes True only once a profiler stores identity histograms (see self.set_identity_histogram_meta)
        self.db.set_meta_value('identity_histograms_stored', False)


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
            stored_coverage = utils.convert_numpy_array_to_binary_blob(coverage_array.astype(self.coverage_dtype))

        self.coverage_entries.append((split_name, sample_name, stored_coverage, ))


    def append_identity_histogram(self, split_name, sample_name, num_below_min, histogram):
        """Buffers the read identity histogram of a split in a sample for storage (see self.store)

        Parameters
        ==========
        split_name : str
        sample_name : str
        num_below_min : int
            The number of reads with a percent identity below the first bin.
        histogram : 1D numpy array
            The number of reads in each bin (see bamops.ReadIdentityBins).
        """

        if not histogram.any():
            # just like split coverages (see self.append), a histogram of zeros is stored as an integer,
            # which here is the number of bins
            stored_histogram = len(histogram)
        else:
            stored_histogram = utils.convert_numpy_array_to_binary_blob(histogram.astype(IDENTITY_HISTOGRAM_DTYPE))

        self.identity_histogram_entries.append((split_name, sample_name, int(num_below_min), stored_histogram, ))


    def store(self):
        self.db.insert_many(t.split_coverages_table_name, entries=self.coverage_entries)
        self.coverage_entries = []

        self.db.insert_many(t.split_identity_histograms_table_name, entries=self.identity_histogram_entries)
        self.identity_histogram_entries = []


    def set_identity_histogram_meta(self, bins, num_reads_without_identity):
        """Marks identity histograms as stored, and records what is needed to interpret them

        Parameters
        ==========
        bins : bamops.ReadIdentityBins
            The bins the histograms were built with.
        num_reads_without_identity : dict
            The number of reads per sample whose percent identity could not be determined, which are
            not in any histogram (e.g., {'SAMPLE_01': 12}).
        """

        self.db.set_meta_value('identity_histogram_min', bins.min_identity)
        self.db.set_meta_value('identity_histogram_bin_width', bins.bin_width)
        self.db.set_meta_value('identity_histogram_num_bins', bins.num_bins)
        self.db.set_meta_value('identity_definition', 'gap_excluded')
        self.db.set_meta_value('num_reads_without_identity', json.dumps(num_reads_without_identity))
        self.db.set_meta_value('identity_histograms_stored', True)


    def identity_histograms_stored(self):
        """Returns True if this database contains read identity histograms"""

        return bool(self.db.get_meta_value('identity_histograms_stored', return_none_if_not_in_table=True))


    def get_identity_histogram_bins(self):
        """Returns the bamops.ReadIdentityBins the stored histograms were built with, or None if there are none"""

        if not self.identity_histograms_stored():
            return None

        # imported here, because bamops imports dbops, which imports this module
        from anvio.bamops import ReadIdentityBins

        bins = ReadIdentityBins(min_identity=float(self.db.get_meta_value('identity_histogram_min', try_as_type_int=False)),
                                bin_width=float(self.db.get_meta_value('identity_histogram_bin_width', try_as_type_int=False)))

        num_bins = self.db.get_meta_value('identity_histogram_num_bins')
        if bins.num_bins != num_bins:
            raise AuxiliaryDataError(f"The auxiliary database at '{self.db_path}' says its read identity histograms have "
                                     f"{num_bins} bins, but its bin parameters give {bins.num_bins} bins :/")

        return bins


    def get_num_reads_without_identity(self):
        """Returns the number of reads per sample whose percent identity could not be determined (empty if not stored)"""

        if not self.identity_histograms_stored():
            return {}

        return json.loads(self.db.get_meta_value('num_reads_without_identity', try_as_type_int=False))


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
                coverage_array = utils.convert_binary_blob_to_numpy_array(blob, dtype=self.coverage_dtype)

            split_coverage[sample_name] = coverage_array

        return split_coverage


    def get_identity_histograms(self, split_name):
        """Returns the read identity histograms of a split in all samples

        Returns
        =======
        identity_histograms : dict
            {sample_name: (num_below_min, histogram), ...}, where histogram is a 1D numpy array of read
            counts per bin (see self.get_identity_histogram_bins). Empty if this database has no
            identity histograms (see self.identity_histograms_stored).
        """

        if not self.identity_histograms_stored():
            return {}

        cursor = self.db._exec(f"SELECT sample_name, num_below_min, histogram FROM {t.split_identity_histograms_table_name} WHERE split_name = ?", (split_name, ))
        rows = cursor.fetchall()

        if len(rows) == 0:
            raise AuxiliaryDataError(f"The auxiliary database at '{self.db_path}' does not have read identity histograms "
                                     f"for the split '{split_name}'.")

        identity_histograms = {}
        for sample_name, num_below_min, blob in rows:
            if isinstance(blob, int):
                # a histogram of zeros, stored as the number of bins (see self.append_identity_histogram)
                histogram = np.zeros(blob, dtype=IDENTITY_HISTOGRAM_DTYPE)
            else:
                histogram = utils.convert_binary_blob_to_numpy_array(blob, dtype=IDENTITY_HISTOGRAM_DTYPE)

            identity_histograms[sample_name] = (num_below_min, histogram)

        return identity_histograms


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
            self.db._exec(f"""CREATE INDEX IF NOT EXISTS identity_histograms_index ON {t.split_identity_histograms_table_name}(split_name, sample_name)""")

        self.db.disconnect()

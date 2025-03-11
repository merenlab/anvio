# -*- coding: utf-8
# pylint: disable=line-too-long

"""Classes to make sense of single nucleotide variation"""

import copy
import numpy as np

import anvio

from anvio.errors import ConfigError


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "A. Murat Eren"
__email__ = "a.murat.eren@gmail.com"


class VariablityTestFactory:
    def __init__(self, params={'b': 2, 'm': 1.45, 'c': 0.05}):
        self.params = params

        if self.params:
            self.b, self.m, self.c = self.params['b'], self.params['m'], self.params['c']
        else:
            self.b, self.m, self.c = None, None, None


    def get_min_acceptable_departure_from_reference(self, coverage):
        """Get minimum allowable departure from consensus

        Notes
        =====
        - 0 returned if self.params is None
        - https://i.imgur.com/oVoDcbT.png (b=2, m=1.45, c=0.05)
        - https://www.desmos.com/calculator/qwocua4zi5 (interactive plot to tune b, m , and c)
        """

        if self.params is None:
            if hasattr(coverage, '__len__'):
                return np.zeros(len(coverage))
            else:
                return 0

        return (1 / self.b) ** (coverage ** (1 / self.b) - self.m) + self.c


class ProcessAlleleCounts:
    def __init__(self, allele_counts, allele_to_array_index, sequence, sequence_as_index=None, min_coverage_for_variability=1, test_class=None, additional_per_position_data={}):
        """A class to process raw variability information for a given allele counts array

        Creates self.d, a dictionary of equal-length arrays that describes information related to
        variability.

        Parameters
        ==========
        allele_counts : array-like
            An allele counts array. Each column is a position in the sequence, and each row is an
            allele (e.g. A, C, T, G, N if alleles are nucleotides).

        allele_to_array_index : dict
            Which allele belongs at which row index? If A is row 0, C is row 1, etc, the dictionary
            should be {'A': 0, 'C': 1, ...}.

        sequence : str
            What sequence is this for? It should have length equal to number of columns of
            allele_counts

        sequence_as_index : None
            allele_to_array_index provides the means to convert sequence into its index-form.
            However, this requires an expensive list comprehension. If you have already calculated
            the sequence as an index, be sure you provide it here. If you don't provide anything, it
            will be calculated at high cost

        min_coverage_for_variability : int, 1
            positions below this coverage value will be filtered out

        test_class : VariablityTestFactory, None
            If not None, positions will be filtered out if they are deemed not worth reporting

        additional_per_position_data : dict, {}
            This class creates self.d, a dictionary of equal length arrays that describes
            information related to variability. If the user has _other_ data for each position in
            this sequence, they can pass it with parameter. For example, if the user has a
            True/False _array_ (not list) that states whether each position is an outlier position
            relative to a contig, they could pass a dictionary {'cov_outlier_in_contig':
            np.array([True, True, ...])}, where the array is the same length as `sequence`. This
            array will be added to self.d, and will be appropriately filtered alongside the other
            variables

        Notes
        =====
        - Originally self.d was a pandas dataframe. While this approach made the code very
          readable and simple to write, it was extremely slow.
        - If you are analyzing nucleotide, amino acid, or codon variability, you should use the
          inheriting classes ProcessNucleotideCounts, ProcessAminoAcidCounts, and ProcessCodonCounts
        """

        self.d = copy.copy(additional_per_position_data)

        for key in self.d:
            if len(self.d[key]) != allele_counts.shape[1]:
                raise ConfigError("ProcessAlleleCounts :: key '%s' in your passed data dictionary \
                                   has %d positions, but sequence has %d." % (key, len(self.d[key]), len(sequence)))

        if len(sequence) != allele_counts.shape[1]:
            raise ConfigError("ProcessAlleleCounts :: allele_counts has %d positions, but sequence has %d." \
                              % (len(sequence), allele_counts.shape[1]))

        if len(allele_to_array_index) != allele_counts.shape[0]:
            raise ConfigError("ProcessAlleleCounts :: allele_counts has %d rows, but the allele_to_array_index dictionary has %d." \
                              % (allele_counts.shape[0], len(allele_to_array_index)))

        self.min_coverage_for_variability = min_coverage_for_variability
        self.test_class = test_class

        # dictionaries to convert to/from array-row-index and allele
        self.allele_to_array_index = allele_to_array_index
        self.array_index_to_allele = {v: k for k, v in self.allele_to_array_index.items()}

        self.d['pos'] = np.arange(len(sequence))
        self.d['allele_counts'] = allele_counts
        self.d['reference'] = np.array(list(sequence))

        if sequence_as_index is not None:
            self.sequence_as_index_provided = True
            self.d['sequence_as_index'] = sequence_as_index
        else:
            self.sequence_as_index_provided = False

        if self.min_coverage_for_variability < 1:
            raise ConfigError("ProcessAlleleCounts :: self.min_coverage_for_variability must be at least 1, currently %d" % self.min_coverage_for_variability)


    def process(self, skip_competing_items=False):
        """The main function call of this class. Populates self.d"""

        if self.get_data_length() == 0:
            return False

        # remove positions that have non-allowed characters in the sequence
        self.filter_or_dont(self.get_boolean_of_allowable_characters_in_reference(), kind='boolean')
        if self.get_data_length() == 0:
            return False

        if not self.sequence_as_index_provided:
            self.d['sequence_as_index'] = np.array([self.allele_to_array_index[item] for item in self.d['reference']])

        self.d['coverage'] = self.get_coverage()

        # Filter if some positions are not well-covered
        indices_to_keep = self.get_indices_above_coverage_threshold(self.d['coverage'], self.min_coverage_for_variability)
        self.filter_or_dont(indices_to_keep)
        if self.get_data_length() == 0:
            return False

        self.d['reference_coverage'] = self.get_reference_coverage()
        self.d['departure_from_reference'] = self.get_departure_from_reference(self.d['reference_coverage'], self.d['coverage'])

        # Filter if some positions were not worth reporting
        indices_to_keep = self.get_positions_worth_reporting(self.d['coverage'], self.d['departure_from_reference'])
        self.filter_or_dont(indices_to_keep)
        if self.get_data_length() == 0:
            return False

        if not skip_competing_items:
            self.d['competing_items'] = self.get_competing_items(self.d['reference_coverage'], self.d['coverage'])

            # Filter if any competing items are None
            indices_to_keep = self.get_positions_with_competing_items(self.d['competing_items'])
            self.filter_or_dont(indices_to_keep)
            if self.get_data_length() == 0:
                return False

        # each allele gets its own key in self.d
        for index, item in self.array_index_to_allele.items():
            self.d[item] = self.d['allele_counts'][index, :]

        # Delete intermediate keys
        del self.d['allele_counts']
        del self.d['sequence_as_index']
        del self.d['reference_coverage']

        return True


    def get_data_length(self):
        """Get the length of data (number of positions with alleles)"""

        return len(self.d['reference'])


    def filter(self, keys):
        """Filters self.d. keys can be an array-like of indices or array-like of booleans"""

        for key in self.d:
            if self.d[key].ndim == 1:
                self.d[key] = self.d[key][keys]
            else:
                self.d[key] = self.d[key][:, keys]


    def filter_or_dont(self, keys, kind='indices'):
        """Filter self.d if it is required

        Parameters
        ==========
        keys : array-like
            What should be used to filter? If kind == 'indices', it should be an array of indices to
            keep, If kind == 'boolean', it should be an array of booleans

        kind : str, 'indices'
            Either 'indices' or 'boolean'. See keys for what each means
        """

        if kind == 'indices':
            if len(keys) == len(self.d['pos']):
                # Nothing to filter
                return

        elif kind == 'boolean':
            if sum(keys) == len(keys):
                # Nothing to filter
                return

        self.filter(keys)


    def get_coverage(self):
        return np.sum(self.d['allele_counts'], axis=0)


    def get_reference_coverage(self):
        return self.d['allele_counts'][self.d['sequence_as_index'], np.arange(self.d['allele_counts'].shape[1])]


    def get_departure_from_reference(self, reference_coverage=None, coverage=None):
        if reference_coverage is None:
            reference_coverage = self.get_reference_coverage()

        if coverage is None:
            coverage = self.get_coverage()

        return 1 - reference_coverage/coverage


    def get_competing_items(self, reference_coverage=None, coverage=None):
        if reference_coverage is None:
            reference_coverage = self.get_reference_coverage()

        if coverage is None:
            coverage = self.get_coverage()

        n = self.d['allele_counts'].shape[1]

        # as a first pass, sort the row indices (-allele_counts_array is used to sort from highest -> lowest)
        competing_items_as_index = np.argsort(-self.d['allele_counts'], axis=0)

        # take the top 2 items
        competing_items_as_index = competing_items_as_index[:2, :]

        # get the coverage of the second item
        coverage_second_item = self.d['allele_counts'][competing_items_as_index[1, :], np.arange(n)]

        # if the coverage of the second item is 0, set the second index equal to the first
        competing_items_as_index[1, :] = np.where(coverage_second_item == 0, competing_items_as_index[0, :], competing_items_as_index[1, :])

        # sort the competing nts
        competing_items_as_index = np.sort(competing_items_as_index, axis=0)

        # make the competing nts list
        nts_1 = [self.array_index_to_allele[index_1] for index_1 in competing_items_as_index[0, :]]
        nts_2 = [self.array_index_to_allele[index_2] for index_2 in competing_items_as_index[1, :]]
        competing_items = np.fromiter((nt_1 + nt_2 for nt_1, nt_2 in zip(nts_1, nts_2)), np.dtype('<U2'), count=n)

        # If the second item is 0, and the reference is the first item, set competing_items to None.
        # This can easily be checked by seeing if reference_coverage == coverage
        competing_items = np.where(reference_coverage == coverage, None, competing_items)

        return competing_items


    def get_boolean_of_allowable_characters_in_reference(self, sequence=None):
        if sequence is None:
            sequence = self.d['reference']

        items_in_sequence = set(sequence)
        for item in items_in_sequence:
            if item not in self.allele_to_array_index:
                return [item in self.allele_to_array_index for item in sequence]

        return [True] * len(sequence)


    def get_indices_above_coverage_threshold(self, coverage=None, threshold=None):
        if coverage is None:
            coverage = self.get_coverage()

        if threshold is None:
            threshold = self.min_coverage_for_variability

        return np.where(coverage >= threshold)[0]


    def get_positions_worth_reporting(self, coverage, departure_from_reference):
        worth_reporting = np.array([True] * len(coverage))

        if not self.test_class:
            return worth_reporting

        threshold = self.test_class.get_min_acceptable_departure_from_reference(coverage)

        return np.where(departure_from_reference >= threshold)[0]


    def get_positions_with_competing_items(self, competing_items):

        return np.where(competing_items != None)[0]


    def rename_key(self, from_this, to_that):
        if from_this in self.d:
            self.d[to_that] = self.d.pop(from_this)


class ProcessNucleotideCounts(ProcessAlleleCounts):
    def __init__(self, *args, **kwargs):
        ProcessAlleleCounts.__init__(self, *args, **kwargs)

    def process(self, *args, **kwargs):
        p = ProcessAlleleCounts.process(self, *args, **kwargs)
        self.rename_key('competing_items', 'competing_nts')
        return p


class ProcessAminoAcidCounts(ProcessAlleleCounts):
    def __init__(self, *args, **kwargs):
        ProcessAlleleCounts.__init__(self, *args, **kwargs)

    def process(self, *args, **kwargs):
        p = ProcessAlleleCounts.process(self, *args, **kwargs)
        self.rename_key('competing_items', 'competing_aas')
        return p


class ProcessCodonCounts(ProcessAlleleCounts):
    def __init__(self, *args, **kwargs):
        ProcessAlleleCounts.__init__(self, *args, **kwargs)

    def process(self, *args, **kwargs):
        p = ProcessAlleleCounts.process(self, *args, **kwargs)
        self.rename_key('competing_items', 'competing_codons')
        self.rename_key('pos', 'codon_order_in_gene')
        return p


class ProcessIndelCounts(object):
    def __init__(self, indels, coverage, min_indel_fraction=0, test_class=None, min_coverage_for_variability=1):
        """A class to process raw variability information for a given allele counts array

        Creates self.d, a dictionary of equal-length arrays that describes information related to
        variability.

        Parameters
        ==========
        indels : dictionary
            A dictionary that looks like this:

                {
                    6279666787066445523: OrderedDict([
                        ('split_name', 'IGD_000000000532_split_00001'),
                        ('pos', 2),
                        ('pos_in_contig', 2),
                        ('corresponding_gene_call', 25396),
                        ('in_noncoding_gene_call', 0),
                        ('in_coding_gene_call', 1),
                        ('base_pos_in_codon', 3),
                        ('codon_order_in_gene', 0),
                        ('cov_outlier_in_split', 1),
                        ('cov_outlier_in_contig', 1),
                        ('reference', 'T'),
                        ('type', 'INS'),
                        ('sequence', 'CTGACGGCT'),
                        ('length', 9),
                        ('count', 1)
                    ]),
                    -5035942137885303221: OrderedDict([
                        ('split_name', 'IGD_000000000532_split_00001'),
                        ('pos', 0),
                        ('pos_in_contig', 0),
                        ('corresponding_gene_call', 25396),
                        ('in_noncoding_gene_call', 0),
                        ('in_coding_gene_call', 1),
                        ('base_pos_in_codon', 1),
                        ('codon_order_in_gene', 0),
                        ('cov_outlier_in_split', 1),
                        ('cov_outlier_in_contig', 1),
                        ('reference', 'G'),
                        ('type', 'INS'),
                        ('sequence', 'CTCACGG'),
                        ('length', 7),
                        ('count', 1)
                    ]),
                    ...
                }

            The keys are unique identifiers. The OrderedDicts should have at least the key `pos`,
            but there are no restrictions on what other keys it may have.

        coverage : array
            What is the coverage for the sequence this is for? This should have length equal to sequence

        test_class : VariablityTestFactory, None
            If not None, indels will be filtered out if they are deemed not worth reporting

        min_coverage_for_variability : int, 1
            positions below this coverage value will be filtered out
        """

        self.indels = indels
        self.coverage = coverage
        self.test_class = test_class if test_class is not None else VariablityTestFactory(params=None)
        self.min_coverage_for_variability = min_coverage_for_variability


    def process(self):
        """Modify self.indels"""

        indel_hashes_to_remove = set()
        for indel_hash in self.indels:

            indel = self.indels[indel_hash]
            pos = indel['pos']

            # Calculate coverage
            if indel['type'] == 'INS':
                if pos == len(self.coverage)-1:
                    # This is the last position in the sequence. so coverage based off only the
                    # NT left of the indel
                    cov = self.coverage[pos]
                else:
                    # The coverage is the average of the coverage left and right of the
                    # insertion
                    cov = (self.coverage[pos] + self.coverage[pos+1])/2
            else:
                # The coverage is the average of the NT coverages that the deletion occurs over
                # plus the number of reads that have the deletion at that position
                cov = np.mean(self.coverage[pos:pos+indel['length']]) + indel['count']

            # Filter the entry if need be
            if cov < self.min_coverage_for_variability:
                # coverage of corresponding position is not high enough
                indel_hashes_to_remove.add(indel_hash)
                continue

            # NOTE We call get_min_acceptable_departure_from_reference, yet the threshold value it
            # spits out is being compared to count/coverage, since there is no "departure from
            # reference" for indels.
            if indel['count']/cov < self.test_class.get_min_acceptable_departure_from_reference(cov):
                # indel count is not high enough compared to coverage value
                indel_hashes_to_remove.add(indel_hash)
                continue

            self.indels[indel_hash]['coverage'] = cov

        self.indels = {k: v for k, v in self.indels.items() if k not in indel_hashes_to_remove}



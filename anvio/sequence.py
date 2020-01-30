# -*- coding: utf-8
# pylint: disable=line-too-long

'''Primitive classes for basic DNA sequence properties.'''

import numpy
import collections

from itertools import permutations

import anvio
import anvio.constants as constants

from anvio.errors import ConfigError

__author__ = "Developers of anvi'o (see AUTHORS.txt)"
__copyright__ = "Copyleft 2015-2018, the Meren Lab (http://merenlab.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "A. Murat Eren"
__email__ = "a.murat.eren@gmail.com"
__status__ = "Development"


class Codon:
    def __init__(self):
        pass

    def get_codon_to_codon_sequence_trajectory(self, start_codon, end_codon, as_amino_acids=False):
        '''
        this returns a list of all possible sequence trajectories to get from one codon to another
        assuming the least amount of mutations necessary. if as_amino_acids, the trajectories will
        be converted into amino acid space
        '''
        indices_of_variation = []
        for i in range(3):
            if start_codon[i] != end_codon[i]:
                indices_of_variation.append(i)
        index_trajectories = list(permutations(indices_of_variation))

        all_trajectories = []
        for index_trajectory in index_trajectories:
            sequence_trajectory = [start_codon]
            mutate = list(start_codon)
            for index in index_trajectory:
                mutate[index] = end_codon[index]
                sequence_trajectory.append(''.join(mutate))
            all_trajectories.append(sequence_trajectory)

        if as_amino_acids:
            # each codon is converted to an amino acid. if two adjacent codons are the same amino
            # acid then only one is kept
            for i, trajectory in enumerate(all_trajectories):
                for j, node in enumerate(trajectory):
                    trajectory[j] = constants.codon_to_AA[node]

                # gets rid of duplicates
                all_trajectories[i] = list(dict.fromkeys(trajectory))

        return all_trajectories


    def get_codon_to_codon_dist_dictionary(self):
        """
        Returns a dictionary containing the number the nucleotide difference
        between two codons, and the number of transitions & transversions required to
        mutate from one codon to another. 

        USAGE:
        =====
        >>> dist['TTC']['AAA']
        (3, 2, 1)

        Here, 3 is the total number of nucleotide differences, 2 is the number
        of transitions, and 1 is the number of transversions. Of course, the number of
        transitions added to the number of transversions is equal to the number of
        nucleotide differences.
        """

        codons = list(constants.codon_to_AA.keys())
        dist = {}

        mutation_type = {"AT": "transition",
                         "CG": "transition",
                         "AG": "transversion",
                         "AC": "transversion",
                         "CT": "transversion",
                         "GT": "transversion"}

        for start_codon in codons:
            dist[start_codon] = {}

            for end_codon in codons:

                # s = number transitions 
                # v = number transversions
                s = 0; v = 0
                for nt_pos in range(3):

                    pair = ''.join(sorted(start_codon[nt_pos] + end_codon[nt_pos]))

                    if pair[0] == pair[1]:
                        continue
                    if mutation_type[pair] == "transition":
                        s += 1
                    if mutation_type[pair] == "transversion":
                        v += 1

                dist[start_codon][end_codon] = (s+v, s, v)

        return dist


class Read:
    def __init__(self, read, _run_test_class=False):
        """Class for manipulating reads

        Some of these class methods parse and manipulate cigar strings. You can read up on cigar
        operations here:
        https://pysam.readthedocs.io/en/latest/api.html#pysam.AlignedSegment.cigartuples

        Parameters
        ==========
        read : pysam.AlignedSegment

        _run_test_class : bool, False
            For developers. During instantiation of class instance ReadTestClass will be run
        """
        if _run_test_class:
            test = ReadTestClass()
            test.test_trim()

        self.r = read

        # redefine all properties of interest explicitly from pysam.AlignedSegment object as
        # attributes of this class. The reason for this is that some of the AlignedSegment
        # attributes have no __set__ methods, so are read only. Since this class is designed to
        # modify some of these attributes, and since we want to maintain consistency across
        # attributes, all attributes of interest are redefined here
        self.query_alignment_sequence = self.r.query_alignment_sequence
        self.cigartuples = self.r.cigartuples
        self.reference_positions = self.r.get_reference_positions()
        self.reference_start = self.r.reference_start
        self.reference_end = self.r.reference_end


    def get_blocks(self):
        """Mimic the get_blocks function from AlignedSegment.

        Calculates directly from self.reference_positions

        Modified from:
        https://stackoverflow.com/questions/7352684/how-to-find-the-groups-of-consecutive-elements-from-an-array-in-numpy/7353335#7353335

        Examples
        ========
        def consecutive(data, stepsize=1):
            return [(x[0], x[-1]+1) for x in numpy.split(data, numpy.where(numpy.diff(data) != stepsize)[0]+1)]

        a = numpy.array([0, 47, 48, 49, 50, 97, 98, 99])
        consecutive(a)

        >>> [(0, 1), (47, 51), (97, 100)]
        """

        return [(x[0], x[-1] + 1) for x in numpy.split(self.reference_positions, numpy.where(numpy.diff(self.reference_positions) != 1)[0] + 1)]


    def get_aligned_sequence(self):
        """Get the aligned sequence at each position in self.reference_positions

        Notes
        =====
        - This method exists because self.read.query_alignment_sequence is the read sequence with
          soft clipping removed, but it is otherwise 'unaligned' to the reference. For example,
          self.read.query_alignment_sequence does not even necessarily have the same length as
          self.read.get_reference_positions() due to indels. To get the aligned sequence, we have to
          parse the cigar string to build `aligned_sequence`, which gives us the base
          contributed by this read at each of its aligned positions.
        """
        sequence = self.query_alignment_sequence
        cigar_tuples = self.cigartuples
        aligned_sequence = ''

        read_pos = 0
        for operation, length in cigar_tuples:
            if operation == 0:
                # there is a mapping segment
                aligned_sequence += sequence[read_pos:(read_pos + length)]
                read_pos += length
            elif operation == 1:
                # there is an insertion in the read
                read_pos += length
            elif operation == 2:
                # there is a gap in the read
                pass
            else:
                # FIXME
                pass

        return aligned_sequence


    def trim(self, trim_by, side='left'):
        """Trims self.read by either the left or right

        Modifies the attributes:

            query_alignment_sequence
            cigartuples
            reference_positions

        Do not expect more than this!

        Parameters
        ==========
        trim_by : int
            The number of REFERENCE bases you would like to trim the read by

        side : str, 'left'
            Either 'left' or 'right' side.
        """
        cigar_tuples = self.cigartuples
        read_sequence = self.query_alignment_sequence
        reference_positions = self.reference_positions

        tuple_indices_to_remove = []
        trimmed_tuple = None
        count = trim_by
        m, n = 0, 0

        if side == 'right':
            # flip the read
            read_sequence = read_sequence[::-1]
            cigar_tuples = cigar_tuples[::-1]
            reference_positions = reference_positions[::-1]

        for i, cigar_tuple in enumerate(cigar_tuples):
            operation, length = cigar_tuple
            tuple_indices_to_remove.append(i)

            if operation == 0:
                if length > count:
                    trimmed_tuple = (operation, length - count)
                    count = 0
                    break
                else:
                    count -= length

            elif operation == 1:
                m += length

            elif operation == 2:
                if length > count:
                    trimmed_tuple = (operation, length - count)
                    n += length - count
                    count = 0
                    break
                else:
                    count -= length
                    n += length

            if count == 0:
                break

        cigar_tuples = [cigar_tuple for i, cigar_tuple in enumerate(cigar_tuples) if i not in tuple_indices_to_remove]
        if trimmed_tuple:
            cigar_tuples.insert(0, trimmed_tuple)

        read_sequence = read_sequence[trim_by + m - n:]
        reference_positions = reference_positions[trim_by - n:]

        if side == 'right':
            # flip the read back
            read_sequence = read_sequence[::-1]
            cigar_tuples = cigar_tuples[::-1]
            reference_positions = reference_positions[::-1]

        # overwrite the attributes of self
        self.cigartuples = cigar_tuples
        self.query_alignment_sequence = read_sequence
        self.reference_positions = reference_positions
        self.reference_start = reference_positions[0]
        self.reference_end = reference_positions[-1]


class ReadTestClass:
    """Small test class for Read"""

    def make_read(self, cigartuples, reference_positions):
        class Read:
            def get_reference_positions(self):
                return reference_positions

        read = Read()
        read.query_alignment_sequence = 'AACCTTGG'
        read.cigartuples = cigartuples
        read.reference_start = 0 # unused
        read.reference_end = 0 # unused

        return read


    def test_trim(self):
        """Tests Read.trim for a number of cases

        Asserts that the trimmed cigartuples, reference_positions, and sequences are trimmed correctly
        for the following cases. Expected answers were created manually based on these diagrams:

        CASE #1
        =======
        A A C C T T G G
        A C T G A C T G A C T G = reference
        [(0,8)]

        CASE #2
        =======
        A A C C T T G G
        A C - - - - T G A C T G A C T G = reference
        [(0,2), (1,4), (0,2)]

        CASE #3
        =======
        A A C C - - T T G G
        A - - - C T G A C T G A C T G = reference
        [(0,1), (1,3), (2,2), (0,4)]

        CASE #4
        =======
        A - A - C - C - T T G G
        A C T G A C T G A C - T G = reference
        [(0,1), (2,1), (0,1), (2,1), (0,1), (2,1), (0,1), (2,1), (0,2), (1,1), (0,1)]
        """

        test_sets = [
            {
                'input_cigartuples': [(0,8)],
                'input_reference_positions': [0,1,2,3,4,5,6,7],
                'output_tuples_left': [(0,5)],
                'output_sequence_left': 'CTTGG',
                'output_reference_positions_left': [3,4,5,6,7],
                'output_tuples_right': [(0,5)],
                'output_sequence_right': 'AACCT',
                'output_reference_positions_right': [0,1,2,3,4],
            },
            {
                'input_cigartuples': [(0,2), (1,4), (0,2)],
                'input_reference_positions': [0,1,2,3],
                'output_tuples_left': [(0,1)],
                'output_sequence_left': 'G',
                'output_reference_positions_left': [3],
                'output_tuples_right': [(0,1)],
                'output_sequence_right': 'A',
                'output_reference_positions_right': [0],
            },
            {
                'input_cigartuples': [(0,1), (1,3), (2,2), (0,4)],
                'input_reference_positions': [0,3,4,5,6],
                'output_tuples_left': [(0,4)],
                'output_sequence_left': 'TTGG',
                'output_reference_positions_left': [3,4,5,6],
                'output_tuples_right': [(0,1), (1,3), (2,2), (0,1)],
                'output_sequence_right': 'AACCT',
                'output_reference_positions_right': [0,3],
            },
            {
                'input_cigartuples': [(0,1), (2,1), (0,1), (2,1), (0,1), (2,1), (0,1), (2,1), (0,2), (1,1), (0,1)],
                'input_reference_positions': [0,2,4,6,8,9,10],
                'output_tuples_left': [(2,1), (0,1), (2,1), (0,1), (2,1), (0,2), (1,1), (0,1)],
                'output_sequence_left': 'CCTTGG',
                'output_reference_positions_left': [4,6,8,9,10],
                'output_tuples_right': [(0,1), (2,1), (0,1), (2,1), (0,1), (2,1), (0,1), (2,1)],
                'output_sequence_right': 'AACC',
                'output_reference_positions_right': [0,2,4,6],
            },
        ]

        for test_set in test_sets:
            read = Read(self.make_read(
                test_set['input_cigartuples'],
                test_set['input_reference_positions']
            ))
            read.trim(trim_by=3, side='left')
            assert read.query_alignment_sequence == test_set['output_sequence_left']
            assert read.cigartuples == test_set['output_tuples_left']
            assert read.reference_positions == test_set['output_reference_positions_left']

            read = Read(self.make_read(
                test_set['input_cigartuples'],
                test_set['input_reference_positions']
            ))
            read.trim(trim_by=3, side='right')
            assert read.query_alignment_sequence == test_set['output_sequence_right']
            assert read.cigartuples == test_set['output_tuples_right']
            assert read.reference_positions == test_set['output_reference_positions_right']


class Composition:
    def __init__(self, sequence):
        self.sequence = sequence
        self.GC_content = 0.0

        self.report()


    def report(self):
        s = self.sequence
        raw_length = len(s)

        self.A = s.count('A')
        self.T = s.count('T')
        self.C = s.count('C')
        self.G = s.count('G')
        self.N = raw_length - (self.A + self.T + self.C + self.G)
        length = raw_length - self.N

        if not length:
            # sequence is composed of only N's
            self.GC_content = 0.0
        else:
            self.GC_content = (self.G + self.C) * 1.0 / length


class Coverage:
    def __init__(self):
        self.c = None # becomes a numpy array of coverage values
        self.outlier_positions = set([]) # set of positions along the sequence, coverage values of which
                                         # are classified as outliers; see `get_indices_for_outlier_values`
        self.min = 0
        self.max = 0
        self.std = 0.0
        self.mean = 0.0
        self.median = 0.0
        self.detection = 0.0
        self.mean_Q2Q3 = 0.0

        self.routine_dict = {
            'approximate': self._approximate_routine,
            'accurate': self._accurate_routine,
        }


    def run(self, bam, contig_or_split, start=None, end=None, method='accurate', **kwargs):
        """Loop through the bam pileup and calculate coverage over a defined region of a contig or split

        Parameters
        ==========
        bam : bamops.BAMFileObject
            Init such an object the way you would a pysam.AlignmentFile, i.e. bam =
            bamops.BAMFileObject(path_to_bam)

        contig_or_split : anvio.contigops.Split or anvio.contigops.Contig or str
            If Split object is passed, and `start` or `end` are None, they are automatically set to
            contig_or_split.start and contig_or_split.end. If str object is passed, it is assumed to
            be a contig name

        start : int
            The index start of where coverage is calculated. Relative to the contig, even when
            `contig_or_split` is a Split object. 

        end : int
            The index end of where coverage is calculated. Relative to the contig, even when
            `contig_or_split` is a Split object.

        method : string
            How do you want to calculate? Options: ('accurate', 'approximate'). 'accurate' accounts
            for gaps in the alignment, 'approximate' does not. For others, see associated methods
            and pass special parameters they take through **kwargs
        """

        # if there are defined start and ends we have to trim reads so their ranges fit inside self.c
        iterator = bam.fetch if (start is None and end is None) else bam.fetch_and_trim

        if isinstance(contig_or_split, anvio.contigops.Split):
            contig_name = contig_or_split.parent
            start = contig_or_split.start if not start else start
            end = contig_or_split.end if not end else end

        elif isinstance(contig_or_split, anvio.contigops.Contig):
            contig_name = contig_or_split.name
            start = 0 if not start else start
            end = contig_or_split.length if not end else end

        elif isinstance(contig_or_split, str):
            contig_name = contig_or_split
            start = 0 if not start else start
            end = bam.get_reference_length(contig_name) if not end else end

        else:
            raise ConfigError("Coverage.run :: You can't pass an object of type %s as contig_or_split" % type(contig_or_split))

        # a coverage array the size of the defined range is allocated in memory
        c = numpy.zeros(end - start).astype(int)

        routine = self.routine_dict.get(method)
        if not routine:
            raise ConfigError("Coverage :: %s is not a valid method.")

        self.c = routine(c, bam, contig_name, start, end, iterator, **kwargs)

        if len(self.c):
            try:
                contig_or_split.explicit_length = len(self.c)
            except AttributeError:
                pass
            self.process_c(self.c)


    def _approximate_routine(self, c, bam, contig_name, start, end, iterator):
        """Routine that does not account for gaps in alignment

        Notes
        =====
        - Should typically not be called explicitly. Use run instead
        """

        for read in iterator(contig_name, start, end):
            c[read.reference_start:read.reference_end] += 1

        return c


    def _accurate_routine(self, c, bam, contig_name, start, end, iterator):
        """Routine that accounts for gaps in the alignment

        Notes
        =====
        - Should typically not be called explicitly. Use run instead
        - This strategy was also considered, but is much slower because it uses fancy-indexing
          https://jakevdp.github.io/PythonDataScienceHandbook/02.07-fancy-indexing.html:

          for read in bam.fetch(contig_name, start, end):
              r = read.get_reference_positions()
              c[r] += 1
        """

        for read in iterator(contig_name, start, end):
            for block in read.get_blocks():
                c[block[0]:block[1]] += 1

        return c


    def process_c(self, c):
        self.min = numpy.amin(c)
        self.max = numpy.amax(c)
        self.median = numpy.median(c)
        self.mean = numpy.mean(c)
        self.std = numpy.std(c)
        self.detection = 1 - (float(collections.Counter(c)[0]) / len(c))

        self.outlier_positions = get_indices_for_outlier_values(c)

        if c.size < 4:
            self.mean_Q2Q3 = self.mean
        else:
            sorted_c = sorted(c)
            Q = int(c.size * 0.25)
            Q2Q3 = sorted_c[Q:-Q]
            self.mean_Q2Q3 = numpy.mean(Q2Q3)


def get_indices_for_outlier_values(c):
    is_outlier = get_list_of_outliers(c)
    return set([p for p in range(0, c.size) if is_outlier[p]])


def get_list_of_outliers(values, threshold=None, zeros_are_outliers=False):
    """
    Returns a boolean array with True if values are outliers and False
    otherwise.

    Modified from Joe Kington's (https://stackoverflow.com/users/325565/joe-kington)
    implementation computing absolute deviation around the median.

    Parameters:
    -----------
        values    : An numobservations by numdimensions array of observations
        threshold : The modified z-score to use as a thresholdold. Observations with
                    a modified z-score (based on the median absolute deviation) greater
                    than this value will be classified as outliers.

    Returns:
    --------
        mask : A numobservations-length boolean array.

    References:
    ----------
        Boris Iglewicz and David Hoaglin (1993), "Volume 16: How to Detect and
        Handle Outliers", The ASQC Basic References in Quality Control:
        Statistical Techniques, Edward F. Mykytka, Ph.D., Editor.

        http://www.sciencedirect.com/science/article/pii/S0022103113000668
    """

    if threshold is None:
        threshold = 1.5

    if len(values.shape) == 1:
        values = values[:, None]

    median = numpy.median(values, axis=0)

    diff = numpy.sum((values - median) ** 2, axis=-1)
    diff = numpy.sqrt(diff)
    median_absolute_deviation = numpy.median(diff)

    if not median_absolute_deviation:
       if values[0] == 0:
            # A vector of all zeros is considered "all outliers"
            return numpy.array([True] * values.size)
       else:
            # A vector of uniform non-zero values is "all non-outliers"
            # This could be important for silly cases (like in megahit) in which there is a maximum value for coverage
            return numpy.array([False] * values.size)

    modified_z_score = 0.6745 * diff / median_absolute_deviation
    non_outliers = modified_z_score > threshold

    if not zeros_are_outliers:
        return non_outliers
    else:
        zero_positions = [x for x in range(len(values)) if values[x] == 0]
        for i in zero_positions:
            non_outliers[i] = True
        return non_outliers

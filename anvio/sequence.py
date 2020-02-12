# -*- coding: utf-8
# pylint: disable=line-too-long

'''Primitive classes for basic DNA sequence properties.'''

import copy
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


class Cigar:
    def __init__(self):
        """Gives some useful information about cigar string and their operations

        Here are the possible bam operations.

            M       BAM_CMATCH      0
            I       BAM_CINS        1
            D       BAM_CDEL        2
            N       BAM_CREF_SKIP   3
            S       BAM_CSOFT_CLIP  4
            H       BAM_CHARD_CLIP  5
            P       BAM_CPAD        6
            =       BAM_CEQUAL      7
            X       BAM_CDIFF       8

        Notes
        =====
        - A description of what possible cigar operations are possible, see
          https://imgur.com/a/fiQZXNg, which comes from here:
          https://samtools.github.io/hts-specs/SAMv1.pdf
        """

        # key = cigar operation, val[0]/val[1] is whether operation consumes the read/reference
        self.consumes = {
            0: (True, True),
            1: (True, False),
            2: (False, True),
            3: (False, True),
            4: (True, False),
            5: (False, False),
            6: (False, False),
            7: (True, True),
            8: (True, True),
        }


    def iterate(self, cigar_tuples):
        """Iterate through a list of cigar tuples

        Parameters
        ==========
        cigar_tuples : list
            Each element is a (operation, length) tuple. Operations are integers, see class
            docstring

        Yields
        ======
        output : tuple
            (operation, length, consumes_read, consumes_ref) -> (int, int, bool, bool)
        """

        for operation, length in cigar_tuples:
            yield (operation, length, *self.consumes[operation])


class Read:
    def __init__(self, read):
        """Class for manipulating reads

        Parameters
        ==========
        read : pysam.AlignedSegment
        """

        self.r = read
        self.cigarops = Cigar()

        # redefine all properties of interest explicitly from pysam.AlignedSegment object as
        # attributes of this class. The reason for this is that some of the AlignedSegment
        # attributes have no __set__ methods, so are read only. Since this class is designed to
        # modify some of these attributes, and since we want to maintain consistency across
        # attributes, all attributes of interest are redefined here
        self.query_sequence = self.r.query_sequence
        self.cigartuples = self.r.cigartuples
        self.reference_positions = self.r.get_reference_positions()
        self.reference_start = self.r.reference_start
        self.reference_end = self.r.reference_end


    def get_blocks(self):
        """Mimic the get_blocks function from AlignedSegment.

        Notes
        =====
        - Takes roughly 200us
        """

        blocks = []
        block_start = self.reference_positions[0]
        block_length = 0

        for _, length, consumes_read, consumes_ref in self.cigarops.iterate(self.cigartuples):

            if consumes_read and consumes_ref:
                block_length += length

            elif consumes_read and not consumes_ref:
                if block_length:
                    blocks.append((block_start, block_start + block_length))

                block_start = block_start + block_length
                block_length = 0

            elif not consumes_read and consumes_ref:
                if block_length:
                    blocks.append((block_start, block_start + block_length))

                block_start = block_start + block_length + length
                block_length = 0

            else:
                pass

        if block_length:
            blocks.append((block_start, block_start + block_length))

        return blocks


    def get_reference_positions(self):
        """Mimic the get_reference_positions function from AlignedSegment."""

        return self.reference_positions


    def get_aligned_sequence(self):
        """Get the aligned sequence at each position in self.reference_positions

        Notes
        =====
        - This method exists because self.r.query_alignment_sequence doesn't return the read's
          nucleotides at the positions it aligns, it only returns the read after removing
          softclipping. It is otherwise unaligned.  To get the aligned sequence, we have to parse
          the cigar string to build `aligned_sequence`, which gives us the base contributed by this
          read at each of its aligned positions.
        - Takes anywhere from 150-450us
        """

        aligned_sequence = ''

        read_pos = 0
        for _, length, consumes_read, consumes_ref in self.cigarops.iterate(self.cigartuples):

            if consumes_read:
                if consumes_ref:
                    aligned_sequence += self.query_sequence[read_pos:(read_pos + length)]

                read_pos += length

        return aligned_sequence


    def __getitem__(self, key):
        """Splice the read based on reference positions

        Makes a new copy, current instance remains unmodified.

        Parameters
        ==========
        key : slice
            See Examples

        Returns
        =======
        output : Read
            A new Read object

        Examples
        ========

        You have a read with 100 length and self.reference_start = 4500 and self.reference_end =
        4601. You want to splice the segment [4550, 4575).

        >>> type(read)
        <class 'anvio.sequence.Read'>
        >>> spliced_segment = read[4550:4575]
        >>> type(spliced_segment)
        <class 'anvio.sequence.Read'>
        >>> spliced_segment.reference_start
        4550
        >>> spliced_segment.reference_end
        4575

        Specify only one bound:

        >>> print(read[:4575].reference_start, read[:4575].reference_end)
        (4500, 4575)
        >>> print(read[4575:].reference_start, read[4575:].reference_end)
        (4575, 4600)

        Specifying outside read range:

        >>> print(read[4400:4700].reference_start, read[4400:4700].reference_end)
        (4500, 4601)
        """

        if not isinstance(key, slice) or key.step is not None:
            raise ValueError("Read class only supports basic slicing for indexing, e.g. read[start:stop], read[:stop]")

        segment = copy.deepcopy(self)

        start = key.start if key.start is not None else segment.reference_start
        end = key.stop if key.stop is not None else segment.reference_end

        segment.trim(start - segment.reference_start, side='left')
        segment.trim(segment.reference_end - end, side='right')

        return segment


    def trim(self, trim_by, side='left'):
        """Trims self.read by either the left or right

        Modifies the attributes:

            query_sequence
            cigartuples
            reference_positions
            reference_start
            reference_end

        Do not expect more than this!

        Parameters
        ==========
        trim_by : int
            The number of REFERENCE bases you would like to trim the read by. If the trim leaves
            operations that are consumed by the reference but not the read, or the read but not the
            reference, these are trimmed AS WELL. For example, if after trimming by `trim_by`, the
            final cigar string is [(2,2),(0,4)], this will be further trimmed to [(0,4)], since
            there is no useful information held in a terminal read gap.

        side : str, 'left'
            Either 'left' or 'right' side.

        Notes
        =====
        - Takes roughly 250us
        """

        if trim_by == 0:
            return

        elif trim_by < 0:
            raise ConfigError("Read.trim :: Requesting to trim an amount %d, which is negative." % trim_by)

        elif trim_by > self.reference_end - self.reference_start:
            raise ConfigError("Read.trim :: Requesting to trim an amount %d that exceeds the alignment"
                              " range of %d" % (trim_by, self.reference_end - self.reference_start))

        cigar_tuples = self.cigartuples

        if side == 'right':
            # flip the cigar tuple
            cigar_tuples = cigar_tuples[::-1]

        trimmed_cigar_tuples = cigar_tuples.copy()

        ref_positions_trimmed = 0
        read_positions_trimmed = 0
        terminate, terminate_next = (False, False)
        for operation, length, consumes_read, consumes_ref in self.cigarops.iterate(cigar_tuples):

            if consumes_ref and consumes_read:
                if terminate_next:
                    break

                remaining = trim_by - ref_positions_trimmed

                if length > remaining:
                    # the length of the operation exceeds the required trim amount. So we will
                    # terminate this iteration. To trim the cigar tuple, we replace it with a
                    # truncated length
                    trimmed_cigar_tuples[0] = (operation, length - remaining)
                    ref_positions_trimmed += remaining
                    read_positions_trimmed += remaining
                    break

                elif length == remaining:
                    # This is a rare case: The are no more reference positions that need to be
                    # trimmed, but it could be the case that the next few cigartuples consume the
                    # read and not the reference, and these still need to be trimmed. Hence,
                    # `terminate_next` is set to True so the next time a reference-consuming
                    # operation occurs, the loop is immediately terminated
                    terminate_next = True

                ref_positions_trimmed += length
                read_positions_trimmed += length

            else:
                ref_positions_trimmed += length if consumes_ref else 0
                read_positions_trimmed += length if consumes_read else 0

            trimmed_cigar_tuples.pop(0)

        # set new cigartuples
        self.cigartuples = trimmed_cigar_tuples if side == 'left' else trimmed_cigar_tuples[::-1]

        # set new query_sequence
        self.query_sequence = self.query_sequence[read_positions_trimmed:] if side == 'left' else self.query_sequence[:-read_positions_trimmed]

        # set new reference_positions
        if side == 'right':
            num_pos = len(self.reference_positions)
            last_pos = self.reference_positions[-1]
            cutoff_index = next(num_pos - i for i, pos in enumerate(self.reference_positions[::-1]) if last_pos - pos >= trim_by)
            self.reference_positions = self.reference_positions[:cutoff_index]
        else:
            cutoff_index = next(i for i, pos in enumerate(self.reference_positions) if pos - self.reference_start >= trim_by)
            self.reference_positions = self.reference_positions[cutoff_index:]

        self.reference_start = self.reference_positions[0]
        self.reference_end = self.reference_positions[-1] + 1


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


    def run(self, bam, contig_or_split, start=None, end=None, method='accurate', max_coverage=None, **kwargs):
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

        if max_coverage is not None:
            if numpy.max(self.c) > max_coverage:
                self.c[self.c > max_coverage] = max_coverage

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
        - The following strategy was also considered, but is much slower because it uses fancy-indexing
          https://jakevdp.github.io/PythonDataScienceHandbook/02.07-fancy-indexing.html:

          for read in bam.fetch(contig_name, start, end):
              r = read.get_reference_positions()
              c[r] += 1
        """

        for read in iterator(contig_name, start, end):
            for start, end in read.get_blocks():
                c[start:end] += 1

        return c


    def process_c(self, c):
        self.min = numpy.amin(c)
        self.max = numpy.amax(c)
        self.median = numpy.median(c)
        self.mean = numpy.mean(c)
        self.std = numpy.std(c)
        self.detection = 1 - (float(collections.Counter(c)[0]) / len(c))

        self.is_outlier = get_list_of_outliers(c) # this is an array not a list

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

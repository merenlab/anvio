# -*- coding: utf-8
# pylint: disable=line-too-long

'''Primitive classes for basic DNA sequence properties.'''

import numpy
import collections

import anvio

__author__ = "A. Murat Eren"
__copyright__ = "Copyright 2015, The anvio Project"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "A. Murat Eren"
__email__ = "a.murat.eren@gmail.com"
__status__ = "Development"


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
        self.c = [] # list of coverage values
        self.outlier_positions = set([]) # set of positions along the sequence, coverage values of which
                                         # are classified as outliers; see `get_indices_for_outlier_values`
        self.min = 0
        self.max = 0
        self.std = 0.0
        self.mean = 0.0
        self.median = 0.0
        self.detection = 0.0
        self.mean_Q2Q3 = 0.0


    def run(self, bam, split):
        coverage_profile = {}
        for pileupcolumn in bam.pileup(split.parent, split.start, split.end):
            if pileupcolumn.pos < split.start or pileupcolumn.pos >= split.end:
                continue

            coverage_profile[pileupcolumn.pos] = pileupcolumn.n

        for i in range(split.start, split.end):
            if i in coverage_profile:
                self.c.append(coverage_profile[i])
            else:
                self.c.append(0)

        if self.c:
            split.explicit_length = len(self.c)
            self.process_c(self.c)

    def process_c(self, c):
        c = numpy.asarray(c)
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


def get_list_of_outliers(values, threshold=1.5):
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

    if len(values.shape) == 1:
        values = values[:, None]

    median = numpy.median(values, axis=0)

    diff = numpy.sum((values - median) ** 2, axis=-1)
    diff = numpy.sqrt(diff)
    median_absolute_deviation = numpy.median(diff)

    if not median_absolute_deviation:
        return [True] * values.size

    modified_z_score = 0.6745 * diff / median_absolute_deviation

    return modified_z_score > threshold

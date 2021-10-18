#! /usr/bin/env python
# -*- coding: utf-8
"""Unit tests for the Read class.

FIXME This class should deal with all the features of the Read class, such as vectorize,
get_aligned_sequence_and_positions, iterate_cigartuples, etc.
"""

import copy
import numpy as np
import unittest

import anvio

from anvio.bamops import Read

__author__ = "Developers of anvi'o (see AUTHORS.txt)"
__copyright__ = "Copyleft 2015-2019, the Meren Lab (http://merenlab.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "Evan Kiefl"
__email__ = "kiefl.evan@gmail.com"


class FakePySamAlignedSegment:
    def __init__(self, **output):
        for key, v in output.items():
            setattr(self, key, v)


class TestTrim(unittest.TestCase):
    """
    CASE #1
    =======
    <anvio.bamops.Read object at 0x12b042470>
     ├── start, end : [0, 8)
     ├── cigartuple : [(0, 8)]
     ├── read       : AACCTTGG
     └── reference  : XXXXXXXX

    CASE #2
    =======
    <anvio.bamops.Read object at 0x12b042358>
     ├── start, end : [0, 4)
     ├── cigartuple : [(0, 2), (1, 4), (0, 2)]
     ├── read       : AACCTTGG
     └── reference  : XX----XX

    CASE #3
    =======
    <anvio.bamops.Read object at 0x12b042390>
     ├── start, end : [0, 7)
     ├── cigartuple : [(0, 1), (1, 3), (2, 2), (0, 4)]
     ├── read       : AACC--TTGG
     └── reference  : X---XXXXXX

    CASE #4
    =======
    <anvio.bamops.Read object at 0x12b042400>
     ├── start, end : [0, 11)
     ├── cigartuple : [(0, 1), (2, 1), (0, 1), (2, 1), (0, 1), (2, 1), (0, 1), (2, 1), (0, 2), (1, 1), (0, 1)]
     ├── read       : A-A-C-C-TTGG
     └── reference  : XXXXXXXXXX-X
    """

    def setUp(self):
        self.read_attrs = [
            {
                'query_sequence' : 'AACCTTGG',
                'reference_positions': np.array([0,1,2,3,4,5,6,7]),
                'reference_start': 0,
                'reference_end': 8,
                'cigartuples': np.array([(0,8)]),
                'trim_by': 3,
                'output_tuples_left': np.array([(0,5)]),
                'output_sequence_left': 'CTTGG',
                'output_reference_positions_left': np.array([3,4,5,6,7]),
                'output_tuples_right': np.array([(0,5)]),
                'output_sequence_right': 'AACCT',
                'output_reference_positions_right': np.array([0,1,2,3,4]),
            },
            {
                'query_sequence' : 'AACCTTGG',
                'reference_positions': np.array([0,1,2,3]),
                'reference_start': 0,
                'reference_end': 4,
                'cigartuples': np.array([(0,2), (1,4), (0,2)]),
                'trim_by': 3,
                'output_tuples_left': np.array([(0,1)]),
                'output_sequence_left': 'G',
                'output_reference_positions_left': np.array([3]),
                'output_tuples_right': np.array([(0,1)]),
                'output_sequence_right': 'A',
                'output_reference_positions_right': np.array([0]),
            },
            {
                'query_sequence' : 'AACCTTGG',
                'reference_positions': [0,3,4,5,6],
                'reference_start': 0,
                'reference_end': 7,
                'cigartuples': np.array([(0,1), (1,3), (2,2), (0,4)]),
                'trim_by': 3,
                'output_tuples_left': np.array([(0,4)]),
                'output_sequence_left': 'TTGG',
                'output_reference_positions_left': np.array([3,4,5,6]),
                'output_tuples_right': np.array([(0,1), (1,3), (2,2), (0,1)]),
                'output_sequence_right': 'AACCT',
                'output_reference_positions_right': np.array([0,3]),
            },
            {
                'query_sequence' : 'AACCTTGG',
                'reference_positions': np.array([0,2,4,6,8,9,10]),
                'reference_start': 0,
                'reference_end': 11,
                'cigartuples': np.array([(0,1), (2,1), (0,1), (2,1), (0,1), (2,1), (0,1), (2,1), (0,2), (1,1), (0,1)]),
                'trim_by': 3,
                'output_tuples_left': np.array([(0,1), (2,1), (0,1), (2,1), (0,2), (1,1), (0,1)]),
                'output_sequence_left': 'CCTTGG',
                'output_reference_positions_left': np.array([4,6,8,9,10]),
                'output_tuples_right': np.array([(0,1), (2,1), (0,1), (2,1), (0,1), (2,1), (0,1)]),
                'output_sequence_right': 'AACC',
                'output_reference_positions_right': np.array([0,2,4,6]),
            },
        ]

        self.reads = {}
        for i, attributes in enumerate(self.read_attrs):
            pysam_read = FakePySamAlignedSegment(**attributes)
            self.reads[i] = Read(pysam_read)


    def test_trim(self):
        for target in range(len(self.reads)):
            print("\n\nWorking on target %d:\n%s" % (target, self.reads[target]))
            attrs = self.read_attrs[target]

            print("Trimming left ...")
            read = copy.deepcopy(self.reads[target])
            read.trim(attrs['trim_by'], side='left')
            print("After trimming left by %d:\n%s" % (attrs['trim_by'], read))
            self.assertEqual(read.cigartuples.tolist(), attrs['output_tuples_left'].tolist())

            print("Trimming right ...")
            read = copy.deepcopy(self.reads[target])
            read.trim(attrs['trim_by'], side='right')
            print("After trimming right by %d:\n%s" % (attrs['trim_by'], read))
            self.assertEqual(read.cigartuples.tolist(), attrs['output_tuples_right'].tolist())


class TestCaseStudy1(unittest.TestCase):
    def test_case_study_1(self):
        d = {'gene_overlap_start': 9460, 'gene_overlap_end': 9481, 'cigartuples': [(0, 20), (2, 3), (0, 14), (1, 1), (0, 19), (1, 2), (0, 45)], 'reference_start': 9460, 'reference_end': 9561, 'reference_positions': [9460, 9461, 9462, 9463, 9464, 9465, 9466, 9467, 9468, 9469, 9470, 9471, 9472, 9473, 9474, 9475, 9476, 9477, 9478, 9479, 9483, 9484, 9485, 9486, 9487, 9488, 9489, 9490, 9491, 9492, 9493, 9494, 9495, 9496, 9497, 9498, 9499, 9500, 9501, 9502, 9503, 9504, 9505, 9506, 9507, 9508, 9509, 9510, 9511, 9512, 9513, 9514, 9515, 9516, 9517, 9518, 9519, 9520, 9521, 9522, 9523, 9524, 9525, 9526, 9527, 9528, 9529, 9530, 9531, 9532, 9533, 9534, 9535, 9536, 9537, 9538, 9539, 9540, 9541, 9542, 9543, 9544, 9545, 9546, 9547, 9548, 9549, 9550, 9551, 9552, 9553, 9554, 9555, 9556, 9557, 9558, 9559, 9560], 'query_sequence': 'ACTAGTTTCTTACCTCTATAATTCATAGAGAAAGAAAAATTTAATATGCGCCAATTTTTAAAAAAATTGGTGCCTATTTTTTTAACCAAAATTCTAATATA'}
        read = Read(FakePySamAlignedSegment(**d))

        read2 = copy.copy(read)
        read2.trim(21, side='left')

        self.assertEqual(read2.cigartuples.tolist(), np.array([(0, 14), (1, 1), (0, 19), (1, 2), (0, 45)]).tolist())


if __name__ == '__main__':
    unittest.main()

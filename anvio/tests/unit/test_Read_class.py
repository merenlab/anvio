#! /usr/bin/env python
# -*- coding: utf-8
"""Unit tests for the Read class."""

import copy
import unittest

import anvio

from anvio.sequence import Read
from anvio.errors import ConfigError

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

    def get_reference_positions(self):
        return self.reference_positions


class TestRead(unittest.TestCase):
    """
    CASE #1
    =======
    A A C C T T G G
    A C T G T C T G A C T G = reference
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

    def setUp(self):
        self.read_attrs = [
            {
                'query_sequence' : 'AACCTTGG',
                'reference_positions': [0,1,2,3,4,5,6,7],
                'reference_start': 0,
                'reference_end': 7,
                'cigartuples': [(0,8)],
                'trim_by': 3,
                'output_tuples_left': [(0,5)],
                'output_sequence_left': 'CTTGG',
                'output_reference_positions_left': [3,4,5,6,7],
                'output_tuples_right': [(0,5)],
                'output_sequence_right': 'AACCT',
                'output_reference_positions_right': [0,1,2,3,4],
            },
            {
                'query_sequence' : 'AACCTTGG',
                'reference_positions': [0,1,2,3],
                'reference_start': 0,
                'reference_end': 3,
                'cigartuples': [(0,2), (1,4), (0,2)],
                'trim_by': 3,
                'output_tuples_left': [(0,1)],
                'output_sequence_left': 'G',
                'output_reference_positions_left': [3],
                'output_tuples_right': [(0,1)],
                'output_sequence_right': 'A',
                'output_reference_positions_right': [0],
            },
            {
                'query_sequence' : 'AACCTTGG',
                'reference_positions': [0,3,4,5,6],
                'reference_start': 0,
                'reference_end': 6,
                'cigartuples': [(0,1), (1,3), (2,2), (0,4)],
                'trim_by': 3,
                'output_tuples_left': [(0,4)],
                'output_sequence_left': 'TTGG',
                'output_reference_positions_left': [3,4,5,6],
                'output_tuples_right': [(0,1), (1,3), (2,2), (0,1)],
                'output_sequence_right': 'AACCT',
                'output_reference_positions_right': [0,3],
            },
            {
                'query_sequence' : 'AACCTTGG',
                'reference_positions': [0,2,4,6,8,9,10],
                'reference_start': 0,
                'reference_end': 10,
                'cigartuples': [(0,1), (2,1), (0,1), (2,1), (0,1), (2,1), (0,1), (2,1), (0,2), (1,1), (0,1)],
                'trim_by': 3,
                'output_tuples_left': [(0,1), (2,1), (0,1), (2,1), (0,2), (1,1), (0,1)],
                'output_sequence_left': 'CCTTGG',
                'output_reference_positions_left': [4,6,8,9,10],
                'output_tuples_right': [(0,1), (2,1), (0,1), (2,1), (0,1), (2,1), (0,1)],
                'output_sequence_right': 'AACC',
                'output_reference_positions_right': [0,2,4,6],
            },
        ]

        self.reads = {}
        for i, attributes in enumerate(self.read_attrs):
            pysam_read = FakePySamAlignedSegment(**attributes)
            self.reads[i] = Read(pysam_read)


    def test_trim(self):
        print('=========================')
        for target in range(len(self.reads)):
            attrs = self.read_attrs[target]

            read = copy.copy(self.reads[target])
            read.trim(attrs['trim_by'], side='left')
            self.assertEqual(read.cigartuples, attrs['output_tuples_left'])

            read = copy.copy(self.reads[target])
            read.trim(attrs['trim_by'], side='right')
            self.assertEqual(read.cigartuples, attrs['output_tuples_right'])


class TestCaseStudy1(unittest.TestCase):
    def test_case_study_1(self):
        d = {'gene_overlap_start': 9460, 'gene_overlap_end': 9481, 'cigartuples': [(0, 20), (2, 3), (0, 14), (1, 1), (0, 19), (1, 2), (0, 45)], 'reference_start': 9460, 'reference_end': 9561, 'reference_positions': [9460, 9461, 9462, 9463, 9464, 9465, 9466, 9467, 9468, 9469, 9470, 9471, 9472, 9473, 9474, 9475, 9476, 9477, 9478, 9479, 9483, 9484, 9485, 9486, 9487, 9488, 9489, 9490, 9491, 9492, 9493, 9494, 9495, 9496, 9497, 9498, 9499, 9500, 9501, 9502, 9503, 9504, 9505, 9506, 9507, 9508, 9509, 9510, 9511, 9512, 9513, 9514, 9515, 9516, 9517, 9518, 9519, 9520, 9521, 9522, 9523, 9524, 9525, 9526, 9527, 9528, 9529, 9530, 9531, 9532, 9533, 9534, 9535, 9536, 9537, 9538, 9539, 9540, 9541, 9542, 9543, 9544, 9545, 9546, 9547, 9548, 9549, 9550, 9551, 9552, 9553, 9554, 9555, 9556, 9557, 9558, 9559, 9560], 'query_sequence': 'ACTAGTTTCTTACCTCTATAATTCATAGAGAAAGAAAAATTTAATATGCGCCAATTTTTAAAAAAATTGGTGCCTATTTTTTTAACCAAAATTCTAATATA'}
        read = Read(FakePySamAlignedSegment(**d))

        print('=========================')
        read2 = copy.copy(read)
        read2.trim(21, side='left')

        self.assertEqual(read2.cigartuples, [(0, 14), (1, 1), (0, 19), (1, 2), (0, 45)])


if __name__ == '__main__':
    unittest.main()

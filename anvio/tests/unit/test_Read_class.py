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
        target = 0

        read = copy.copy(self.reads[target])
        attrs = self.read_attrs[target]

        read.trim(attrs['trim_by'], side='left')
        self.assertEqual(read.cigartuples, attrs['output_tuples_left'])

        read.trim(attrs['trim_by'], side='right')
        self.assertEqual(read.cigartuples, attrs['output_tuples_right'])



if __name__ == '__main__':
    unittest.main()

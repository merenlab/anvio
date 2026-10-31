# pylint: disable=line-too-long
"""Unit tests for bamops.Read.iterate_deletions_starting_in: deletions are taken from the untrimmed alignment, so a
deletion that crosses a split boundary is recorded once, whole, in the split where it starts."""

import random
import unittest

import pysam

import anvio

from anvio.bamops import Read

__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__license__ = "GPL 3.0"
__version__ = anvio.__version__


HEADER = pysam.AlignmentHeader.from_dict({'HD': {'VN': '1.6', 'SO': 'coordinate'}, 'SQ': [{'SN': 'c1', 'LN': 100000}]})


def aligned_segment(start, cigar):
    a = pysam.AlignedSegment(HEADER)
    a.query_name = 'r'
    a.reference_id = 0
    a.reference_start = start
    a.cigarstring = cigar
    a.query_sequence = 'A' * a.infer_query_length()
    a.query_qualities = pysam.qualitystring_to_array('I' * a.infer_query_length())
    a.mapping_quality = 60
    return a


def deletions_from_blocks(read):
    """What the profiler recorded before: deletion blocks of the vectorized read"""

    read.vectorize()
    return [(int(b[0, 0]), b.shape[0]) for b in read.iterate_blocks_by_mapping_type(mapping_type=2)]


class DeletionsAcrossRegions(unittest.TestCase):
    def test_deletion_inside_the_region(self):
        r = Read(aligned_segment(1000, '500M300D500M'))
        self.assertEqual(list(r.iterate_deletions_starting_in(0, 20000)), [(1500, 300)])


    def test_deletion_crossing_the_region_end_is_whole_in_the_region_where_it_starts(self):
        # 1,500 bp deletion from 19,500 to 21,000 across a boundary at 20,000
        r = Read(aligned_segment(18000, '1500M1500D1500M'))
        self.assertEqual(list(r.iterate_deletions_starting_in(0, 20000)), [(19500, 1500)])
        self.assertEqual(list(r.iterate_deletions_starting_in(20000, 40000)), [])


    def test_trimming_does_not_change_the_deletions(self):
        a = aligned_segment(18000, '1500M1500D1500M')
        r = Read(a)
        r.trim(trim_by=r.reference_end - 20000, side='right')
        # the trimmed read holds only the first 500 bp of the deletion ...
        self.assertEqual(deletions_from_blocks(r), [(19500, 500)])
        # ... but the deletion is reported whole
        self.assertEqual(list(r.iterate_deletions_starting_in(0, 20000)), [(19500, 1500)])


    def test_deletion_longer_than_a_region(self):
        r = Read(aligned_segment(15000, '1000M25000D1000M'))
        self.assertEqual(list(r.iterate_deletions_starting_in(0, 20000)), [(16000, 25000)])
        self.assertEqual(list(r.iterate_deletions_starting_in(20000, 40000)), [])
        self.assertEqual(list(r.iterate_deletions_starting_in(40000, 60000)), [])


    def test_adjacent_d_and_n_are_one_deletion_and_insertions_split_them(self):
        r = Read(aligned_segment(0, '100M50D20N100M60D5I40D100M'))
        self.assertEqual(list(r.iterate_deletions_starting_in(0, 1000)), [(100, 70), (270, 60), (330, 40)])


    def test_same_as_vectorized_blocks_when_no_deletion_crosses_the_region(self):
        rng = random.Random(1)
        for _ in range(300):
            ops = []
            for _ in range(rng.randint(1, 12)):
                ops.append(f"{rng.randint(1, 60)}M")
                ops.append(f"{rng.randint(1, 80)}{rng.choice('DDDNI')}")
            ops.append(f"{rng.randint(1, 60)}M")
            if rng.random() < 0.5:
                ops.insert(0, f"{rng.randint(1, 50)}S")
            r = Read(aligned_segment(rng.randint(0, 5000), ''.join(ops)))
            expected = deletions_from_blocks(r)
            self.assertEqual(list(r.iterate_deletions_starting_in(0, 100000)), expected, msg=''.join(ops))


if __name__ == '__main__':
    unittest.main()

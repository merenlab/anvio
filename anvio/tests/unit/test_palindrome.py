# -*- coding: utf-8
# pylint: disable=line-too-long
"""Tests for Palindromes"""

import anvio
from anvio.sequencefeatures import Palindromes

import unittest
import argparse

__copyright__ = "Copyleft 2015-2019, the Meren Lab (http://merenlab.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "Evan Kiefl"
__email__ = "kiefl.evan@gmail.com"


class TestFindPalindrome(unittest.TestCase):
    def setUp(self):
        self.nt_30 = 'GNAAANCNNTTTTNTAGAAGNCCAAGTGNN'
        self.nt_10 = 'AACTGGAGCT'


    def test_fixedlength_mismatchless(self):
        seq = self.seq_with_palindromes(
            template = self.nt_30,
            start_stops = [
                (5, 15, 18, 28),
            ],
            palindrome_seqs = [
                (self.nt_10, anvio.utils.rev_comp(self.nt_10)),
            ],
        )

        p = Palindromes(argparse.Namespace(
            min_palindrome_length = 10,
        ))
        palindromes = p._find_numba(seq, coords_only=True)

        for_start, for_stop, rev_start, rev_stop = palindromes[0]

        self.assertEqual(for_start, 5)
        self.assertEqual(for_stop, 15)
        self.assertEqual(rev_start, 18)
        self.assertEqual(rev_stop, 28)

        self.assertEqual(seq[for_start:for_stop], seq[5:15])
        self.assertEqual(seq[rev_start:rev_stop], seq[18:28])


    def test_fixedlength_mismatchless_touching(self):
        seq = self.seq_with_palindromes(
            template = self.nt_30,
            start_stops = [
                (5, 15, 15, 25),
            ],
            palindrome_seqs = [
                (self.nt_10, anvio.utils.rev_comp(self.nt_10)),
            ],
        )

        p = Palindromes(argparse.Namespace(
            min_palindrome_length = 10,
        ))
        palindromes = p._find_numba(seq, coords_only=True)

        for_start, for_stop, rev_start, rev_stop = palindromes[0]

        self.assertEqual(for_start, 5)
        self.assertEqual(for_stop, 15)
        self.assertEqual(rev_start, 15)
        self.assertEqual(rev_stop, 25)

        self.assertEqual(seq[for_start:for_stop], seq[5:15])
        self.assertEqual(seq[rev_start:rev_stop], seq[15:25])


    def test_fixedlength_mismatchless_1_apart(self):
        delta = 1
        seq = self.seq_with_palindromes(
            template = self.nt_30,
            start_stops = [
                (4, 14, 14+delta, 24+delta),
            ],
            palindrome_seqs = [
                (self.nt_10, anvio.utils.rev_comp(self.nt_10)),
            ],
        )

        p = Palindromes(argparse.Namespace(
            min_palindrome_length = 10,
        ))
        palindromes = p._find_numba(seq, coords_only=True)

        for_start, for_stop, rev_start, rev_stop = palindromes[0]

        self.assertEqual(for_start, 4)
        self.assertEqual(for_stop, 14)
        self.assertEqual(rev_start, 14+delta)
        self.assertEqual(rev_stop, 24+delta)

        self.assertEqual(seq[for_start:for_stop], seq[4:14])
        self.assertEqual(seq[rev_start:rev_stop], seq[14+delta:24+delta])


    def test_fixedlength_mismatchless_multimatch(self):
        seq = self.seq_with_palindromes(
            template = self.nt_30*2,
            start_stops = [
                (5, 15, 18, 28),
                (5, 15, 38, 48),
            ],
            palindrome_seqs = [
                (self.nt_10, anvio.utils.rev_comp(self.nt_10)),
                (self.nt_10, anvio.utils.rev_comp(self.nt_10)),
            ],
        )

        p = Palindromes(argparse.Namespace(
            min_palindrome_length = 10,
        ))
        palindromes = p._find_numba(seq, coords_only=True)

        for_start, for_stop, rev_start, rev_stop = palindromes[0]
        self.assertEqual(for_start, 5)
        self.assertEqual(for_stop, 15)
        self.assertEqual(rev_start, 38)
        self.assertEqual(rev_stop, 48)
        self.assertEqual(seq[for_start:for_stop], seq[5:15])
        self.assertEqual(seq[rev_start:rev_stop], seq[38:48])

        for_start, for_stop, rev_start, rev_stop = palindromes[1]
        self.assertEqual(for_start, 5)
        self.assertEqual(for_stop, 15)
        self.assertEqual(rev_start, 18)
        self.assertEqual(rev_stop, 28)
        self.assertEqual(seq[for_start:for_stop], seq[5:15])
        self.assertEqual(seq[rev_start:rev_stop], seq[18:28])


    def test_fixedlength_mismatch(self):
        pal_for = 'AACTGGAGCT'
        #          ||| ||| ||
        pal_rev = 'AACAGGAACT'
        pal_rev = anvio.utils.rev_comp(pal_rev)

        seq = self.seq_with_palindromes(
            template = self.nt_30,
            start_stops = [
                (5, 15, 18, 28),
            ],
            palindrome_seqs = [
                (pal_for, pal_rev),
            ],
        )

        # No palindrome should be found
        MISMATCH_TOL = 1
        p = Palindromes(argparse.Namespace(
            min_palindrome_length = 10,
            max_num_mismatches=MISMATCH_TOL,
        ))
        palindromes = p._find_numba(seq, coords_only=True)
        self.assertEqual(palindromes, [])

        # Palindrome should be found
        MISMATCH_TOL = 2
        p = Palindromes(argparse.Namespace(
            min_palindrome_length = 10,
            max_num_mismatches=MISMATCH_TOL,
        ))
        palindromes = p._find_numba(seq, coords_only=True)
        for_start, for_stop, rev_start, rev_stop = palindromes[0]
        self.assertEqual(for_start, 5)
        self.assertEqual(for_stop, 15)
        self.assertEqual(rev_start, 18)
        self.assertEqual(rev_stop, 28)
        self.assertEqual(seq[for_start:for_stop], seq[5:15])
        self.assertEqual(seq[rev_start:rev_stop], seq[18:28])

        # Palindrome should be found
        MISMATCH_TOL = 3
        p = Palindromes(argparse.Namespace(
            min_palindrome_length = 10,
            max_num_mismatches=MISMATCH_TOL,
        ))
        palindromes = p._find_numba(seq, coords_only=True)
        for_start, for_stop, rev_start, rev_stop = palindromes[0]
        self.assertEqual(for_start, 5)
        self.assertEqual(for_stop, 15)
        self.assertEqual(rev_start, 18)
        self.assertEqual(rev_stop, 28)
        self.assertEqual(seq[for_start:for_stop], seq[5:15])
        self.assertEqual(seq[rev_start:rev_stop], seq[18:28])


    def test_variedlength_mismatchless(self):
        for m in range(4, 10):
            seq = self.seq_with_palindromes(
                template = self.nt_30,
                start_stops = [
                    (5, 15, 18, 28),
                ],
                palindrome_seqs = [
                    (self.nt_10, anvio.utils.rev_comp(self.nt_10)),
                ],
            )
            p = Palindromes(argparse.Namespace(
                min_palindrome_length = m,
            ))
            palindromes = p._find_numba(seq, coords_only=True)

            for_start, for_stop, rev_start, rev_stop = palindromes[0]

            self.assertEqual(for_start, 5)
            self.assertEqual(for_stop, 15)
            self.assertEqual(rev_start, 18)
            self.assertEqual(rev_stop, 28)

            self.assertEqual(seq[for_start:for_stop], seq[5:15])
            self.assertEqual(seq[rev_start:rev_stop], seq[18:28])


    def test_variedlength_mismatchless_touching(self):
        seq = self.seq_with_palindromes(
            template = self.nt_30,
            start_stops = [
                (5, 15, 15, 25),
            ],
            palindrome_seqs = [
                (self.nt_10, anvio.utils.rev_comp(self.nt_10)),
            ],
        )

        p = Palindromes(argparse.Namespace(
            min_palindrome_length = 4,
        ))
        palindromes = p._find_numba(seq, coords_only=True)

        for_start, for_stop, rev_start, rev_stop = palindromes[0]

        self.assertEqual(for_start, 5)
        self.assertEqual(for_stop, 15)
        self.assertEqual(rev_start, 15)
        self.assertEqual(rev_stop, 25)

        self.assertEqual(seq[for_start:for_stop], seq[5:15])
        self.assertEqual(seq[rev_start:rev_stop], seq[15:25])


    def test_variedlength_mismatchless_multimatch(self):
        seq = self.seq_with_palindromes(
            template = self.nt_30*2,
            start_stops = [
                (5, 15, 18, 28),
                (5, 15, 38, 48),
            ],
            palindrome_seqs = [
                (self.nt_10, anvio.utils.rev_comp(self.nt_10)),
                (self.nt_10, anvio.utils.rev_comp(self.nt_10)),
            ],
        )

        p = Palindromes(argparse.Namespace(
            min_palindrome_length = 6,
        ))
        palindromes = p._find_numba(seq, coords_only=True)

        for_start, for_stop, rev_start, rev_stop = palindromes[0]
        self.assertEqual(for_start, 5)
        self.assertEqual(for_stop, 15)
        self.assertEqual(rev_start, 38)
        self.assertEqual(rev_stop, 48)
        self.assertEqual(seq[for_start:for_stop], seq[5:15])
        self.assertEqual(seq[rev_start:rev_stop], seq[38:48])

        for_start, for_stop, rev_start, rev_stop = palindromes[1]
        self.assertEqual(for_start, 5)
        self.assertEqual(for_stop, 15)
        self.assertEqual(rev_start, 18)
        self.assertEqual(rev_stop, 28)
        self.assertEqual(seq[for_start:for_stop], seq[5:15])
        self.assertEqual(seq[rev_start:rev_stop], seq[18:28])


    def test_split(self):
        pal_for = 'AACTGGAGCTCGAACTG'
        #          |||||||| |||||||| 
        pal_rev = 'AACTGGAGGTCGAACTG'
        pal_rev = anvio.utils.rev_comp(pal_rev)

        x0, x1 = 5, 30

        seq = self.seq_with_palindromes(
            template = self.nt_30*2,
            start_stops = [
                (x0, x0+len(pal_for), x1, x1+len(pal_for)),
            ],
            palindrome_seqs = [
                (pal_for, pal_rev),
            ],
        )

        # One long palindrome should be found
        MISMATCH_TOL = 1
        p = Palindromes(argparse.Namespace(
            min_palindrome_length = 8,
            max_num_mismatches=MISMATCH_TOL,
        ))
        palindromes = p._find_numba(seq, coords_only=True)

        # Two shorter palindrome should be found
        MISMATCH_TOL = 0
        p = Palindromes(argparse.Namespace(
            min_palindrome_length = 7,
            max_num_mismatches=MISMATCH_TOL,
        ))
        palindromes = p._find_numba(seq, coords_only=True)

        for_start, for_stop, rev_start, rev_stop = palindromes[0]
        self.assertEqual(for_start, x0)
        self.assertEqual(for_stop, 13)
        self.assertEqual(rev_start, 39)
        self.assertEqual(rev_stop, 47)
        self.assertEqual(seq[for_start:for_stop], seq[5:13])
        self.assertEqual(seq[rev_start:rev_stop], seq[39:47])

        for_start, for_stop, rev_start, rev_stop = palindromes[1]
        self.assertEqual(for_start, 14)
        self.assertEqual(for_stop, 22)
        self.assertEqual(rev_start, 30)
        self.assertEqual(rev_stop, 38)
        self.assertEqual(seq[for_start:for_stop], seq[14:22])
        self.assertEqual(seq[rev_start:rev_stop], seq[30:38])


    def test_min_distance(self):
        delta = 3
        seq = self.seq_with_palindromes(
            template = self.nt_30*2,
            start_stops = [
                (4, 14, 14+delta, 24+delta),
            ],
            palindrome_seqs = [
                (self.nt_10, anvio.utils.rev_comp(self.nt_10)),
            ],
        )

        for D in range(10):
            expected = 0 if D > delta else 1
            p = Palindromes(argparse.Namespace(
                min_palindrome_length = 10,
                min_distance=D,
            ))
            palindromes = p._find_numba(seq, coords_only=True)

            self.assertEqual(len(palindromes), expected)


    def test_min_distance2(self):
        seq = 'TTTCAAGGGGGGGGGTTGAAA'
        for d in range(13):
            p = Palindromes(argparse.Namespace(
                min_palindrome_length = 4,
                min_distance = d,
            ))
            palindromes = p._find_numba(seq, coords_only=True)

            if d <= 9:
                # the gap is 9, so this should produce 1 palindrome
                self.assertEqual(len(palindromes), 1)
            else:
                # the gap is 9, so this should not produce palindrome
                self.assertEqual(len(palindromes), 0)


    def test_triple_inversion(self):
        seq = 'TTTCAGGTGAAACTGAA'

        p = Palindromes(argparse.Namespace(
            min_palindrome_length = 5,
            min_distance = 3,
        ))
        palindromes = p._find_numba(seq, coords_only=True)
        self.assertEqual(palindromes, [(1,6,12,17)])

        p = Palindromes(argparse.Namespace(
            min_palindrome_length = 5,
            min_distance = 2,
        ))
        palindromes = p._find_numba(seq, coords_only=True)
        self.assertEqual(palindromes, [(0, 5, 7, 12), (1, 6, 12, 17)])


    def seq_with_palindromes(self, template, start_stops, palindrome_seqs):
        seq = list(template)

        for start_stop, palindrome_seq in zip(start_stops, palindrome_seqs):
            for_start, for_stop, rev_start, rev_stop = start_stop
            subseq, revseq = palindrome_seq

            seq[for_start:for_stop] = subseq
            seq[rev_start:rev_stop] = revseq

        return ''.join(seq)


if __name__ == '__main__':
    unittest.main()

# -*- coding: utf-8
# pylint: disable=line-too-long
"""Tests for FindPalindrome"""

import anvio
from anvio.sequencefeatures import FindPalindrome

import unittest

__author__ = "Developers of anvi'o (see AUTHORS.txt)"
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

        p = FindPalindrome(min_len=10)
        palindromes = p.find(seq)

        for_start, for_stop, rev_start, rev_stop = palindromes[0]

        self.assertEqual(for_start, 5)
        self.assertEqual(for_stop, 15)
        self.assertEqual(rev_start, 18)
        self.assertEqual(rev_stop, 28)

        self.assertEqual(seq[for_start:for_stop], seq[5:15])
        self.assertEqual(seq[rev_start:rev_stop], seq[18:28])


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

        p = FindPalindrome(min_len=10)
        palindromes = p.find(seq)

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


    #def test_variedlength_mismatchless(self):
    #    seq = 'GNAAANCNNTTTTNTAGAAGNCCAAGTGNN'

    #    for_start_expected = 5
    #    for_stop_expected = 15

    #    rev_start_expected = 18
    #    rev_stop_expected = 28

    #    for_pal = self.nt_10
    #    rev_pal = anvio.utils.rev_comp(for_pal)

    #    # Insert palindrome into seq
    #    for_pal = list(for_pal)
    #    rev_pal = list(rev_pal)
    #    seq = list(seq)
    #    seq[for_start_expected:for_stop_expected] = for_pal
    #    seq[rev_start_expected:rev_stop_expected] = rev_pal
    #    seq = ''.join(seq)

    #    p = FindPalindrome(seq, min_len=3)
    #    palindromes = p.find()

    #    for_start, for_stop, rev_start, rev_stop = palindromes[0]

    #    self.assertEqual(for_start, for_start_expected)
    #    self.assertEqual(for_stop, for_stop_expected)
    #    self.assertEqual(rev_start, rev_start_expected)
    #    self.assertEqual(rev_stop, rev_stop_expected)

    #    self.assertEqual(seq[for_start:for_stop], seq[for_start_expected:for_stop_expected])
    #    self.assertEqual(seq[rev_start:rev_stop], seq[rev_start_expected:rev_stop_expected])


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

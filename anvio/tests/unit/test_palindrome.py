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
        pass


    def test_fixedlength_mismatchless(self):
        seq = 'GNAAANCNNTTTTNTAGAAGNCCAAGTGNN'

        for_start_expected = 5
        for_stop_expected = 15

        rev_start_expected = 18
        rev_stop_expected = 28

        for_pal = 'AACTGGAGCT'
        rev_pal = anvio.utils.rev_comp(for_pal)

        # Insert palindrome into seq
        for_pal = list(for_pal)
        rev_pal = list(rev_pal)
        seq = list(seq)
        seq[for_start_expected:for_stop_expected] = for_pal
        seq[rev_start_expected:rev_stop_expected] = rev_pal
        seq = ''.join(seq)

        p = FindPalindrome(seq, min_len=len(for_pal))
        palindromes = p.find()

        for_start, for_stop, rev_start, rev_stop = palindromes[0]

        self.assertEqual(for_start, for_start_expected)
        self.assertEqual(for_stop, for_stop_expected)
        self.assertEqual(rev_start, rev_start_expected)
        self.assertEqual(rev_stop, rev_stop_expected)

        self.assertEqual(seq[for_start:for_stop], seq[for_start_expected:for_stop_expected])
        self.assertEqual(seq[rev_start:rev_stop], seq[rev_start_expected:rev_stop_expected])


    def test_mismatchless(self):
        seq = 'GNAAANCNNTTTTNTAGAAGNCCAAGTGNN'

        for_start_expected = 5
        for_stop_expected = 15

        rev_start_expected = 18
        rev_stop_expected = 28

        for_pal = 'AACTGGAGCT'
        rev_pal = anvio.utils.rev_comp(for_pal)

        # Insert palindrome into seq
        for_pal = list(for_pal)
        rev_pal = list(rev_pal)
        seq = list(seq)
        seq[for_start_expected:for_stop_expected] = for_pal
        seq[rev_start_expected:rev_stop_expected] = rev_pal
        seq = ''.join(seq)

        p = FindPalindrome(seq, min_len=8)
        palindromes = p.find()

        for_start, for_stop, rev_start, rev_stop = palindromes[0]

        self.assertEqual(for_start, for_start_expected)
        self.assertEqual(for_stop, for_stop_expected)
        self.assertEqual(rev_start, rev_start_expected)
        self.assertEqual(rev_stop, rev_stop_expected)

        self.assertEqual(seq[for_start:for_stop], seq[for_start_expected:for_stop_expected])
        self.assertEqual(seq[rev_start:rev_stop], seq[rev_start_expected:rev_stop_expected])




if __name__ == '__main__':
    unittest.main()

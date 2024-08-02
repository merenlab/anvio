# -*- coding: utf-8
# pylint: disable=line-too-long
"""
Unit tests for the split_fasta function.
"""

import unittest as ut

import os

import anvio
from anvio.errors import FilesNPathsError
from anvio.utils import split_fasta
from anvio.fastalib import ReadFasta

__copyright__ = "Copyleft 2015-2019, the Meren Lab (http://merenlab.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "Ryan Moore"
__email__ = "moorer@udel.edu"

class SplitFastaTestCase(ut.TestCase):
    def setUp(self):
        self.this_dir = os.path.dirname(os.path.realpath(__file__))
        self.test_files = os.path.join(self.this_dir, "test_files")
        self.not_fasta = os.path.join(self.test_files, "not_a_fasta.txt")
        self.empty_fasta = os.path.join(self.test_files, "empty.fasta")
        self.single_seq_fasta = os.path.join(self.test_files, "one_sequence.fasta")
        self.five_seq_fasta = os.path.join(self.test_files, "five_sequences.fasta")

    def test_non_existent_file_raises_error(self):
        self.assertRaises(FilesNPathsError, split_fasta, "arstoien.fasta")

    def test_non_fasta_file_raises_error(self):
        self.assertRaises(FilesNPathsError, split_fasta, self.not_fasta)

    def test_empty_fasta_file_raises_error(self):
        self.assertRaises(FilesNPathsError, split_fasta, self.empty_fasta)

    def test_single_fasta_gives_one_split(self):
        out_files = split_fasta(self.single_seq_fasta)

        expected_out_file = os.path.join(self.test_files, f'{self.single_seq_fasta}.0')

        self.assertEqual(out_files, [expected_out_file])

        self.assertTrue(os.path.exists(expected_out_file))

        fasta = ReadFasta(expected_out_file)

        self.assertEqual(fasta.ids, ['seq1 apple'])
        self.assertEqual(fasta.sequences, ['AA'])

        fasta.close()

        os.remove(expected_out_file)

    def test_fasta_splitting(self):
        parts = 2
        expected_out_files = [os.path.join(self.test_files, f'{self.five_seq_fasta}.{i}') for i in range(parts)]

        out_files = split_fasta(self.five_seq_fasta, parts=parts)

        self.assertEqual(out_files, expected_out_files)

        fasta = ReadFasta(out_files[0])
        self.assertEqual(fasta.ids, ['seq1 apple', 'seq2 banana'])
        self.assertEqual(fasta.sequences, ['AA', 'ACAC'])
        fasta.close()

        fasta = ReadFasta(out_files[1])
        self.assertEqual(fasta.ids, ['seq3 cat', 'seq4 dog', 'seq5 extra'])
        self.assertEqual(fasta.sequences, ['ACTACT', 'ACTGACTG', 'ACTGAACTGA'])
        fasta.close()

        for f in out_files:
            os.remove(f)

    def test_more_parts_than_sequences(self):
        parts = 10
        num_sequences = 5
        expected_out_files = [os.path.join(self.test_files, f'{self.five_seq_fasta}.{i}') for i in range(num_sequences)]

        out_files = split_fasta(self.five_seq_fasta, parts=parts)

        self.assertEqual(out_files, expected_out_files)

        for f in out_files:
            os.remove(f)

    def test_custom_prefix(self):
        parts = 1
        file_name_prefix = 'silly'

        out_files = split_fasta(self.five_seq_fasta, parts=parts, file_name_prefix=file_name_prefix, output_dir=self.this_dir)
        expected_out_files = [os.path.join(self.this_dir, 'silly.0')]

        self.assertEqual(out_files, expected_out_files)

        for f in out_files:
            os.remove(f)

    def test_shuffle_mode(self):
        parts = 2

        out_files = split_fasta(self.five_seq_fasta, parts=parts, shuffle=True)

        fasta = ReadFasta(out_files[0])
        self.assertEqual(fasta.ids, ['seq1 apple', 'seq3 cat', 'seq5 extra'])
        self.assertEqual(fasta.sequences, ['AA', 'ACTACT', 'ACTGAACTGA'])
        fasta.close()

        fasta = ReadFasta(out_files[1])
        self.assertEqual(fasta.ids, ['seq2 banana', 'seq4 dog'])
        self.assertEqual(fasta.sequences, ['ACAC', 'ACTGACTG'])
        fasta.close()

        for f in out_files:
            os.remove(f)






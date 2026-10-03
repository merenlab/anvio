"""CLI defaults for genome similarity calculations."""

import sys
import unittest
from unittest import mock

from anvio.cli.compute_genome_similarity import get_args


class ComputeGenomeSimilarityCLITestCase(unittest.TestCase):
    def test_fastani_is_the_default_program(self):
        argv = [
            'anvi-compute-genome-similarity',
            '-e', 'external-genomes.txt',
            '-o', 'similarity-output',
        ]
        with mock.patch.object(sys, 'argv', argv):
            args = get_args()

        self.assertEqual(args.program, 'fastANI')

    def test_pyani_remains_selectable(self):
        argv = [
            'anvi-compute-genome-similarity',
            '-e', 'external-genomes.txt',
            '-o', 'similarity-output',
            '--program', 'pyANI',
        ]
        with mock.patch.object(sys, 'argv', argv):
            args = get_args()

        self.assertEqual(args.program, 'pyANI')


if __name__ == '__main__':
    unittest.main()

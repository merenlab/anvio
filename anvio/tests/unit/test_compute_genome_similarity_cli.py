"""CLI defaults for genome similarity calculations."""

import sys
import unittest
from unittest import mock

from anvio.cli.compute_genome_similarity import get_args, validate_program_options
from anvio.errors import ConfigError


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

    def test_pyani_options_are_rejected_for_fastani_when_explicit(self):
        option_values = {
            '--method': 'ANIm',
            '--min-alignment-fraction': '0.8',
            '--significant-alignment-length': '5000',
            '--min-full-percent-identity': '90',
        }
        for program in ('fastANI', 'sourmash'):
            for option, value in option_values.items():
                for spelling in ([option, value], [f'{option}={value}']):
                    with self.subTest(program=program, spelling=spelling):
                        args = mock.Mock(program=program)
                        with self.assertRaises(ConfigError):
                            validate_program_options(args, argv=spelling)

    def test_omitted_pyani_defaults_do_not_reject_fastani_or_sourmash(self):
        for program in ('fastANI', 'sourmash'):
            with self.subTest(program=program):
                validate_program_options(mock.Mock(program=program), argv=[])

    def test_pyani_options_remain_valid_for_pyani(self):
        validate_program_options(mock.Mock(program='pyANI'), argv=['--method=ANIm'])


if __name__ == '__main__':
    unittest.main()

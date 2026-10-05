"""CLI defaults for genome similarity calculations."""

import contextlib
import io
import sys
import unittest
from argparse import Namespace
from unittest import mock

from anvio.cli.compute_genome_similarity import get_args, main, validate_program_options
from anvio.errors import ConfigError


class ComputeGenomeSimilarityCLITestCase(unittest.TestCase):
    def test_help_describes_pyani_plus_cpu_allocation(self):
        output = io.StringIO()
        with mock.patch.object(sys, 'argv', ['anvi-compute-genome-similarity', '--help']):
            with contextlib.redirect_stdout(output):
                with self.assertRaises(SystemExit) as error:
                    get_args()
        self.assertEqual(error.exception.code, 0)
        self.assertIn('CPU affinity', output.getvalue())
        self.assertIn('taskset', output.getvalue())

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

    def test_pyani_backend_defaults_to_plus_and_legacy_remains_selectable(self):
        argv = [
            'anvi-compute-genome-similarity',
            '-e', 'external-genomes.txt',
            '-o', 'similarity-output',
            '--program', 'pyANI',
        ]
        with mock.patch.object(sys, 'argv', argv):
            args = get_args()

        self.assertIsNone(args.ani_backend)
        self.assertIsNone(args.pyani_plus_program)

        argv += ['--ani-backend', 'legacy']
        with mock.patch.object(sys, 'argv', argv):
            args = get_args()

        self.assertEqual(args.ani_backend, 'legacy')

    def test_backend_options_are_declared_once_in_shared_args(self):
        from anvio import A, K

        self.assertEqual(A('ani-backend'), ['--ani-backend'])
        self.assertEqual(K('ani-backend')['default'], None)
        self.assertEqual(A('pyani-plus-program'), ['--pyani-plus-program'])
        self.assertEqual(K('pyani-plus-program')['default'], None)

    def test_pyani_plus_program_is_rejected_for_fastani_before_processing(self):
        args = Namespace(program='fastANI', ani_backend=None, pyani_plus_program='/custom/pyani-plus')
        with mock.patch('anvio.cli.compute_genome_similarity.get_args', return_value=args), \
                mock.patch('anvio.cli.compute_genome_similarity.genomesimilarity.program_class_dictionary') as programs:
            with self.assertRaises(SystemExit):
                main()
        programs.__getitem__.assert_not_called()

    def test_pyani_option_is_rejected_before_processing(self):
        explicit_pyani_options = [
            ['--method', 'ANIb'],
            ['--min-alignment-fraction', '0.0'],
            ['--significant-alignment-length', '5000'],
            ['--min-full-percent-identity', '0.0'],
            ['--meth=ANIb'],
        ]
        for program in ('fastANI', 'sourmash'):
            for option in explicit_pyani_options:
                with self.subTest(program=program, option=option):
                    args = Namespace(program=program)
                    argv = ['anvi-compute-genome-similarity', *option]
                    with mock.patch.object(sys, 'argv', argv), \
                            mock.patch('anvio.cli.compute_genome_similarity.get_args', return_value=args), \
                            mock.patch('anvio.cli.compute_genome_similarity.genomesimilarity.program_class_dictionary') as programs:
                        with self.assertRaises(SystemExit):
                            main()
                    programs.__getitem__.assert_not_called()

    def test_pyani_options_are_rejected_for_fastani_when_explicit(self):
        option_values = {
            '--method': 'ANIb',
            '--min-alignment-fraction': '0.0',
            '--significant-alignment-length': '5000',
            '--min-full-percent-identity': '0.0',
        }
        for program in ('fastANI', 'sourmash'):
            for option, value in option_values.items():
                for spelling in ([option, value], [f'{option}={value}']):
                    with self.subTest(program=program, spelling=spelling):
                        args = mock.Mock(program=program)
                        with self.assertRaises(ConfigError):
                            validate_program_options(args, argv=spelling)

    def test_abbreviated_pyani_options_are_also_detected(self):
        abbreviations = [
            ['--meth=ANIb'],
            ['--min-alignment-fra=0.8'],
            ['--significant-alignment-len=5000'],
            ['--min-full-percent-ident=90'],
        ]
        for argv in abbreviations:
            with self.subTest(argv=argv):
                with self.assertRaises(ConfigError):
                    validate_program_options(mock.Mock(program='fastANI'), argv=argv)

    def test_pyani_option_is_rejected_before_processing(self):
        explicit_pyani_options = [
            ['--method', 'ANIb'],
            ['--min-alignment-fraction', '0.0'],
            ['--significant-alignment-length', '5000'],
            ['--min-full-percent-identity', '0.0'],
            ['--meth=ANIb'],
            ['--min-alignment-fra=0.8'],
            ['--significant-alignment-len=5000'],
            ['--min-full-percent-ident=90'],
        ]
        for program in ('fastANI', 'sourmash'):
            for option in explicit_pyani_options:
                with self.subTest(program=program, option=option):
                    args = Namespace(program=program)
                    argv = ['anvi-compute-genome-similarity', *option]
                    with mock.patch.object(sys, 'argv', argv), \
                            mock.patch('anvio.cli.compute_genome_similarity.get_args', return_value=args), \
                            mock.patch('anvio.cli.compute_genome_similarity.genomesimilarity.program_class_dictionary') as programs:
                        with self.assertRaises(SystemExit):
                            main()
                programs.__getitem__.assert_not_called()

    def test_omitted_pyani_defaults_do_not_reject_fastani_or_sourmash(self):
        for program in ('fastANI', 'sourmash'):
            with self.subTest(program=program):
                validate_program_options(mock.Mock(program=program), argv=[])

    def test_pyani_options_remain_valid_for_pyani(self):
        validate_program_options(mock.Mock(program='pyANI'), argv=['--method=ANIm'])


if __name__ == '__main__':
    unittest.main()

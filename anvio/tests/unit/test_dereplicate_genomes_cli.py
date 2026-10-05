"""CLI options for pyANI-backed dereplication."""

import contextlib
import io
import sys
import unittest
from unittest import mock

from anvio.cli.dereplicate_genomes import get_args


class DereplicateGenomesCLITestCase(unittest.TestCase):
    def test_help_describes_pyani_plus_cpu_allocation(self):
        output = io.StringIO()
        with mock.patch.object(sys, 'argv', ['anvi-dereplicate-genomes', '--help']):
            with contextlib.redirect_stdout(output):
                with self.assertRaises(SystemExit) as error:
                    get_args()
        self.assertEqual(error.exception.code, 0)
        self.assertIn('CPU affinity', output.getvalue())
        self.assertIn('taskset', output.getvalue())

    def test_pyani_plus_backend_options_are_available_for_matrix_import(self):
        argv = [
            'anvi-dereplicate-genomes',
            '--ani-dir', 'ani-results',
            '--program', 'pyANI',
            '--similarity-threshold', '0.99',
            '-o', 'derep-output',
        ]
        with mock.patch.object(sys, 'argv', argv):
            args = get_args()

        self.assertIsNone(args.ani_backend)
        self.assertIsNone(args.pyani_plus_program)

        argv += ['--ani-backend', 'legacy', '--pyani-plus-program', '/opt/ani/bin/pyani-plus']
        with mock.patch.object(sys, 'argv', argv):
            args = get_args()
        self.assertEqual(args.ani_backend, 'legacy')
        self.assertEqual(args.pyani_plus_program, '/opt/ani/bin/pyani-plus')

    def test_incompatible_pyani_options_fail_before_dereplication(self):
        from anvio.cli import dereplicate_genomes

        options = (
            ('--method', '--meth', 'ANIb'),
            ('--min-alignment-fraction', '--min-alignment-fra', '0.25'),
            ('--significant-alignment-length', '--significant-alignment-len', '10'),
            ('--min-full-percent-identity', '--min-full-percent-iden', '20.0'),
        )
        spellings = [
            spelling
            for full_option, abbreviated_option, value in options
            for spelling in (
                [full_option, value],
                [f'{full_option}={value}'],
                [abbreviated_option, value],
                [f'{abbreviated_option}={value}'],
            )
        ]

        for program in (None, 'fastANI', 'sourmash'):
            for spelling in spellings:
                argv = ['anvi-dereplicate-genomes']
                if program is not None:
                    argv += ['--program', program]
                argv += ['--similarity-threshold', '0.99', '-o', 'output', *spelling]
                with self.subTest(program=program or '(default)', spelling=spelling):
                    with mock.patch.object(sys, 'argv', argv), \
                         mock.patch.object(dereplicate_genomes, 'Dereplicate') as dereplicate, \
                         contextlib.redirect_stdout(io.StringIO()):
                        with self.assertRaises(SystemExit):
                            dereplicate_genomes.main()
                    dereplicate.assert_not_called()

    def test_omitted_pyani_options_and_compatible_methods_reach_dereplication(self):
        from anvio.cli import dereplicate_genomes

        cases = (
            (None, []),
            ('fastANI', []),
            ('sourmash', []),
            ('pyANI', []),
            ('pyANI', ['--method', 'ANIb']),
            ('pyANI', ['--method=ANIm']),
            ('pyANI', ['--min-alignment-fraction=0.25', '--min-full-percent-identity=20.0']),
        )
        for program, options in cases:
            argv = ['anvi-dereplicate-genomes']
            if program is not None:
                argv += ['--program', program]
            argv += ['--similarity-threshold', '0.99', '-o', 'output', *options]
            with self.subTest(program=program, options=options):
                with mock.patch.object(sys, 'argv', argv), \
                     mock.patch.object(dereplicate_genomes, 'Dereplicate') as dereplicate:
                    dereplicate_genomes.main()
                dereplicate.assert_called_once()
                args = dereplicate.call_args.args[0]
                if not options:
                    self.assertEqual(args.method, 'ANIb')
                    self.assertEqual(args.min_alignment_fraction, 0.25)
                    self.assertIsNone(args.significant_alignment_length)
                    self.assertEqual(args.min_full_percent_identity, 20.0)
                if program == 'pyANI' and options:
                    self.assertEqual(args.method, 'ANIm' if '--method=ANIm' in options else 'ANIb')
                dereplicate.return_value.process.assert_called_once()

    def test_pyani_computation_options_fail_when_importing_matrices(self):
        from anvio.cli import dereplicate_genomes

        options = (
            ('--method', '--meth', 'ANIb'),
            ('--min-alignment-fraction', '--min-alignment-fra', '0.25'),
            ('--significant-alignment-length', '--significant-alignment-len', '10'),
            ('--min-full-percent-identity', '--min-full-percent-iden', '20.0'),
        )
        spellings = [
            spelling
            for full_option, abbreviated_option, value in options
            for spelling in (
                [full_option, value],
                [f'{full_option}={value}'],
                [abbreviated_option, value],
                [f'{abbreviated_option}={value}'],
            )
        ]
        import_modes = (
            (['--program', 'pyANI', '--ani-dir', 'ani-results'], '--ani-dir'),
            (['--program', 'sourmash', '--mash-dir', 'mash-results'], '--mash-dir'),
        )

        for import_args, import_option in import_modes:
            for spelling in spellings:
                argv = [
                    'anvi-dereplicate-genomes', *import_args,
                    '--similarity-threshold', '0.99', '-o', 'output', *spelling,
                ]
                with self.subTest(import_option=import_option, spelling=spelling):
                    output = io.StringIO()
                    with mock.patch.object(sys, 'argv', argv), \
                         mock.patch.object(dereplicate_genomes, 'Dereplicate') as dereplicate, \
                         contextlib.redirect_stdout(output):
                        with self.assertRaises(SystemExit):
                            dereplicate_genomes.main()
                    dereplicate.assert_not_called()
                    self.assertIn('cannot be used while importing existing', output.getvalue())
                    self.assertIn(import_option, output.getvalue())

    def test_imported_matrices_preserve_omitted_pyani_defaults(self):
        from anvio.cli import dereplicate_genomes

        cases = (
            (['--program', 'pyANI', '--ani-dir', 'ani-results'], 'ani_dir'),
            (['--program', 'sourmash', '--mash-dir', 'mash-results'], 'mash_dir'),
        )
        for import_args, import_attribute in cases:
            argv = [
                'anvi-dereplicate-genomes', *import_args,
                '--similarity-threshold', '0.99', '-o', 'output',
            ]
            with self.subTest(import_attribute=import_attribute):
                with mock.patch.object(sys, 'argv', argv), \
                     mock.patch.object(dereplicate_genomes, 'Dereplicate') as dereplicate:
                    dereplicate_genomes.main()
                dereplicate.assert_called_once()
                args = dereplicate.call_args.args[0]
                self.assertEqual(getattr(args, import_attribute), import_args[-1])
                self.assertEqual(args.method, 'ANIb')
                self.assertEqual(args.min_alignment_fraction, 0.25)
                self.assertIsNone(args.significant_alignment_length)
                self.assertEqual(args.min_full_percent_identity, 20.0)
                dereplicate.return_value.process.assert_called_once()

    def test_direct_api_rejects_irrelevant_pyani_plus_program(self):
        from argparse import Namespace
        from anvio.errors import ConfigError
        from anvio.genomesimilarity import Dereplicate

        with self.assertRaisesRegex(ConfigError, 'only be used with --program pyANI'):
            Dereplicate(Namespace(program='sourmash', pyani_plus_program='/custom/pyani-plus'))

        with self.assertRaisesRegex(ConfigError, 'cannot be used when importing'):
            Dereplicate(Namespace(program='pyANI', ani_dir='existing-results',
                                  pyani_plus_program='/custom/pyani-plus'))

        with self.assertRaisesRegex(ConfigError, 'cannot be used when importing'):
            Dereplicate(Namespace(program='pyANI', ani_dir='existing-results', ani_backend='legacy'))


if __name__ == '__main__':
    unittest.main()

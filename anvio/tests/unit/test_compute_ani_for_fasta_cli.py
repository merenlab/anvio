"""CLI defaults for the single-FASTA ANI command."""

import contextlib
import io
import sys
import unittest
from unittest import mock

from anvio.cli.compute_ani_for_fasta import get_args, run_program
from anvio.errors import ConfigError


class ComputeANIForFastaCLITestCase(unittest.TestCase):
    def test_help_describes_pyani_plus_cpu_allocation(self):
        output = io.StringIO()
        with mock.patch.object(sys, 'argv', ['anvi-script-compute-ani-for-fasta', '--help']):
            with contextlib.redirect_stdout(output):
                with self.assertRaises(SystemExit) as error:
                    get_args()
        self.assertEqual(error.exception.code, 0)
        self.assertIn('CPU affinity', output.getvalue())
        self.assertIn('taskset', output.getvalue())

    def test_pyani_plus_is_default_and_legacy_remains_selectable(self):
        argv = [
            'anvi-script-compute-ani-for-fasta',
            '-f', 'contigs.fa',
            '-o', 'ani-output',
        ]
        with mock.patch.object(sys, 'argv', argv):
            args = get_args()
        self.assertEqual(args.ani_backend, 'pyani-plus')
        self.assertIsNone(args.pyani_plus_program)

        argv += ['--ani-backend', 'legacy']
        with mock.patch.object(sys, 'argv', argv):
            args = get_args()
        self.assertEqual(args.ani_backend, 'legacy')

        argv += ['--pyani-plus-program', '/opt/ani/bin/pyani-plus']
        with mock.patch.object(sys, 'argv', argv):
            args = get_args()
        self.assertEqual(args.pyani_plus_program, '/opt/ani/bin/pyani-plus')

    def test_legacy_backend_rejects_pyani_plus_program(self):
        args = mock.Mock(method='ANIb', ani_backend='legacy', pyani_plus_program='/custom/pyani-plus')
        with mock.patch('anvio.cli.compute_ani_for_fasta.get_args', return_value=args):
            with self.assertRaisesRegex(ConfigError, '--pyani-plus-program.*pyANI-plus backend'):
                run_program()

    def test_retired_method_reports_migration_guidance_before_running(self):
        args = mock.Mock(method='ANIblastall')
        with mock.patch('anvio.cli.compute_ani_for_fasta.get_args', return_value=args):
            with self.assertRaises(ConfigError) as error:
                run_program()
        self.assertIn('ANIblastall', str(error.exception))
        self.assertIn('retired', str(error.exception))
        self.assertIn('not guaranteed', str(error.exception))
        self.assertIn('match ANIblastall', str(error.exception))


if __name__ == '__main__':
    unittest.main()

"""Tests for anvi'o argparse helpers."""

import contextlib
import io
import os
import unittest

import anvio

from anvio.argparse import ArgumentParser
from anvio.errors import ConfigError


class AnvioThreadsEnvironmentTestCase(unittest.TestCase):
    def setUp(self):
        self.old_anvio_threads = os.environ.pop('ANVIO_THREADS', None)


    def tearDown(self):
        if self.old_anvio_threads is None:
            os.environ.pop('ANVIO_THREADS', None)
        else:
            os.environ['ANVIO_THREADS'] = self.old_anvio_threads


    def get_parser(self, include_num_threads=True):
        parser = ArgumentParser(description="Test parser")

        if include_num_threads:
            parser.add_argument(*anvio.A('num-threads'), **anvio.K('num-threads'))

        return parser


    def test_num_threads_default_without_environment(self):
        args, unknown = self.get_parser().parse_known_args([])

        self.assertEqual(args.num_threads, 1)
        self.assertEqual(unknown, [])


    def test_anvio_threads_sets_default(self):
        os.environ['ANVIO_THREADS'] = '8'

        args, unknown = self.get_parser().parse_known_args([])

        self.assertEqual(args.num_threads, 8)
        self.assertEqual(unknown, [])


    def test_short_num_threads_arg_overrides_environment(self):
        os.environ['ANVIO_THREADS'] = '8'

        args, unknown = self.get_parser().parse_known_args(['-T', '2'])

        self.assertEqual(args.num_threads, 2)
        self.assertEqual(unknown, [])


    def test_long_num_threads_arg_overrides_environment(self):
        os.environ['ANVIO_THREADS'] = '8'

        args, unknown = self.get_parser().parse_known_args(['--num-threads', '3'])

        self.assertEqual(args.num_threads, 3)
        self.assertEqual(unknown, [])


    def test_invalid_anvio_threads_fails_when_num_threads_is_not_provided(self):
        for invalid_value in ('zero', '0', '-1', ''):
            with self.subTest(invalid_value=invalid_value):
                os.environ['ANVIO_THREADS'] = invalid_value

                with self.assertRaises(ConfigError):
                    self.get_parser().parse_known_args([])


    def test_invalid_anvio_threads_is_ignored_when_num_threads_is_provided(self):
        os.environ['ANVIO_THREADS'] = 'zero'

        args, unknown = self.get_parser().parse_known_args(['-T', '2'])

        self.assertEqual(args.num_threads, 2)
        self.assertEqual(unknown, [])


    def test_anvio_threads_is_ignored_without_num_threads_arg(self):
        os.environ['ANVIO_THREADS'] = 'zero'

        args, unknown = self.get_parser(include_num_threads=False).parse_known_args([])

        self.assertFalse(hasattr(args, 'num_threads'))
        self.assertEqual(unknown, [])


    def test_anvio_threads_is_ignored_for_help(self):
        os.environ['ANVIO_THREADS'] = 'zero'

        with contextlib.redirect_stdout(io.StringIO()):
            with self.assertRaises(SystemExit) as context:
                self.get_parser().parse_known_args(['--help'])

        self.assertEqual(context.exception.code, 0)


if __name__ == '__main__':
    unittest.main()

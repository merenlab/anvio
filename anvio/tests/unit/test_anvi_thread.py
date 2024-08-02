# -*- coding: utf-8
"""
    Unit tests for the AnviThread class.
"""

import subprocess
import unittest

import anvio

from anvio.threadingops import AnviThread

__copyright__ = "Copyleft 2015-2020, the Meren Lab (http://merenlab.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "Ryan Moore"
__email__ = "moorer@udel.edu"


def run_script_in_thread(runner, script):
    args = (['python', '-c', script],)
    thread = AnviThread(target=runner, args=args)
    thread.start()
    thread.join()

    return thread


class AnviThreadSavesReturnCodesTestCase(unittest.TestCase):
    def setUp(self):
        self.command_runner = subprocess.call

        self.successful_script = "import sys; sys.exit(0)"
        self.successful_exit_code = 0

        self.failing_script = "import sys; sys.exit(1)"
        self.failing_exit_code = 1

    def test_successful_program(self):
        thread = run_script_in_thread(self.command_runner, self.successful_script)

        self.assertEqual(thread.target_return_value, self.successful_exit_code)

    def test_failing_program(self):
        thread = run_script_in_thread(self.command_runner, self.failing_script)

        self.assertEqual(thread.target_return_value, self.failing_exit_code)

if __name__ == '__main__':
    unittest.main()

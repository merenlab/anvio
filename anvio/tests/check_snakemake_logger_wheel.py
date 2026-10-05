"""Run with a wheel-installed Python, outside the checkout with PYTHONPATH unset."""

import importlib.metadata
import json
import os
from pathlib import Path
import sys
import unittest

import snakemake_logger_plugin_anvio
from snakemake_interface_logger_plugins.registry import LoggerPluginRegistry


class WheelLoggerDiscoveryTestCase(unittest.TestCase):
    def test_wheel_is_automatically_discovered(self):
        self.assertNotIn('PYTHONPATH', os.environ)
        self.assertTrue(Path(snakemake_logger_plugin_anvio.__file__).resolve().is_relative_to(Path(sys.prefix).resolve()))
        direct_url = importlib.metadata.distribution('anvio').read_text('direct_url.json')
        if direct_url:
            self.assertFalse(json.loads(direct_url).get('dir_info', {}).get('editable', False))
        # Do not call the registration helper: this checks Snakemake's package scan.
        registry = LoggerPluginRegistry()
        self.assertTrue(registry.is_installed('anvio'))
        self.assertIs(registry.get_plugin('anvio').log_handler, snakemake_logger_plugin_anvio.LogHandler)


if __name__ == '__main__':
    unittest.main()

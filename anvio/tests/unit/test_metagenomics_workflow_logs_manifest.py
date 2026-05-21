# pylint: disable=line-too-long
"""
    Unit tests for the metagenomics workflow logs manifest.
"""

import os
import shutil
import tempfile
import unittest

import anvio

from anvio.workflows.metagenomics import MetagenomicsWorkflow

__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "OpenAI"
__email__ = "n/a"


class MetagenomicsWorkflowLogsManifestTestCase(unittest.TestCase):
    def setUp(self):
        self.temp_dir = tempfile.mkdtemp()
        self.logs_dir = os.path.join(self.temp_dir, "00_LOGS")
        os.makedirs(self.logs_dir)

        self.workflow = MetagenomicsWorkflow.__new__(MetagenomicsWorkflow)
        self.workflow.dirs_dict = {"LOGS_DIR": self.logs_dir}
        self.workflow.config = {
            "all_against_all": False,
            "gzip_fastqs": {"run": True},
            "import_percent_of_reads_mapped": {"run": True},
            "idba_ud": {"run": False},
            "metaspades": {"run": False},
            "megahit": {"run": True},
            "flye": {"run": True},
            "anvi_cluster_contigs": {"run": False},
        }
        self.workflow.has_sr = True
        self.workflow.has_lr = True
        self.workflow.run_qc = True
        self.workflow.references_for_removal_txt = None
        self.workflow.run_krakenuniq = False
        self.workflow.references_mode = False
        self.workflow.group_names = ["G_SR", "G_LR"]
        self.workflow.group_sizes = {"G_SR": 1, "G_LR": 1}
        self.workflow.assembly_types = {"G_SR": "SR", "G_LR": "LR"}
        self.workflow.collections = None
        self.workflow.run_summary = False
        self.workflow.run_split = False
        self.workflow.readsets = [
            {"id": "S1_SR", "type": "SR", "group": "G", "reads": {}},
            {"id": "S1_LR", "type": "LR", "group": "G", "reads": {}},
        ]
        self.workflow.readset_ids = [readset["id"] for readset in self.workflow.readsets]
        self.workflow._readsets_by_id = {readset["id"]: readset for readset in self.workflow.readsets}


    def tearDown(self):
        shutil.rmtree(self.temp_dir)


    def test_logs_manifest_rows_include_expected_scopes_and_configured_rules(self):
        rows = self.workflow.get_logs_manifest_rows()

        self.assertIn({"rule": "megahit",
                       "scope": "group",
                       "group": "G_SR",
                       "readset": "",
                       "driver": "",
                       "status": "expected",
                       "log_path": os.path.join(self.logs_dir, "G_SR-megahit.log")}, rows)

        self.assertIn({"rule": "flye",
                       "scope": "group",
                       "group": "G_LR",
                       "readset": "",
                       "driver": "",
                       "status": "expected",
                       "log_path": os.path.join(self.logs_dir, "G_LR-flye.log")}, rows)

        self.assertIn({"rule": "bowtie",
                       "scope": "group-readset",
                       "group": "G_SR",
                       "readset": "S1_SR",
                       "driver": "",
                       "status": "expected",
                       "log_path": os.path.join(self.logs_dir, "G_SR-S1_SR-bowtie.log")}, rows)

        self.assertIn({"rule": "minimap2",
                       "scope": "group-readset",
                       "group": "G_LR",
                       "readset": "S1_LR",
                       "driver": "",
                       "status": "expected",
                       "log_path": os.path.join(self.logs_dir, "G_LR-S1_LR-minimap2.log")}, rows)

        self.assertNotIn("krakenuniq", {row["rule"] for row in rows})


    def test_write_logs_manifest_writes_header_and_rows(self):
        self.workflow.write_logs_manifest()

        manifest_path = os.path.join(self.logs_dir, "LOGS-MANIFEST.txt")
        self.assertTrue(os.path.exists(manifest_path))

        with open(manifest_path) as manifest:
            lines = [line.strip() for line in manifest.readlines()]

        self.assertEqual(lines[2], "rule\tscope\tgroup\treadset\tdriver\tstatus\tlog_path")
        self.assertIn("megahit\tgroup\tG_SR\t\t\texpected\t%s" % os.path.join(self.logs_dir, "G_SR-megahit.log"), lines)
        self.assertIn("flye\tgroup\tG_LR\t\t\texpected\t%s" % os.path.join(self.logs_dir, "G_LR-flye.log"), lines)


if __name__ == '__main__':
    unittest.main()

"""Snakemake freshness checks for the mapped-read profile-layer rule."""

import os
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[4]
RULE_SOURCE = REPO_ROOT / "anvio/workflows/read_recruitment/rules/main.smk"
RULE_TEXT = RULE_SOURCE.read_text()
RULE_START = RULE_TEXT.index("rule import_percent_of_reads_mapped:")
RULE_BLOCK = RULE_TEXT[RULE_START:]


@unittest.skipUnless(shutil.which("snakemake"), "Snakemake is not on PATH")
class MappedReadRerunFreshnessTest(unittest.TestCase):
    def setUp(self):
        self.tempdir = tempfile.TemporaryDirectory(prefix="anvio-mapped-rerun-")
        self.root = Path(self.tempdir.name)
        self.profile_dir = self.root / "profile" / "GroupA" / "DAY_15A"
        self.qc_dir = self.root / "qc"
        self.mapping_dir = self.root / "mapping" / "GroupA"
        self.logs_dir = self.root / "logs"
        for directory in (self.profile_dir, self.qc_dir, self.mapping_dir, self.logs_dir):
            directory.mkdir(parents=True)

        self.profile = self.profile_dir / "PROFILE.db"
        self.total_reads = self.qc_dir / "DAY_15A-total_num_reads.txt"
        self.bam = self.mapping_dir / "DAY_15A.bam"
        self.output = self.profile_dir / "layers-additional-data.txt"
        self.profile.touch()
        self.total_reads.write_text("100\n")
        self.bam.touch()
        self.output.write_text("existing layer TSV\n")

        prefix = f'''from pathlib import Path
DIR = {str(self.root)!r}
LOGS = str(Path(DIR) / "logs")
dirs_dict = {{
    "QC_DIR": str(Path(DIR) / "qc"),
    "PROFILE_DIR": str(Path(DIR) / "profile"),
    "MAPPING_DIR": str(Path(DIR) / "mapping"),
}}
ALL_RS_RE = "DAY_15A"
def rule_log(name, wildcards):
    return str(Path(LOGS) / (name + "-" + wildcards + ".log"))

'''
        (self.root / "Snakefile").write_text(prefix + RULE_BLOCK + "\n")
        self.env = os.environ.copy()
        self.env["XDG_CACHE_HOME"] = str(self.root / "cache")

    def tearDown(self):
        self.tempdir.cleanup()

    def dryrun(self):
        result = subprocess.run(
            [
                shutil.which("snakemake"),
                "--snakefile",
                str(self.root / "Snakefile"),
                "--dry-run",
                "--cores",
                "1",
                str(self.output),
            ],
            cwd=self.root,
            env=self.env,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            check=False,
        )
        self.assertEqual(result.returncode, 0, result.stdout)
        return result.stdout

    @staticmethod
    def set_mtime(path, timestamp):
        os.utime(path, (timestamp, timestamp))

    def test_mutable_profile_does_not_hide_read_input_changes(self):
        # Both real data inputs predate the layer TSV. PROFILE.db is newer but
        # mutable, so it must not make an otherwise-current layer import stale.
        self.set_mtime(self.total_reads, 1700000000)
        self.set_mtime(self.bam, 1700000000)
        self.set_mtime(self.output, 1700000100)
        self.set_mtime(self.profile, 1700000200)

        self.assertIn("Nothing to be done", self.dryrun())

        # Mapping output and read totals are real dependencies; changing either
        # one must schedule the layer import even when PROFILE.db is ancient.
        self.set_mtime(self.bam, 1700000300)
        self.assertIn("rule import_percent_of_reads_mapped:", self.dryrun())

        self.set_mtime(self.bam, 1700000000)
        self.set_mtime(self.total_reads, 1700000400)
        self.assertIn("rule import_percent_of_reads_mapped:", self.dryrun())

    def test_rule_keeps_importer_overwrite_scoped_to_layers(self):
        self.assertIn("profiledb=ancient(", RULE_BLOCK)
        self.assertIn('bam=dirs_dict["MAPPING_DIR"]', RULE_BLOCK)
        self.assertIn('total_reads=dirs_dict["QC_DIR"]', RULE_BLOCK)
        self.assertIn('"--target-data-table layers "', RULE_BLOCK)
        self.assertIn('"--just-do-it >> {log} 2>&1"', RULE_BLOCK)


if __name__ == "__main__":
    unittest.main()

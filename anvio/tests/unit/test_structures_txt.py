# pylint: disable=line-too-long
"""Unit tests for the `structures-txt` artifact and the structure-mode branch of
`Pangenome.gen_mcl_input`."""

import os
import shutil
import tempfile
import unittest
from pathlib import Path

import anvio
import anvio.terminal as terminal

from anvio.artifacts.structures_txt import StructuresTxt
from anvio.errors import ConfigError, FilesNPathsError
from anvio.panops import Pangenome


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__


def _make_pdb(dir_, name):
    """Write a placeholder PDB file under dir_ and return its absolute path."""
    p = os.path.join(dir_, name)
    with open(p, 'w') as f:
        f.write('ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C\n')
    return p


def _write_txt(dir_, name, content):
    p = os.path.join(dir_, name)
    with open(p, 'w') as f:
        f.write(content)
    return p


class TestStructuresTxt(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.pdb1 = _make_pdb(self.tmp, 'gc1.pdb')
        self.pdb2 = _make_pdb(self.tmp, 'gc2.pdb')

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def _resolved(self, p):
        return str(Path(p).resolve())

    def test_canonical_header(self):
        path = _write_txt(self.tmp, 's.txt', f"gene_id\tpath\nGC_00000001\t{self.pdb1}\nGC_00000002\t{self.pdb2}\n")
        s = StructuresTxt(path)
        self.assertEqual(s.gene_ids(), ['GC_00000001', 'GC_00000002'])
        self.assertEqual(s.get_path('GC_00000001'), self._resolved(self.pdb1))
        self.assertEqual(len(s), 2)
        self.assertIn('GC_00000002', s)

    def test_legacy_gene_callers_id_header_accepted(self):
        path = _write_txt(self.tmp, 's.txt', f"gene_callers_id\tpath\n1\t{self.pdb1}\n2\t{self.pdb2}\n")
        s = StructuresTxt(path)
        self.assertEqual(s.gene_ids(), ['1', '2'])

    def test_unknown_id_header_rejected(self):
        path = _write_txt(self.tmp, 's.txt', f"protein_id\tpath\nX\t{self.pdb1}\n")
        with self.assertRaises(ConfigError):
            StructuresTxt(path)

    def test_wrong_path_header_rejected(self):
        path = _write_txt(self.tmp, 's.txt', f"gene_id\tfile\nX\t{self.pdb1}\n")
        with self.assertRaises(ConfigError):
            StructuresTxt(path)

    def test_wrong_column_count_rejected(self):
        path = _write_txt(self.tmp, 's.txt', f"gene_id\tpath\textra\nX\t{self.pdb1}\tY\n")
        with self.assertRaises(ConfigError):
            StructuresTxt(path)

    def test_duplicate_gene_id_rejected(self):
        path = _write_txt(self.tmp, 's.txt', f"gene_id\tpath\nGC_1\t{self.pdb1}\nGC_1\t{self.pdb2}\n")
        with self.assertRaises(ConfigError) as cm:
            StructuresTxt(path)
        self.assertIn('GC_1', str(cm.exception))

    def test_duplicate_message_uses_legacy_label_when_legacy_header_used(self):
        path = _write_txt(self.tmp, 's.txt', f"gene_callers_id\tpath\n7\t{self.pdb1}\n7\t{self.pdb2}\n")
        with self.assertRaises(ConfigError) as cm:
            StructuresTxt(path)
        self.assertIn('gene_callers_id', str(cm.exception))

    def test_missing_file_rejected(self):
        path = _write_txt(self.tmp, 's.txt', "gene_id\tpath\nGC_1\t/nonexistent/foo.pdb\n")
        with self.assertRaises(FilesNPathsError):
            StructuresTxt(path)

    def test_disallowed_extension_rejected(self):
        bad = _write_txt(self.tmp, 'x.txt', 'not a structure')
        path = _write_txt(self.tmp, 's.txt', f"gene_id\tpath\nGC_1\t{bad}\n")
        with self.assertRaises(ConfigError):
            StructuresTxt(path)

    def test_relative_paths_resolved_against_txt_dir(self):
        # use 'gc1.pdb' (a relative path) in the TSV; should resolve against self.tmp
        path = _write_txt(self.tmp, 's.txt', "gene_id\tpath\nGC_1\tgc1.pdb\n")
        s = StructuresTxt(path)
        self.assertEqual(s.get_path('GC_1'), self._resolved(self.pdb1))

    def test_multi_suffix_extensions_accepted(self):
        gz = os.path.join(self.tmp, 'foo.cif.gz')
        with open(gz, 'w') as f:
            f.write('placeholder')
        path = _write_txt(self.tmp, 's.txt', f"gene_id\tpath\nGC_1\t{gz}\n")
        self.assertEqual(StructuresTxt(path).get_path('GC_1'), self._resolved(gz))

    def test_empty_body_rejected(self):
        path = _write_txt(self.tmp, 's.txt', "gene_id\tpath\n")
        with self.assertRaises(ConfigError):
            StructuresTxt(path)


class TestStructureModeGenMclInput(unittest.TestCase):
    """Exercise the structure-mode branch of Pangenome.gen_mcl_input against a synthetic
    14-column foldseek output, without instantiating a real Pangenome (the relevant code
    only touches a handful of attributes).

    Contract: structure-mode gen_mcl_input reads exactly these attributes off `self`:
        - self.run, self.progress  (anvi'o terminal helpers)
        - self.min_tm_score        (filter threshold)
        - self.get_output_file_path (factory that returns a writable path under output_dir)

    If a future change to gen_mcl_input starts reading a new attribute, these tests will
    fail with an AttributeError rather than a meaningful diagnostic. Update setUp below to
    bind whatever new attribute is being added.
    """

    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        p = Pangenome.__new__(Pangenome)
        p.run = terminal.Run(verbose=False)
        p.progress = terminal.Progress(verbose=False)
        p.min_tm_score = 0.5
        p.output_dir = self.tmp
        p.get_output_file_path = lambda name: os.path.join(self.tmp, name)
        self.p = p

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def _write_foldseek_output(self, rows):
        path = os.path.join(self.tmp, 'fs.tsv')
        with open(path, 'w') as f:
            for r in rows:
                f.write('\t'.join(str(x) for x in r) + '\n')
        return path

    @staticmethod
    def _row(query, target, qtm, ttm):
        # 14-column foldseek easy-search output: query, target, fident, alnlen, mismatch,
        # gapopen, qstart, qend, tstart, tend, evalue, bits, qtmscore, ttmscore
        return (query, target, 1.0, 100, 0, 0, 1, 100, 1, 100, 0.0, 200.0, qtm, ttm)

    def _read_mcl(self, path):
        edges = []
        for line in open(path):
            q, t, w = line.rstrip().split('\t')
            edges.append((q, t, float(w)))
        return edges

    def test_tm_min_used_as_weight(self):
        # qtm=0.7, ttm=0.6 -> min=0.6 kept and used as weight
        out = self._write_foldseek_output([
            self._row('GC_001', 'GC_001', 1.0, 1.0),
            self._row('GC_002', 'GC_002', 1.0, 1.0),
            self._row('GC_001', 'GC_002', 0.7, 0.6),
        ])
        mcl_path = self.p.gen_mcl_input(out, mode='structure')
        edges = self._read_mcl(mcl_path)
        self.assertIn(('GC_001', 'GC_002', 0.6), edges)

    def test_below_threshold_dropped(self):
        # qtm=0.4, ttm=0.8 -> min=0.4 < 0.5 default, dropped
        out = self._write_foldseek_output([
            self._row('GC_001', 'GC_001', 1.0, 1.0),
            self._row('GC_002', 'GC_002', 1.0, 1.0),
            self._row('GC_001', 'GC_002', 0.4, 0.8),
        ])
        mcl_path = self.p.gen_mcl_input(out, mode='structure')
        edges = self._read_mcl(mcl_path)
        # only the two self-hits survive
        self.assertEqual(len(edges), 2)
        for q, t, _ in edges:
            self.assertEqual(q, t)

    def test_missing_self_hit_gets_self_loop(self):
        # GC_002 has no self-hit row, so the heuristic should inject one with weight 1.0
        out = self._write_foldseek_output([
            self._row('GC_001', 'GC_001', 1.0, 1.0),
            self._row('GC_001', 'GC_002', 0.8, 0.9),
        ])
        mcl_path = self.p.gen_mcl_input(out, mode='structure')
        edges = self._read_mcl(mcl_path)
        self.assertIn(('GC_002', 'GC_002', 1.0), edges)

    def test_expected_ids_catches_fully_absent_structure(self):
        # GC_003 doesn't appear anywhere in the foldseek output. Without expected_ids it would
        # be invisible. With expected_ids, the heuristic should add a self-loop for it.
        out = self._write_foldseek_output([
            self._row('GC_001', 'GC_001', 1.0, 1.0),
            self._row('GC_002', 'GC_002', 1.0, 1.0),
        ])
        mcl_path = self.p.gen_mcl_input(out, mode='structure',
                                        expected_ids={'GC_001', 'GC_002', 'GC_003'})
        edges = self._read_mcl(mcl_path)
        self.assertIn(('GC_003', 'GC_003', 1.0), edges)

    def test_threshold_is_inclusive(self):
        # exactly equal to threshold is kept
        out = self._write_foldseek_output([
            self._row('GC_001', 'GC_001', 1.0, 1.0),
            self._row('GC_002', 'GC_002', 1.0, 1.0),
            self._row('GC_001', 'GC_002', 0.5, 0.5),
        ])
        mcl_path = self.p.gen_mcl_input(out, mode='structure')
        edges = self._read_mcl(mcl_path)
        self.assertIn(('GC_001', 'GC_002', 0.5), edges)

    def test_malformed_row_raises_config_error(self):
        bad = os.path.join(self.tmp, 'bad.tsv')
        with open(bad, 'w') as f:
            f.write('GC_001\tGC_002\tnot_a_number\n')
        with self.assertRaises(ConfigError):
            self.p.gen_mcl_input(bad, mode='structure')

    def test_unknown_mode_raises_config_error(self):
        out = self._write_foldseek_output([self._row('a', 'a', 1, 1)])
        with self.assertRaises(ConfigError):
            self.p.gen_mcl_input(out, mode='banana')

    def test_query_ids_with_pdb_extension_get_stripped(self):
        # Some foldseek builds emit 'GC_001.pdb' in the query/target columns even though
        # we hand them gene_id-named symlinks. The defensive _strip_foldseek_id_ext call in
        # _parse_search_row should normalize them.
        out = self._write_foldseek_output([
            self._row('GC_001.pdb', 'GC_001.pdb', 1.0, 1.0),
            self._row('GC_002.pdb', 'GC_002.pdb', 1.0, 1.0),
            self._row('GC_001.pdb', 'GC_002.pdb', 0.8, 0.9),
        ])
        mcl_path = self.p.gen_mcl_input(out, mode='structure',
                                        expected_ids={'GC_001', 'GC_002'})
        edges = self._read_mcl(mcl_path)
        # IDs in the MCL output should be the stripped gene_ids
        for q, t, _ in edges:
            self.assertFalse(q.endswith('.pdb'))
            self.assertFalse(t.endswith('.pdb'))
        # And expected_ids comparison should succeed - no spurious "missing structure" warning
        # (we'd see it if .pdb-suffixed IDs failed to match expected_ids 'GC_001'/'GC_002')
        self.assertEqual(len(edges), 3)  # 2 self-hits + 1 pair


class TestStripFoldseekIdExt(unittest.TestCase):
    def test_strips_pdb(self):
        self.assertEqual(Pangenome._strip_foldseek_id_ext('GC_001.pdb'), 'GC_001')

    def test_strips_case_insensitive(self):
        self.assertEqual(Pangenome._strip_foldseek_id_ext('GC_001.PDB'), 'GC_001')

    def test_strips_multi_suffix(self):
        # Longest matching suffix wins: '.cif.gz' before '.gz'
        self.assertEqual(Pangenome._strip_foldseek_id_ext('foo.cif.gz'), 'foo')

    def test_passthrough_when_no_suffix(self):
        self.assertEqual(Pangenome._strip_foldseek_id_ext('GC_001'), 'GC_001')

    def test_passthrough_when_unknown_suffix(self):
        self.assertEqual(Pangenome._strip_foldseek_id_ext('GC_001.xyz'), 'GC_001.xyz')


if __name__ == '__main__':
    unittest.main()

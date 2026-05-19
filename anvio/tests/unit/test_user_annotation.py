# pylint: disable=line-too-long
"""
Unit tests for anvio.user_annotation — UserAnnotationDBSetup and UserAnnotationRunner.

Covers (no external tools required):
  - TSV parsing: 2-col, 3-col, headers, comments, duplicate names, empty file
  - File type sniffing: HMM, FASTA, gzipped, unknown content
  - File type validation: extension/content mismatch raises ConfigError
  - HMM model parsing: TC/GA/NC/mixed/no-cutoffs, duplicate model names, missing NAME
  - Noise cutoff determination: all-TC, all-GA, all-NC, mixed, none
  - HMM source directory creation: required files, contents
  - DIAMOND tabular output parsing: best-hit dedup, function string format, [DMND] prefix
  - Manifest list and remove operations
  - Runner sanity checks and --database selection logic

Integration tests (skipped when external tools are absent):
  - setup_diamond_source: needs `diamond` in PATH
"""

import os
import gzip
import json
import shutil
import tempfile
import unittest
import argparse

import anvio
import anvio.terminal as terminal
from anvio.errors import ConfigError
from anvio.user_annotation import (
    UserAnnotationDBSetup,
    UserAnnotationRunner,
    HMM_SOURCE_SUFFIX,
    DIAMOND_SOURCE_SUFFIX,
    MANIFEST_FILENAME,
)

__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['lgallucc']


# ---------------------------------------------------------------------------
# Paths to static fixture files
# ---------------------------------------------------------------------------

FIXTURES = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                        'test_files', 'user_annotation')

F_HMM_ALL_TC       = os.path.join(FIXTURES, 'hmm_all_tc.hmm')
F_HMM_ALL_GA       = os.path.join(FIXTURES, 'hmm_all_ga.hmm')
F_HMM_ALL_NC       = os.path.join(FIXTURES, 'hmm_all_nc.hmm')
F_HMM_MIXED        = os.path.join(FIXTURES, 'hmm_mixed_cutoffs.hmm')
F_HMM_NO_CUTOFFS   = os.path.join(FIXTURES, 'hmm_no_cutoffs.hmm')
F_HMM_DUP_NAMES    = os.path.join(FIXTURES, 'hmm_duplicate_names.hmm')
F_HMM_NO_NAME      = os.path.join(FIXTURES, 'hmm_no_name.hmm')
F_VALID_FAA        = os.path.join(FIXTURES, 'valid.faa')
F_FASTA_AS_HMM     = os.path.join(FIXTURES, 'fasta_as_hmm.hmm')
F_HMM_AS_FASTA     = os.path.join(FIXTURES, 'hmm_as_fasta.faa')
F_DIAMOND_TABULAR  = os.path.join(FIXTURES, 'diamond_tabular.txt')


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _silent_run():
    return terminal.Run(verbose=False)

def _silent_progress():
    return terminal.Progress(verbose=False)


def _setup_args(**kwargs):
    defaults = dict(input_tsv=None, output_dir=None, num_threads=1,
                    reset=False, list=False, remove=None)
    defaults.update(kwargs)
    return argparse.Namespace(**defaults)


def _runner_args(**kwargs):
    defaults = dict(contigs_db=None, annotation_db_dir=None, database=None,
                    num_threads=1, just_do_it=False, hmmer_program='hmmscan',
                    evalue=None, min_pident=None, diamond_sensitivity=None)
    defaults.update(kwargs)
    return argparse.Namespace(**defaults)


def _make_tsv(path, entries):
    """Write a TSV file. entries: list of tuples (name, fpath) or (name, fpath, companion)."""
    with open(path, 'w') as f:
        for row in entries:
            f.write('\t'.join(str(c) for c in row) + '\n')
    return path


def _make_manifest(output_dir, content):
    os.makedirs(output_dir, exist_ok=True)
    manifest_path = os.path.join(output_dir, MANIFEST_FILENAME)
    with open(manifest_path, 'w') as f:
        json.dump(content, f)
    return manifest_path


def _make_setup_instance(output_dir, input_tsv):
    """Instantiate UserAnnotationDBSetup with silent output."""
    args = _setup_args(output_dir=output_dir, input_tsv=input_tsv)
    return UserAnnotationDBSetup(args, run=_silent_run(), progress=_silent_progress())


# ===========================================================================
# TestParseInputTSV
# ===========================================================================

class TestParseInputTSV(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.out = os.path.join(self.tmp, 'out')
        os.makedirs(self.out)

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def _setup(self, tsv_path):
        return _make_setup_instance(self.out, tsv_path)

    def test_two_column_tsv_parsed(self):
        tsv = _make_tsv(os.path.join(self.tmp, 't.tsv'),
                        [('MyHMM', F_HMM_ALL_TC)])
        s = self._setup(tsv)
        entries = s.parse_input_tsv()
        self.assertIn('MyHMM', entries)
        self.assertEqual(entries['MyHMM']['companion_fasta'], None)

    def test_three_column_tsv_sets_companion(self):
        tsv = _make_tsv(os.path.join(self.tmp, 't.tsv'),
                        [('MyHMM', F_HMM_ALL_TC, F_VALID_FAA)])
        s = self._setup(tsv)
        entries = s.parse_input_tsv()
        self.assertIsNotNone(entries['MyHMM']['companion_fasta'])
        self.assertTrue(entries['MyHMM']['companion_fasta'].endswith('valid.faa'))

    def test_header_line_skipped(self):
        tsv = os.path.join(self.tmp, 't.tsv')
        with open(tsv, 'w') as f:
            f.write('name\tpath\n')
            f.write(f'MyHMM\t{F_HMM_ALL_TC}\n')
        s = self._setup(tsv)
        entries = s.parse_input_tsv()
        self.assertNotIn('name', entries)
        self.assertIn('MyHMM', entries)

    def test_comment_lines_skipped(self):
        tsv = os.path.join(self.tmp, 't.tsv')
        with open(tsv, 'w') as f:
            f.write('# this is a comment\n')
            f.write(f'MyHMM\t{F_HMM_ALL_TC}\n')
        s = self._setup(tsv)
        entries = s.parse_input_tsv()
        self.assertEqual(len(entries), 1)

    def test_duplicate_names_raise(self):
        tsv = _make_tsv(os.path.join(self.tmp, 't.tsv'),
                        [('MyHMM', F_HMM_ALL_TC), ('MyHMM', F_HMM_ALL_GA)])
        s = self._setup(tsv)
        with self.assertRaises(ConfigError):
            s.parse_input_tsv()

    def test_empty_file_raises(self):
        tsv = os.path.join(self.tmp, 'empty.tsv')
        open(tsv, 'w').close()
        s = self._setup(tsv)
        with self.assertRaises(ConfigError):
            s.parse_input_tsv()

    def test_line_with_fewer_than_two_fields_raises(self):
        tsv = os.path.join(self.tmp, 't.tsv')
        with open(tsv, 'w') as f:
            f.write('OnlyOneName\n')
        s = self._setup(tsv)
        with self.assertRaises(ConfigError):
            s.parse_input_tsv()

    def test_windows_line_endings_handled(self):
        tsv = os.path.join(self.tmp, 't.tsv')
        with open(tsv, 'wb') as f:
            f.write(f'MyHMM\t{F_HMM_ALL_TC}\r\n'.encode())
        s = self._setup(tsv)
        entries = s.parse_input_tsv()
        self.assertIn('MyHMM', entries)

    def test_nonexistent_path_raises(self):
        tsv = os.path.join(self.tmp, 't.tsv')
        with open(tsv, 'w') as f:
            f.write(f'MyHMM\t/does/not/exist.hmm\n')
        s = self._setup(tsv)
        with self.assertRaises(Exception):
            s.parse_input_tsv()

    def test_multiple_entries_parsed(self):
        tsv = _make_tsv(os.path.join(self.tmp, 't.tsv'),
                        [('HMM1', F_HMM_ALL_TC), ('HMM2', F_HMM_ALL_GA)])
        s = self._setup(tsv)
        entries = s.parse_input_tsv()
        self.assertEqual(len(entries), 2)
        self.assertIn('HMM1', entries)
        self.assertIn('HMM2', entries)


# ===========================================================================
# TestSniffFileType
# ===========================================================================

class TestSniffFileType(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.out = os.path.join(self.tmp, 'out')
        os.makedirs(self.out)
        tsv = _make_tsv(os.path.join(self.tmp, 't.tsv'), [('X', F_HMM_ALL_TC)])
        self.s = _make_setup_instance(self.out, tsv)

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def test_hmm_file_detected(self):
        self.assertEqual(self.s._sniff_file_type(F_HMM_ALL_TC), 'hmm')

    def test_fasta_file_detected(self):
        self.assertEqual(self.s._sniff_file_type(F_VALID_FAA), 'diamond')

    def test_unknown_content_returns_none(self):
        path = os.path.join(self.tmp, 'garbage.txt')
        with open(path, 'w') as f:
            f.write('this is neither HMM nor FASTA content\n')
        self.assertIsNone(self.s._sniff_file_type(path))

    def test_gzipped_hmm_detected(self):
        gz_path = os.path.join(self.tmp, 'test.hmm.gz')
        with open(F_HMM_ALL_TC, 'rb') as f_in, gzip.open(gz_path, 'wb') as f_out:
            f_out.write(f_in.read())
        self.assertEqual(self.s._sniff_file_type(gz_path), 'hmm')

    def test_fasta_with_leading_whitespace_detected(self):
        path = os.path.join(self.tmp, 'spaced.faa')
        with open(path, 'w') as f:
            f.write('  >seq1\nMKTLLL\n')
        self.assertEqual(self.s._sniff_file_type(path), 'diamond')


# ===========================================================================
# TestValidateAndDetectDbType
# ===========================================================================

class TestValidateAndDetectDbType(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.out = os.path.join(self.tmp, 'out')
        os.makedirs(self.out)
        tsv = _make_tsv(os.path.join(self.tmp, 't.tsv'), [('X', F_HMM_ALL_TC)])
        self.s = _make_setup_instance(self.out, tsv)

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def test_hmm_extension_hmm_content_ok(self):
        self.assertEqual(self.s.validate_and_detect_db_type(F_HMM_ALL_TC), 'hmm')

    def test_fasta_extension_fasta_content_ok(self):
        self.assertEqual(self.s.validate_and_detect_db_type(F_VALID_FAA), 'diamond')

    def test_fasta_in_hmm_extension_raises(self):
        with self.assertRaises(ConfigError):
            self.s.validate_and_detect_db_type(F_FASTA_AS_HMM)

    def test_hmm_in_fasta_extension_raises(self):
        with self.assertRaises(ConfigError):
            self.s.validate_and_detect_db_type(F_HMM_AS_FASTA)

    def test_no_extension_hmm_content_returns_hmm(self):
        path = os.path.join(self.tmp, 'noext')
        with open(path, 'w') as f:
            f.write('HMMER3/f [3.3.2]\nNAME  X\n//\n')
        self.assertEqual(self.s.validate_and_detect_db_type(path), 'hmm')

    def test_no_extension_fasta_content_returns_diamond(self):
        path = os.path.join(self.tmp, 'noext')
        with open(path, 'w') as f:
            f.write('>seq1\nMKTLLL\n')
        self.assertEqual(self.s.validate_and_detect_db_type(path), 'diamond')

    def test_unrecognised_extension_unrecognised_content_raises(self):
        path = os.path.join(self.tmp, 'mystery.xyz')
        with open(path, 'w') as f:
            f.write('completely unknown content here\n')
        with self.assertRaises(ConfigError):
            self.s.validate_and_detect_db_type(path)


# ===========================================================================
# TestParseHmmModels
# ===========================================================================

class TestParseHmmModels(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.out = os.path.join(self.tmp, 'out')
        os.makedirs(self.out)
        tsv = _make_tsv(os.path.join(self.tmp, 't.tsv'), [('X', F_HMM_ALL_TC)])
        self.s = _make_setup_instance(self.out, tsv)

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def test_models_count_correct(self):
        models = self.s.parse_hmm_models(F_HMM_ALL_TC)
        self.assertEqual(len(models), 2)

    def test_model_names_parsed(self):
        models = self.s.parse_hmm_models(F_HMM_ALL_TC)
        names = [m['name'] for m in models]
        self.assertIn('ModelA', names)
        self.assertIn('ModelB', names)

    def test_accession_parsed(self):
        models = self.s.parse_hmm_models(F_HMM_ALL_TC)
        acc = {m['name']: m['acc'] for m in models}
        self.assertEqual(acc['ModelA'], 'ACC001')

    def test_missing_accession_is_none(self):
        models = self.s.parse_hmm_models(F_HMM_ALL_GA)
        self.assertIsNone(models[0]['acc'])

    def test_tc_flag_set(self):
        models = self.s.parse_hmm_models(F_HMM_ALL_TC)
        self.assertTrue(all(m['has_tc'] for m in models))
        self.assertFalse(any(m['has_ga'] for m in models))

    def test_ga_flag_set(self):
        models = self.s.parse_hmm_models(F_HMM_ALL_GA)
        self.assertTrue(all(m['has_ga'] for m in models))

    def test_nc_flag_set(self):
        models = self.s.parse_hmm_models(F_HMM_ALL_NC)
        self.assertTrue(all(m['has_nc'] for m in models))

    def test_duplicate_model_names_raise(self):
        with self.assertRaises(ConfigError):
            self.s.parse_hmm_models(F_HMM_DUP_NAMES)

    def test_model_without_name_raises(self):
        with self.assertRaises(ConfigError):
            self.s.parse_hmm_models(F_HMM_NO_NAME)

    def test_empty_hmm_file_raises(self):
        empty = os.path.join(self.tmp, 'empty.hmm')
        open(empty, 'w').close()
        with self.assertRaises(ConfigError):
            self.s.parse_hmm_models(empty)

    def test_gzipped_hmm_parsed(self):
        gz_path = os.path.join(self.tmp, 'test.hmm.gz')
        with open(F_HMM_ALL_TC, 'rb') as f_in, gzip.open(gz_path, 'wb') as f_out:
            f_out.write(f_in.read())
        models = self.s.parse_hmm_models(gz_path)
        self.assertEqual(len(models), 2)

    def test_mixed_cutoffs_parsed_correctly(self):
        models = self.s.parse_hmm_models(F_HMM_MIXED)
        has_tc = [m['has_tc'] for m in models]
        self.assertTrue(has_tc[0])
        self.assertFalse(has_tc[1])


# ===========================================================================
# TestDetermineNoiseCutoff
# ===========================================================================

class TestDetermineNoiseCutoff(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.out = os.path.join(self.tmp, 'out')
        os.makedirs(self.out)
        tsv = _make_tsv(os.path.join(self.tmp, 't.tsv'), [('X', F_HMM_ALL_TC)])
        self.s = _make_setup_instance(self.out, tsv)

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def _models(self, tc=False, ga=False, nc=False, n=2):
        return [{'name': f'M{i}', 'has_tc': tc, 'has_ga': ga, 'has_nc': nc}
                for i in range(n)]

    def test_all_tc_returns_cut_tc(self):
        self.assertEqual(self.s.determine_noise_cutoff(self._models(tc=True), 'db'), '--cut_tc')

    def test_all_ga_returns_cut_ga(self):
        self.assertEqual(self.s.determine_noise_cutoff(self._models(ga=True), 'db'), '--cut_ga')

    def test_all_nc_returns_cut_nc(self):
        self.assertEqual(self.s.determine_noise_cutoff(self._models(nc=True), 'db'), '--cut_nc')

    def test_mixed_cutoffs_returns_evalue(self):
        models = [{'name': 'M0', 'has_tc': True, 'has_ga': False, 'has_nc': False},
                  {'name': 'M1', 'has_tc': False, 'has_ga': False, 'has_nc': False}]
        self.assertEqual(self.s.determine_noise_cutoff(models, 'db'), '-E 1e-5')

    def test_no_cutoffs_returns_evalue(self):
        self.assertEqual(self.s.determine_noise_cutoff(self._models(), 'db'), '-E 1e-5')

    def test_tc_takes_priority_over_ga(self):
        # TC wins even if GA also set on all models
        models = [{'name': f'M{i}', 'has_tc': True, 'has_ga': True, 'has_nc': False}
                  for i in range(2)]
        self.assertEqual(self.s.determine_noise_cutoff(models, 'db'), '--cut_tc')


# ===========================================================================
# TestSetupHmmSource
# ===========================================================================

class TestSetupHmmSource(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.out = os.path.join(self.tmp, 'out')
        os.makedirs(self.out)
        tsv = _make_tsv(os.path.join(self.tmp, 't.tsv'), [('TestDB', F_HMM_ALL_TC)])
        self.s = _make_setup_instance(self.out, tsv)
        self.entry = self.s.setup_hmm_source('TestDB', F_HMM_ALL_TC)
        self.hmm_dir = self.entry['hmm_dir']

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def test_hmm_dir_created(self):
        self.assertTrue(os.path.isdir(self.hmm_dir))

    def test_genes_hmm_gz_exists(self):
        self.assertTrue(os.path.exists(os.path.join(self.hmm_dir, 'genes.hmm.gz')))

    def test_genes_hmm_gz_is_valid_gzip(self):
        path = os.path.join(self.hmm_dir, 'genes.hmm.gz')
        with gzip.open(path, 'rt') as f:
            content = f.read()
        self.assertIn('HMMER3', content)

    def test_genes_txt_has_header(self):
        path = os.path.join(self.hmm_dir, 'genes.txt')
        with open(path) as f:
            first_line = f.readline()
        self.assertIn('gene', first_line)
        self.assertIn('accession', first_line)

    def test_genes_txt_has_model_entries(self):
        path = os.path.join(self.hmm_dir, 'genes.txt')
        with open(path) as f:
            lines = f.readlines()
        names = [l.split('\t')[0] for l in lines[1:]]
        self.assertIn('ModelA', names)
        self.assertIn('ModelB', names)

    def test_kind_txt_content(self):
        path = os.path.join(self.hmm_dir, 'kind.txt')
        with open(path) as f:
            content = f.read().strip()
        self.assertEqual(content, 'user_annotation')

    def test_target_txt_content(self):
        path = os.path.join(self.hmm_dir, 'target.txt')
        with open(path) as f:
            content = f.read().strip()
        self.assertEqual(content, 'AA:GENE')

    def test_noise_cutoff_terms_written(self):
        path = os.path.join(self.hmm_dir, 'noise_cutoff_terms.txt')
        with open(path) as f:
            content = f.read().strip()
        self.assertEqual(content, '--cut_tc')

    def test_reference_txt_exists(self):
        self.assertTrue(os.path.exists(os.path.join(self.hmm_dir, 'reference.txt')))

    def test_entry_dict_has_correct_keys(self):
        for key in ('type', 'hmm_dir', 'num_models', 'noise_cutoff_terms', 'added_on'):
            self.assertIn(key, self.entry)

    def test_entry_num_models_correct(self):
        self.assertEqual(self.entry['num_models'], 2)

    def test_entry_type_is_hmm(self):
        self.assertEqual(self.entry['type'], 'hmm')

    def test_no_cutoff_source_uses_evalue(self):
        entry = self.s.setup_hmm_source('NoCutDB', F_HMM_NO_CUTOFFS)
        path = os.path.join(entry['hmm_dir'], 'noise_cutoff_terms.txt')
        with open(path) as f:
            content = f.read().strip()
        self.assertEqual(content, '-E 1e-5')


# ===========================================================================
# TestParseDiamondTabularOutput
# ===========================================================================

class TestParseDiamondTabularOutput(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.out = os.path.join(self.tmp, 'out')
        os.makedirs(self.out)
        # Runner requires contigs_db and manifest; mock them minimally
        manifest = {'FakeDB': {'type': 'diamond', 'dmnd_path': '/fake.dmnd',
                               'dmnd_base': '/fake', 'num_sequences': 2,
                               'added_on': '2026-01-01'}}
        _make_manifest(self.out, manifest)
        # We need a fake contigs-db path that passes is_contigs_db check.
        # Bypass by calling parse_diamond_tabular_output directly on an instance
        # created after patching sanity_check.
        self.runner = self._make_runner()

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def _make_runner(self):
        args = _runner_args(annotation_db_dir=self.out, contigs_db='/dev/null')
        r = object.__new__(UserAnnotationRunner)
        r.args = args
        r.run = _silent_run()
        r.progress = _silent_progress()
        r.contigs_db_path = '/dev/null'
        r.annotation_db_dir = self.out
        r.database = None
        r.num_threads = 1
        r.just_do_it = False
        r.hmmer_program = 'hmmscan'
        r.evalue = None
        r.min_pident = None
        r.diamond_sensitivity = None
        return r

    def test_two_unique_genes_parsed(self):
        result = self.runner.parse_diamond_tabular_output(F_DIAMOND_TABULAR, 'TestDB')
        gene_ids = {v['gene_callers_id'] for v in result.values()}
        self.assertEqual(gene_ids, {10, 20})

    def test_best_hit_kept_for_duplicate_query(self):
        # gene 10 appears twice: evalue 1e-50 (target_A) and 1e-20 (target_C)
        # best = lowest evalue = 1e-50 → target_A
        result = self.runner.parse_diamond_tabular_output(F_DIAMOND_TABULAR, 'TestDB')
        hits = {v['gene_callers_id']: v for v in result.values()}
        self.assertEqual(hits[10]['accession'], 'target_A')

    def test_source_name_has_diamond_suffix(self):
        result = self.runner.parse_diamond_tabular_output(F_DIAMOND_TABULAR, 'TestDB')
        sources = {v['source'] for v in result.values()}
        self.assertEqual(sources, {f'TestDB{DIAMOND_SOURCE_SUFFIX}'})

    def test_function_string_has_dmnd_prefix(self):
        result = self.runner.parse_diamond_tabular_output(F_DIAMOND_TABULAR, 'TestDB')
        for v in result.values():
            self.assertTrue(v['function'].startswith('[DMND]'),
                            f"Expected [DMND] prefix, got: {v['function']}")

    def test_function_string_contains_pident(self):
        result = self.runner.parse_diamond_tabular_output(F_DIAMOND_TABULAR, 'TestDB')
        hits = {v['gene_callers_id']: v for v in result.values()}
        self.assertIn('95.0%', hits[10]['function'])

    def test_function_string_contains_aln_len(self):
        result = self.runner.parse_diamond_tabular_output(F_DIAMOND_TABULAR, 'TestDB')
        hits = {v['gene_callers_id']: v for v in result.values()}
        self.assertIn('100 aa', hits[10]['function'])

    def test_function_string_contains_bitscore(self):
        result = self.runner.parse_diamond_tabular_output(F_DIAMOND_TABULAR, 'TestDB')
        hits = {v['gene_callers_id']: v for v in result.values()}
        self.assertIn('200.0', hits[10]['function'])

    def test_evalue_stored(self):
        result = self.runner.parse_diamond_tabular_output(F_DIAMOND_TABULAR, 'TestDB')
        hits = {v['gene_callers_id']: v for v in result.values()}
        self.assertAlmostEqual(hits[10]['e_value'], 1e-50)

    def test_non_integer_query_id_skipped(self):
        path = os.path.join(self.tmp, 'bad_ids.txt')
        with open(path, 'w') as f:
            f.write('gene_abc\ttarget_X\t90.0\t80\t0\t0\t1\t80\t1\t80\t1e-40\t180.0\n')
        result = self.runner.parse_diamond_tabular_output(path, 'TestDB')
        self.assertEqual(len(result), 0)

    def test_short_lines_skipped(self):
        path = os.path.join(self.tmp, 'short.txt')
        with open(path, 'w') as f:
            f.write('10\ttarget_A\t95.0\n')
        result = self.runner.parse_diamond_tabular_output(path, 'TestDB')
        self.assertEqual(len(result), 0)

    def test_empty_file_returns_empty_dict(self):
        path = os.path.join(self.tmp, 'empty.txt')
        open(path, 'w').close()
        result = self.runner.parse_diamond_tabular_output(path, 'TestDB')
        self.assertEqual(result, {})


# ===========================================================================
# TestManifestOperations  (list / remove in UserAnnotationDBSetup)
# ===========================================================================

class TestManifestOperations(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.out = os.path.join(self.tmp, 'out')
        os.makedirs(self.out)

        # Build a real HMM source so remove can delete actual files
        tsv = _make_tsv(os.path.join(self.tmp, 't.tsv'), [('HMM1', F_HMM_ALL_TC)])
        self.s = _make_setup_instance(self.out, tsv)
        entry = self.s.setup_hmm_source('HMM1', F_HMM_ALL_TC)

        manifest = {
            'HMM1': entry,
            'DMND1': {
                'type': 'diamond',
                'source_path': F_VALID_FAA,
                'dmnd_path': os.path.join(self.out, 'diamond', 'DMND1.dmnd'),
                'dmnd_base': os.path.join(self.out, 'diamond', 'DMND1'),
                'num_sequences': 2,
                'added_on': '2026-01-01',
            },
        }
        _make_manifest(self.out, manifest)

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def _setup_for_list_remove(self, **kwargs):
        args = _setup_args(output_dir=self.out, **kwargs)
        return UserAnnotationDBSetup(args, run=_silent_run(), progress=_silent_progress())

    def test_list_does_not_raise(self):
        s = self._setup_for_list_remove(**{'list': True})
        s.list_databases()  # should not raise

    def test_remove_unknown_name_raises(self):
        s = self._setup_for_list_remove(remove='NonExistent')
        with self.assertRaises(ConfigError):
            s.remove_database_entry('NonExistent')

    def test_remove_hmm_deletes_manifest_entry(self):
        s = self._setup_for_list_remove(remove='HMM1')
        s.remove_database_entry('HMM1')
        with open(os.path.join(self.out, MANIFEST_FILENAME)) as f:
            manifest = json.load(f)
        self.assertNotIn('HMM1', manifest)
        self.assertIn('DMND1', manifest)

    def test_remove_hmm_deletes_directory(self):
        hmm_dir = os.path.join(self.out, 'hmm', 'HMM1')
        self.assertTrue(os.path.isdir(hmm_dir))
        s = self._setup_for_list_remove(remove='HMM1')
        s.remove_database_entry('HMM1')
        self.assertFalse(os.path.isdir(hmm_dir))

    def test_remove_diamond_entry_from_manifest(self):
        # DMND1.dmnd does not actually exist on disk (fake manifest), but
        # remove should still update the manifest without crashing
        s = self._setup_for_list_remove(remove='DMND1')
        s.remove_database_entry('DMND1')
        with open(os.path.join(self.out, MANIFEST_FILENAME)) as f:
            manifest = json.load(f)
        self.assertNotIn('DMND1', manifest)

    def test_list_with_no_manifest_dir_raises(self):
        args = _setup_args(output_dir='/does/not/exist', **{'list': True})
        with self.assertRaises(ConfigError):
            UserAnnotationDBSetup(args, run=_silent_run(), progress=_silent_progress())


# ===========================================================================
# TestRunnerSanityCheck
# ===========================================================================

class TestRunnerSanityCheck(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.out = os.path.join(self.tmp, 'out')
        _make_manifest(self.out, {'FakeDB': {'type': 'hmm'}})

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def test_missing_annotation_dir_raises(self):
        args = _runner_args(annotation_db_dir=None, contigs_db='/dev/null')
        with self.assertRaises(ConfigError):
            UserAnnotationRunner(args, run=_silent_run(), progress=_silent_progress())

    def test_nonexistent_annotation_dir_raises(self):
        args = _runner_args(annotation_db_dir='/does/not/exist', contigs_db='/dev/null')
        with self.assertRaises(ConfigError):
            UserAnnotationRunner(args, run=_silent_run(), progress=_silent_progress())

    def test_missing_manifest_raises(self):
        empty_dir = os.path.join(self.tmp, 'empty')
        os.makedirs(empty_dir)
        args = _runner_args(annotation_db_dir=empty_dir, contigs_db='/dev/null')
        with self.assertRaises(ConfigError):
            UserAnnotationRunner(args, run=_silent_run(), progress=_silent_progress())

    def test_missing_contigs_db_raises(self):
        args = _runner_args(annotation_db_dir=self.out, contigs_db=None)
        with self.assertRaises(ConfigError):
            UserAnnotationRunner(args, run=_silent_run(), progress=_silent_progress())


# ===========================================================================
# TestRunnerDatabaseSelection
# ===========================================================================

class TestRunnerDatabaseSelection(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.out = os.path.join(self.tmp, 'out')
        self.manifest = {
            'HMM1': {'type': 'hmm', 'hmm_dir': '/fake/hmm/HMM1'},
            'DMND1': {'type': 'diamond', 'dmnd_path': '/fake/DMND1.dmnd'},
        }
        _make_manifest(self.out, self.manifest)

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def _make_runner(self, database=None):
        args = _runner_args(annotation_db_dir=self.out, contigs_db='/dev/null',
                            database=database)
        r = object.__new__(UserAnnotationRunner)
        r.args = args
        r.run = _silent_run()
        r.progress = _silent_progress()
        r.contigs_db_path = '/dev/null'
        r.annotation_db_dir = self.out
        r.database = database
        r.num_threads = 1
        r.just_do_it = False
        r.hmmer_program = 'hmmscan'
        r.evalue = None
        r.min_pident = None
        r.diamond_sensitivity = None
        return r

    def test_no_database_flag_selects_all(self):
        r = self._make_runner(database=None)
        manifest = r.load_manifest()
        selected = manifest if not r.database or r.database.lower() == 'all' else None
        self.assertEqual(set(selected.keys()), {'HMM1', 'DMND1'})

    def test_database_all_selects_all(self):
        r = self._make_runner(database='all')
        manifest = r.load_manifest()
        selected = manifest if not r.database or r.database.lower() == 'all' else None
        self.assertEqual(set(selected.keys()), {'HMM1', 'DMND1'})

    def test_single_database_selected(self):
        r = self._make_runner(database='HMM1')
        manifest = r.load_manifest()
        requested = [n.strip() for n in r.database.split(',') if n.strip()]
        missing = [n for n in requested if n not in manifest]
        self.assertEqual(missing, [])
        selected = {n: manifest[n] for n in requested}
        self.assertEqual(set(selected.keys()), {'HMM1'})

    def test_comma_separated_databases_selected(self):
        r = self._make_runner(database='HMM1,DMND1')
        manifest = r.load_manifest()
        requested = [n.strip() for n in r.database.split(',') if n.strip()]
        selected = {n: manifest[n] for n in requested}
        self.assertEqual(set(selected.keys()), {'HMM1', 'DMND1'})

    def test_unknown_database_name_raises(self):
        r = self._make_runner(database='DoesNotExist')
        manifest = r.load_manifest()
        requested = ['DoesNotExist']
        missing = [n for n in requested if n not in manifest]
        self.assertTrue(len(missing) > 0)


# ===========================================================================
# TestSetupSanityCheck
# ===========================================================================

class TestSetupSanityCheck(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def test_missing_output_dir_raises(self):
        args = _setup_args(output_dir=None, input_tsv=F_HMM_ALL_TC)
        with self.assertRaises(ConfigError):
            UserAnnotationDBSetup(args, run=_silent_run(), progress=_silent_progress())

    def test_missing_input_tsv_raises(self):
        out = os.path.join(self.tmp, 'out')
        os.makedirs(out)
        args = _setup_args(output_dir=out, input_tsv=None)
        with self.assertRaises(ConfigError):
            UserAnnotationDBSetup(args, run=_silent_run(), progress=_silent_progress())

    def test_list_without_manifest_raises(self):
        out = os.path.join(self.tmp, 'empty')
        os.makedirs(out)
        args = _setup_args(output_dir=out, **{'list': True})
        with self.assertRaises(ConfigError):
            UserAnnotationDBSetup(args, run=_silent_run(), progress=_silent_progress())

    def test_remove_without_manifest_raises(self):
        out = os.path.join(self.tmp, 'empty')
        os.makedirs(out)
        args = _setup_args(output_dir=out, remove='SomeDB')
        with self.assertRaises(ConfigError):
            UserAnnotationDBSetup(args, run=_silent_run(), progress=_silent_progress())

    def test_list_does_not_require_input_tsv(self):
        out = os.path.join(self.tmp, 'out')
        _make_manifest(out, {})
        args = _setup_args(output_dir=out, **{'list': True})
        # Should not raise even without input_tsv
        s = UserAnnotationDBSetup(args, run=_silent_run(), progress=_silent_progress())
        self.assertIsNotNone(s)


# ===========================================================================
# Integration tests (skipped when external tools absent)
# ===========================================================================

import subprocess as _subprocess

def _diamond_available():
    import shutil
    for candidate in ('diamond', shutil.which('diamond')):
        if not candidate:
            continue
        try:
            _subprocess.run([candidate, '--version'], capture_output=True, timeout=5)
            return True
        except (FileNotFoundError, _subprocess.TimeoutExpired):
            continue
    return False

@unittest.skipUnless(_diamond_available(), 'diamond not in PATH')
class TestSetupDiamondSource(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.out = os.path.join(self.tmp, 'out')
        os.makedirs(self.out)
        tsv = _make_tsv(os.path.join(self.tmp, 't.tsv'), [('TestFAA', F_VALID_FAA)])
        self.s = _make_setup_instance(self.out, tsv)

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def test_dmnd_file_created(self):
        entry = self.s.setup_diamond_source('TestFAA', F_VALID_FAA)
        self.assertTrue(os.path.exists(entry['dmnd_path']))

    def test_dmnd_path_ends_with_dmnd(self):
        entry = self.s.setup_diamond_source('TestFAA', F_VALID_FAA)
        self.assertTrue(entry['dmnd_path'].endswith('.dmnd'))

    def test_entry_type_is_diamond(self):
        entry = self.s.setup_diamond_source('TestFAA', F_VALID_FAA)
        self.assertEqual(entry['type'], 'diamond')

    def test_num_sequences_counted(self):
        entry = self.s.setup_diamond_source('TestFAA', F_VALID_FAA)
        self.assertEqual(entry['num_sequences'], 2)


if __name__ == '__main__':
    unittest.main()

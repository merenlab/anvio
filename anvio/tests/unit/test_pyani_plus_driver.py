"""Unit coverage for the pyANI-plus matrix adapter and backend policy."""

import os
import tempfile
import unittest
from argparse import Namespace
from unittest import mock

from anvio.drivers.pyani_plus import PyANIPlus
from anvio.drivers.pyani import PyANI
from anvio.drivers.pyani_methods import validate_ani_method
from anvio.errors import ConfigError
from anvio.genomesimilarity import ANI, Dereplicate, FastANI, GenomeSimilarity, SourMash


class PyANIPlusDriverTestCase(unittest.TestCase):
    def test_matrix_reader_keeps_undefined_comparisons_as_none(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            matrix_path = os.path.join(temp_dir, 'identity.tsv')
            with open(matrix_path, 'w') as matrix_file:
                matrix_file.write('\talpha\tbeta\n')
                matrix_file.write('alpha\t1.0\tNA\n')
                matrix_file.write('beta\t0.98\t1.0\n')

            matrix = PyANIPlus._read_matrix(matrix_path)

        self.assertEqual(matrix['alpha']['alpha'], 1.0)
        self.assertIsNone(matrix['alpha']['beta'])
        self.assertEqual(matrix['beta']['alpha'], 0.98)

    def test_matrix_reader_rejects_row_column_label_mismatch(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            matrix_path = os.path.join(temp_dir, 'bad.tsv')
            with open(matrix_path, 'w') as matrix_file:
                matrix_file.write('\talpha\tbeta\n')
                matrix_file.write('alpha\t1.0\t0.5\n')
                matrix_file.write('gamma\t0.5\t1.0\n')

            with self.assertRaises(ConfigError):
                PyANIPlus._read_matrix(matrix_path)

    def test_adapter_runs_isolated_program_and_imports_exported_matrices(self):
        args = Namespace(method='ANIb', num_threads=1, pyani_plus_program='/isolated/bin/pyani-plus',
                         log_file='/tmp/pyani-plus-test.log', quiet=True)
        with mock.patch.object(PyANIPlus, '_get_cpu_allocation', return_value=('/usr/bin/taskset', [0, 1], [0])), \
                mock.patch('anvio.drivers.pyani_plus.utils.is_program_exists', return_value='/isolated/bin/pyani-plus'):
            driver = PyANIPlus(args)

        with tempfile.TemporaryDirectory() as input_dir:
            with open(os.path.join(input_dir, 'alpha.fa'), 'w') as fasta:
                fasta.write('>alpha\nACGT\n')
            with open(os.path.join(input_dir, 'beta.fa'), 'w') as fasta:
                fasta.write('>beta\nACGA\n')

            commands = []

            def fake_run_command(command, log_path, **kwargs):
                commands.append((command, kwargs))
                if 'export-run' in command:
                    export_dir = command[command.index('--outdir') + 1]
                    for suffix, values in {
                            'identity': 'alpha\t1.0\t0.9\nbeta\t0.9\t1.0\n',
                            'query_cov': 'alpha\t1.0\t0.5\nbeta\t0.75\t1.0\n',
                            'aln_lengths': 'alpha\t4\t2\nbeta\t3\t4\n',
                            'sim_errors': 'alpha\t0\t1\nbeta\t1\t0\n',
                            'hadamard': 'alpha\t1.0\t0.45\nbeta\t0.675\t1.0\n',
                    }.items():
                        with open(os.path.join(export_dir, 'ANIb_%s.tsv' % suffix), 'w') as matrix:
                            matrix.write('\talpha\tbeta\n%s' % values)
                return 0

            with mock.patch('anvio.drivers.pyani_plus.utils.is_program_exists', return_value='/isolated/bin/pyani-plus'), \
                    mock.patch('anvio.drivers.pyani_plus.utils.run_command', side_effect=fake_run_command):
                matrices = driver.run_command(input_dir)

        self.assertEqual(matrices['percentage_identity']['alpha']['beta'], 0.9)
        self.assertEqual(matrices['alignment_coverage']['alpha']['beta'], 0.5)
        self.assertEqual(matrices['alignment_coverage']['beta']['alpha'], 0.75)
        self.assertEqual(matrices['alignment_lengths']['alpha']['beta'], 2.0)
        self.assertEqual(matrices['similarity_errors']['alpha']['beta'], 1.0)
        self.assertEqual(commands[0][0][:3], ['/usr/bin/taskset', '-c', '0'])
        self.assertIn('PATH=/isolated/bin:' + os.environ.get('PATH', ''), commands[0][0][4])
        self.assertEqual(commands[1][1]['remove_log_file_if_exists'], False)

    def test_progress_ends_when_command_raises(self):
        args = Namespace(method='ANIb', num_threads=1, pyani_plus_program='/isolated/bin/pyani-plus',
                         log_file='/tmp/pyani-plus-test.log', quiet=True)
        progress = mock.Mock()
        with mock.patch.object(PyANIPlus, '_get_cpu_allocation', return_value=('/usr/bin/taskset', [0], [0])), \
                mock.patch('anvio.drivers.pyani_plus.utils.is_program_exists', return_value='/isolated/bin/pyani-plus'):
            driver = PyANIPlus(args, progress=progress)

        with tempfile.TemporaryDirectory() as input_dir:
            with open(os.path.join(input_dir, 'alpha.fa'), 'w') as fasta:
                fasta.write('>alpha\nACGT\n')
            with mock.patch('anvio.drivers.pyani_plus.utils.run_command', side_effect=RuntimeError('launcher failed')):
                with self.assertRaisesRegex(RuntimeError, 'launcher failed'):
                    driver.run_command(input_dir)

        progress.end.assert_called_once_with()

    def test_unsupported_method_is_rejected(self):
        args = Namespace(method='TETRA', num_threads=1, pyani_plus_program='pyani-plus',
                         log_file='/tmp/pyani-plus-test.log', quiet=True)
        with self.assertRaisesRegex(ConfigError, 'TETRA.*retired'):
            PyANIPlus(args)

        args = Namespace(method='TETRA', num_threads=1, log_file='/tmp/pyani-plus-test.log', quiet=True)
        with self.assertRaisesRegex(ConfigError, 'TETRA.*retired'):
            PyANI(args)

    def test_legacy_driver_rejects_inconsistent_backend_api_options(self):
        with self.assertRaisesRegex(ConfigError, 'Unknown ANI backend'):
            PyANI(Namespace(ani_backend=''))
        with self.assertRaisesRegex(ConfigError, 'cannot be used with ANI backend'):
            PyANI(Namespace(ani_backend='pyani-plus'))
        with self.assertRaisesRegex(ConfigError, '--pyani-plus-program.*legacy PyANI driver'):
            PyANI(Namespace(pyani_plus_program='/custom/pyani-plus'))

    def test_non_pyani_similarity_apis_reject_pyani_backend_options(self):
        with self.assertRaisesRegex(ConfigError, 'only be used together with --program pyANI'):
            FastANI(Namespace(ani_backend='legacy'))
        with self.assertRaisesRegex(ConfigError, 'only be used with --program pyANI'):
            FastANI(Namespace(pyani_plus_program='/custom/pyani-plus'))
        with self.assertRaisesRegex(ConfigError, 'only be used together with --program pyANI'):
            SourMash(Namespace(ani_backend='pyani-plus'))
        with self.assertRaisesRegex(ConfigError, 'only be used with --program pyANI'):
            SourMash(Namespace(pyani_plus_program='/custom/pyani-plus'))

    def test_cpu_allocation_info_reports_requested_and_actual_cpus(self):
        args = Namespace(method='ANIb', num_threads=1, pyani_plus_program='pyani-plus',
                         log_file='/tmp/pyani-plus-test.log', quiet=True)
        run = mock.Mock()
        with mock.patch.object(PyANIPlus, '_get_cpu_allocation', return_value=('/usr/bin/taskset', [3, 8], [3])), \
                mock.patch('anvio.drivers.pyani_plus.utils.is_program_exists', return_value='/isolated/bin/pyani-plus'):
            PyANIPlus(args, run=run, progress=mock.Mock())
        run.warning.assert_not_called()
        self.assertIn('Requested 1; detected 2 allowed CPUs; allocated 1 CPUs 3', str(run.info.call_args_list))
        self.assertIn('does not set an OS thread limit', str(run.info.call_args_list))

    def test_missing_executable_recommends_required_dependency_reinstall(self):
        args = Namespace(method='ANIb', num_threads=1, pyani_plus_program='pyani-plus',
                         log_file='/tmp/pyani-plus-test.log', quiet=True)
        with mock.patch.object(PyANIPlus, '_get_cpu_allocation', return_value=('/usr/bin/taskset', [0], [0])), \
                mock.patch('anvio.drivers.pyani_plus.utils.is_program_exists', return_value=False):
            with self.assertRaises(ConfigError) as error:
                PyANIPlus(args)
        message = str(error.exception)
        self.assertIn('python -m pip install anvio', message)
        self.assertIn('legacy pyANI executable', message)
        self.assertIn('older-Python environment', ' '.join(message.split()))


class PyANIPlusCPUAllocationTestCase(unittest.TestCase):
    def test_api_budget_defaults_only_when_missing_or_none(self):
        for include_budget, raw in ((False, None), (True, None)):
            with self.subTest(raw=raw):
                args = Namespace(method='ANIb', pyani_plus_program='/isolated/bin/pyani-plus',
                                 log_file='/tmp/pyani-plus-test.log', quiet=True)
                if include_budget:
                    args.num_threads = raw
                with mock.patch.object(PyANIPlus, '_get_cpu_allocation', return_value=('/usr/bin/taskset', [0, 1], [0])), \
                        mock.patch('anvio.drivers.pyani_plus.utils.is_program_exists', return_value='/isolated/bin/pyani-plus'):
                    driver = PyANIPlus(args)
                self.assertEqual(driver.num_threads, 1)

    def test_api_rejects_zero_boolean_and_non_integer_budgets(self):
        for raw in (0, -1, True, 1.5, '2'):
            with self.subTest(raw=raw):
                args = Namespace(method='ANIb', num_threads=raw, pyani_plus_program='/isolated/bin/pyani-plus',
                                 log_file='/tmp/pyani-plus-test.log', quiet=True)
                with self.assertRaisesRegex(ConfigError, 'positive integer'):
                    PyANIPlus(args)

    def test_api_rejects_unknown_and_incompatible_backend_values(self):
        for backend in ('', 'unknown'):
            with self.subTest(backend=backend):
                with self.assertRaisesRegex(ConfigError, 'Unknown ANI backend'):
                    ANI(Namespace(ani_backend=backend))

        with self.assertRaisesRegex(ConfigError, 'Unknown ANI backend'):
            PyANIPlus(Namespace(ani_backend=''))

        with self.assertRaisesRegex(ConfigError, 'cannot be used with ANI backend'):
            PyANIPlus(Namespace(ani_backend='legacy'))

    def test_non_contiguous_allowed_cpu_ids_are_selected_and_request_is_capped(self):
        taskset, allowed, allocated = PyANIPlus._get_cpu_allocation(
            2, system='Linux', get_affinity=lambda pid: {3, 8, 12}, taskset_path='/usr/bin/taskset')
        self.assertEqual(taskset, '/usr/bin/taskset')
        self.assertEqual(allowed, [3, 8, 12])
        self.assertEqual(allocated, [3, 8])

        taskset, allowed, allocated = PyANIPlus._get_cpu_allocation(
            20, system='Linux', get_affinity=lambda pid: {3, 8, 12}, taskset_path='/usr/bin/taskset')
        self.assertEqual(allocated, allowed)

    def test_invalid_budgets_fail_before_platform_or_tool_lookup(self):
        for value in (0, -1, 1.5, True, '2'):
            with self.subTest(value=value):
                with self.assertRaisesRegex(ConfigError, 'positive integer'):
                    PyANIPlus._get_cpu_allocation(
                        value, system='Windows', get_affinity=None, taskset_path=None)

    def test_non_linux_and_missing_affinity_support_fail_clearly(self):
        with self.assertRaisesRegex(ConfigError, 'only on Linux'):
            PyANIPlus._get_cpu_allocation(1, system='Darwin', get_affinity=lambda pid: {0}, taskset_path='taskset')

        with self.assertRaisesRegex(ConfigError, 'sched_getaffinity'):
            PyANIPlus._get_cpu_allocation(1, system='Linux', get_affinity=None, taskset_path='taskset')

    def test_missing_taskset_and_empty_allocation_fail_clearly(self):
        with self.assertRaisesRegex(ConfigError, 'not found'):
            PyANIPlus._get_cpu_allocation(1, system='Linux', get_affinity=lambda pid: {0}, taskset_path=False)

        with self.assertRaisesRegex(ConfigError, 'no CPUs'):
            PyANIPlus._get_cpu_allocation(1, system='Linux', get_affinity=lambda pid: set(), taskset_path='taskset')

    def test_retired_methods_have_specific_migration_guidance(self):
        with self.assertRaises(ConfigError) as error:
            validate_ani_method('ANIblastall')
        self.assertIn('ANIblastall', str(error.exception))
        self.assertIn('retired', str(error.exception))
        self.assertIn('ANIb', str(error.exception))
        self.assertIn('not guaranteed', str(error.exception))

        with self.assertRaises(ConfigError) as error:
            validate_ani_method('TETRA')
        self.assertIn('TETRA', str(error.exception))
        self.assertIn('retired', str(error.exception))
        self.assertIn('no equivalent', str(error.exception))
        self.assertIsNone(validate_ani_method('ANIm'))
        self.assertIsNone(validate_ani_method('ANIb'))

    def test_unsupported_method_is_rejected(self):
        with self.assertRaisesRegex(ConfigError, "Unsupported ANI method 'not-a-method'"):
            validate_ani_method('not-a-method')

    def test_dereplication_rejects_imported_ani_matrices_with_undefined_values(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            identity_path = os.path.join(temp_dir, 'percentage_identity.txt')
            coverage_path = os.path.join(temp_dir, 'alignment_coverage.txt')
            matrix_text = 'key\talpha\tbeta\nalpha\t1.0\t\nbeta\t\t1.0\n'
            for path in [identity_path, coverage_path]:
                with open(path, 'w') as matrix_file:
                    matrix_file.write(matrix_text)

            dereplicate = Dereplicate.__new__(Dereplicate)
            dereplicate.program_name = 'pyANI'
            dereplicate.ani_dir = temp_dir
            dereplicate.mash_dir = None
            dereplicate.program_info = {
                'metric_name': 'percentage_identity',
                'necessary_reports': ['percentage_identity', 'alignment_coverage'],
            }
            dereplicate.similarity = type('Similarity', (), {'results': {}})()

            with self.assertRaisesRegex(ConfigError, 'cannot use incomplete ANI matrices'):
                dereplicate.import_similarity_matrix()

            with open(identity_path) as matrix_file:
                self.assertEqual(matrix_file.read(), matrix_text)


    def _generated_dereplicate(self, metric='percentage_identity'):
        dereplicate = Dereplicate.__new__(Dereplicate)
        dereplicate.program_name = 'pyANI'
        dereplicate.program_info = {
            'metric_name': metric,
            'necessary_reports': [metric, 'alignment_coverage'],
        }
        dereplicate.temp_dir = 'temporary-genome-inputs'
        dereplicate.similarity = mock.Mock()
        dereplicate.similarity.results = {
            report: {'alpha': {'alpha': 1.0, 'beta': 0.98},
                     'beta': {'alpha': 0.97, 'beta': 1.0}}
            for report in dereplicate.program_info['necessary_reports']
        }
        return dereplicate


    def test_generated_ani_rejects_undefined_values_in_both_directions_and_diagonal(self):
        for value in (None, float('nan'), float('inf'), float('-inf'), ''):
            for row, column in (('alpha', 'beta'), ('beta', 'alpha'), ('alpha', 'alpha')):
                with self.subTest(value=value, row=row, column=column):
                    dereplicate = self._generated_dereplicate()
                    matrix = dereplicate.similarity.results['percentage_identity']
                    matrix[row][column] = value
                    with self.assertRaisesRegex(ConfigError, '(?s)generated ANI matrix.*cannot use incomplete ANI matrices'):
                        dereplicate.gen_similarity_matrix()
                    dereplicate.similarity.process.assert_called_once_with(dereplicate.temp_dir)
                    self.assertIs(matrix[row][column], value)


    def test_generated_full_identity_and_alignment_coverage_are_validated(self):
        for report in ('full_percentage_identity', 'alignment_coverage'):
            with self.subTest(report=report):
                dereplicate = self._generated_dereplicate('full_percentage_identity')
                dereplicate.similarity.results[report]['beta']['alpha'] = None
                with self.assertRaisesRegex(ConfigError, report):
                    dereplicate.gen_similarity_matrix()


    def test_complete_generated_ani_matrix_is_returned_unchanged(self):
        dereplicate = self._generated_dereplicate()
        matrix = dereplicate.similarity.results['percentage_identity']
        self.assertIs(dereplicate.gen_similarity_matrix(), matrix)


    def test_generated_incomplete_ani_is_rejected_before_cluster_mutation(self):
        dereplicate = self._generated_dereplicate()
        dereplicate.sequence_source_provided = False
        dereplicate.output_dir = 'output'
        dereplicate.similarity.results['percentage_identity']['beta']['alpha'] = None
        dereplicate.import_previous_results = False
        dereplicate.clusters = {'sentinel': {'alpha'}}
        with mock.patch.object(dereplicate, 'init_output_dir'), \
                mock.patch.object(dereplicate, 'init_genome_similarity'), \
                mock.patch.object(dereplicate, 'init_clusters') as init_clusters, \
                mock.patch.object(dereplicate, 'dereplicate') as cluster:
            with self.assertRaises(ConfigError):
                dereplicate.process()
        self.assertEqual(dereplicate.clusters, {'sentinel': {'alpha'}})
        init_clusters.assert_not_called()
        cluster.assert_not_called()


    def test_incomplete_ani_releases_owned_temporary_fasta_directory(self):
        for imported in (False, True):
            with self.subTest(imported=imported), tempfile.TemporaryDirectory() as parent:
                dereplicate = self._generated_dereplicate()
                dereplicate.temp_dir = os.path.join(parent, 'genomes')
                os.mkdir(dereplicate.temp_dir)
                with open(os.path.join(dereplicate.temp_dir, 'alpha.fa'), 'w') as fasta:
                    fasta.write('>alpha\nACGT\n')
                owned_path = dereplicate.temp_dir
                dereplicate.import_previous_results = imported
                dereplicate.similarity.results['percentage_identity']['beta']['alpha'] = None
                importer = mock.patch.object(dereplicate, 'import_similarity_matrix',
                                             side_effect=ConfigError('Incomplete imported ANI matrix'))
                with importer, mock.patch('anvio.DEBUG', False):
                    with self.assertRaises(ConfigError):
                        dereplicate.get_similarity_matrix()
                self.assertFalse(os.path.exists(owned_path))
                self.assertIsNone(dereplicate.temp_dir)
                self.assertTrue(os.path.isdir(parent))


    def test_incomplete_ani_preserves_temporary_fasta_directory_in_debug_mode(self):
        with tempfile.TemporaryDirectory() as parent:
            dereplicate = self._generated_dereplicate()
            dereplicate.temp_dir = os.path.join(parent, 'genomes')
            os.mkdir(dereplicate.temp_dir)
            owned_path = dereplicate.temp_dir
            dereplicate.import_previous_results = False
            dereplicate.similarity.results['percentage_identity']['beta']['alpha'] = None
            with mock.patch('anvio.DEBUG', True):
                with self.assertRaises(ConfigError):
                    dereplicate.get_similarity_matrix()
            self.assertTrue(os.path.isdir(owned_path))
            self.assertEqual(dereplicate.temp_dir, owned_path)


    def test_generated_fastani_does_not_apply_pyani_validation(self):
        dereplicate = self._generated_dereplicate('ani')
        dereplicate.program_name = 'fastANI'
        matrix = dereplicate.similarity.results['ani']
        self.assertIs(dereplicate.gen_similarity_matrix(), matrix)


    def test_pan_db_rejects_undefined_pyani_plus_values_before_mutation(self):
        ani = ANI.__new__(ANI)
        ani.ani_backend = 'pyani-plus'
        ani.pan_db = 'pan.db'
        ani.output_dir = 'similarity-output'
        ani.results = {'percentage_identity': {'alpha': {'alpha': 1.0, 'beta': None},
                                               'beta': {'alpha': None, 'beta': 1.0}}}

        with mock.patch.object(GenomeSimilarity, 'add_to_pan_db') as add_to_pan_db:
            with self.assertRaises(ConfigError) as error:
                ani.add_to_pan_db()

        self.assertIn('similarity-output', str(error.exception))
        add_to_pan_db.assert_not_called()


if __name__ == '__main__':
    unittest.main()

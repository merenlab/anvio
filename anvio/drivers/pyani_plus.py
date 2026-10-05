"""Adapter from pyANI-plus exports to anvi'o's legacy ANI matrix format."""

import csv
import os
import platform
import shutil
import tempfile

import anvio
import anvio.filesnpaths as filesnpaths
import anvio.terminal as terminal
import anvio.utils as utils

from anvio.errors import ConfigError
from anvio.drivers.pyani_methods import validate_ani_method

_PROCESS_CPU_AFFINITY = object()


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://www.anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__


class PyANIPlus:
    """Run pyANI-plus in an isolated executable environment and import its TSVs.

    pyANI-plus writes identity and query coverage as fractions, matching the
    scale expected by the existing anvi'o ANI code. Its query coverage is
    directional, and its alignment lengths are left undefined when no
    alignment is available. Those values are retained as ``None``.
    """

    MATRIX_FILENAMES = {
        'percentage_identity': 'identity',
        'alignment_coverage': 'query_cov',
        'alignment_lengths': 'aln_lengths',
        'similarity_errors': 'sim_errors',
        'hadamard': 'hadamard',
    }

    def __init__(self, args={}, run=terminal.Run(), progress=terminal.Progress()):
        self.run = run
        self.progress = progress

        A = lambda x: args.__dict__.get(x)
        self.method = A('method') or 'ANIb'
        raw_cpu_budget = A('num_threads')
        self.num_threads = 1 if raw_cpu_budget is None else raw_cpu_budget
        raw_backend = A('ani_backend')
        if raw_backend is not None and raw_backend not in ('legacy', 'pyani-plus'):
            raise ConfigError(f"Unknown ANI backend '{raw_backend}'. Choose 'legacy' or 'pyani-plus'.")
        if raw_backend == 'legacy':
            raise ConfigError(f"The pyANI-plus driver cannot be used with ANI backend '{raw_backend}'.")
        raw_program_name = A('pyani_plus_program')
        if raw_program_name is not None and not raw_program_name:
            raise ConfigError("The pyANI-plus executable name must not be empty.")
        self.program_name = raw_program_name or 'pyani-plus'
        self.log_file_path = os.path.abspath(A('log_file') or filesnpaths.get_temp_file_path())
        self.quiet = A('quiet')

        validate_ani_method(self.method)

        self.taskset_path, self.allowed_cpu_ids, self.allocated_cpu_ids = self._get_cpu_allocation(self.num_threads)

        self.program_path = utils.is_program_exists(self.program_name, dont_raise=True)
        if not self.program_path:
            raise ConfigError("The pyani-plus backend needs the 'pyani-plus' executable, but '%s' was not found. "
                              "Correct `--pyani-plus-program` or reinstall anvi'o and its required dependencies with "
                              "`python -m pip install anvio`, "
                              "or use `--ani-backend legacy` only with a separately validated legacy pyANI executable in an older-Python environment; legacy PyANI is not supported in this Python 3.13 environment." % self.program_name)
        self.program_bin_dir = os.path.dirname(self.program_path)

        self.run.info('[pyANI-plus] Executable', self.program_path)
        self.run.info('[pyANI-plus] Alignment method', self.method)
        self.run.info(
            '[pyANI-plus] CPU allocation',
            'Requested %d; detected %d allowed CPUs; allocated %d CPUs %s with taskset. '
            'pyANI-plus manages worker concurrency within this CPU allocation; this does not set an OS thread limit.' % (
                self.num_threads, len(self.allowed_cpu_ids), len(self.allocated_cpu_ids),
                ','.join(map(str, self.allocated_cpu_ids)))
        )
        self.run.info('[pyANI-plus] Log file path', self.log_file_path, nl_after=1)


    @staticmethod
    def _get_cpu_allocation(num_threads, system=None, get_affinity=_PROCESS_CPU_AFFINITY, taskset_path=None):
        """Resolve a Linux taskset CPU allocation without changing this process."""
        if isinstance(num_threads, bool) or not isinstance(num_threads, int) or num_threads < 1:
            raise ConfigError("The pyANI-plus CPU allocation must be a positive integer; received %r." % num_threads)

        system = platform.system() if system is None else system
        if system != 'Linux':
            raise ConfigError("pyANI-plus CPU allocation through taskset is supported only on Linux.")

        if get_affinity is _PROCESS_CPU_AFFINITY:
            get_affinity = getattr(os, 'sched_getaffinity', None)
        if get_affinity is None:
            raise ConfigError("Cannot allocate CPUs for pyANI-plus: this Python build does not provide os.sched_getaffinity.")

        try:
            allowed_cpu_ids = sorted(get_affinity(0))
        except (OSError, ProcessLookupError) as e:
            raise ConfigError("Cannot read the CPUs currently available to pyANI-plus: %s" % e)
        if not allowed_cpu_ids:
            raise ConfigError("Cannot allocate CPUs for pyANI-plus: no CPUs are currently available to this process.")

        taskset_path = shutil.which('taskset') if taskset_path is None else taskset_path
        if not taskset_path:
            raise ConfigError("Cannot allocate CPUs for pyANI-plus because the Linux `taskset` executable was not found. "
                              "Install the system package that provides `taskset` (commonly util-linux).")

        allocated_cpu_ids = allowed_cpu_ids[:min(num_threads, len(allowed_cpu_ids))]
        return taskset_path, allowed_cpu_ids, allocated_cpu_ids

    def _command(self, *args):
        """Use the selected executable's bin directory for its helper scripts."""
        current_path = os.environ.get('PATH', '')
        isolated_path = os.pathsep.join([self.program_bin_dir, current_path])
        cpu_list = ','.join(map(str, self.allocated_cpu_ids))
        return [self.taskset_path, '-c', cpu_list, 'env', 'PATH=' + isolated_path, self.program_path, *args]

    @staticmethod
    def _read_matrix(path):
        """Read a pyANI-plus matrix while preserving undefined values as None."""
        matrix = {}
        try:
            with open(path, newline='') as matrix_file:
                rows = csv.reader(matrix_file, delimiter='\t')
                header = next(rows)
                if len(header) < 2:
                    raise ValueError("matrix has no column labels")
                column_names = header[1:]
                if len(column_names) != len(set(column_names)):
                    raise ValueError("matrix has duplicate column labels")

                for row in rows:
                    if len(row) != len(column_names) + 1:
                        raise ValueError("matrix row has an unexpected number of values")
                    name = row[0]
                    if name in matrix:
                        raise ValueError("matrix has duplicate row labels")
                    matrix[name] = {}
                    for column, value in zip(column_names, row[1:]):
                        if value.strip().lower() in ['', 'na', 'nan', 'none']:
                            matrix[name][column] = None
                        else:
                            matrix[name][column] = float(value)

                if set(matrix) != set(column_names):
                    raise ValueError("matrix row and column labels do not match")
        except (OSError, StopIteration, ValueError) as e:
            raise ConfigError("Could not import pyANI-plus matrix '%s': %s" % (path, e))

        return matrix

    def run_command(self, input_path):
        """Run one all-vs-all comparison and return anvi'o-compatible matrices."""
        input_path = os.path.abspath(input_path)
        fasta_paths = [os.path.join(input_path, f) for f in os.listdir(input_path)
                       if f.lower().endswith(('.fa', '.fas', '.fasta', '.fna'))]
        if not fasta_paths:
            raise ConfigError("No FASTA files were found for pyANI-plus in '%s'." % input_path)

        with tempfile.TemporaryDirectory(prefix='anvio-pyani-plus-', dir=input_path) as work_dir:
            database_path = os.path.join(work_dir, 'pyani-plus.sqlite')
            export_dir = os.path.join(work_dir, 'export')
            os.mkdir(export_dir)

            compute_command = self._command(
                self.method.lower(), input_path,
                '--database', database_path,
                '--create-db',
                '--executor', 'local',
            )
            self.progress.new('pyANI-plus')
            try:
                self.progress.update(f'Running {self.method} ...')
                exit_code = utils.run_command(compute_command, self.log_file_path)
                if int(exit_code):
                    raise ConfigError("pyANI-plus returned with non-zero exit code. Please check the log file for details.")

                export_command = self._command(
                    'export-run', '--database', database_path,
                    '--outdir', export_dir, '--run-id', '1', '--label', 'stem',
                )
                exit_code = utils.run_command(export_command, self.log_file_path,
                                              remove_log_file_if_exists=False)
                if int(exit_code):
                    raise ConfigError("pyANI-plus export-run returned with non-zero exit code. Please check the log file for details.")
            finally:
                self.progress.end()

            matrices = {}
            for matrix_name, file_suffix in self.MATRIX_FILENAMES.items():
                matrix_path = os.path.join(export_dir, '%s_%s.tsv' % (self.method, file_suffix))
                if os.path.exists(matrix_path):
                    matrices[matrix_name] = self._read_matrix(matrix_path)

        required = {'percentage_identity', 'alignment_coverage', 'alignment_lengths'}
        missing = required - set(matrices)
        if missing:
            raise ConfigError("pyANI-plus did not export the required matrix file(s): %s." % ', '.join(sorted(missing)))

        expected_labels = {os.path.splitext(os.path.basename(path))[0] for path in fasta_paths}
        for matrix_name, matrix in matrices.items():
            observed_labels = set(matrix)
            if observed_labels != expected_labels:
                missing_labels = sorted(expected_labels - observed_labels)
                extra_labels = sorted(observed_labels - expected_labels)
                raise ConfigError("pyANI-plus %s matrix genome labels do not match its input FASTA files. Missing: %s. Extra: %s." %
                                  (matrix_name, ', '.join(missing_labels) or 'none', ', '.join(extra_labels) or 'none'))

        if not self.quiet:
            self.run.info_single("pyANI-plus result matrices imported: %s." % ', '.join(sorted(matrices)),
                                 nl_before=1, mc='green')

        return matrices

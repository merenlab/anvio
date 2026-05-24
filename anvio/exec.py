"""Centralized execution engine for third-party tool invocations in anvi'o.

Two public objects live here:

  registry       — a lazy-loaded singleton that merges the three YAML files
                   under anvio/registry/ and resolves {version} placeholders
                   in runner strings.

  AnvioExecutor  — dispatches a third-party command to one of four backends
                   (native, conda, docker, singularity), selected by the
                   $ANVIO_RUNNER environment variable or an explicit
                   constructor argument.

Typical usage from a driver class:

    from anvio.exec import AnvioExecutor

    class Diamond:
        def makedb(self, query_fasta, run=run):
            executor = AnvioExecutor('diamond', run=run)
            result = executor.run_command(['diamond', 'makedb', '--in', query_fasta, ...])
            # result.stdout / result.stderr / result.returncode

To query metadata without running anything:

    from anvio.exec import registry

    tool = registry.get_tool('diamond')         # ToolEntry or None
    db   = registry.get_database('kegg-data')   # DatabaseEntry or None
    prog = registry.get_program('anvi-estimate-metabolism')  # ProgramEntry or None

    print(tool.version)                # '2.1.9'
    print(db.citation)                 # 'doi:10.1093/nar/gkw1092'
    print(prog.get_methods_text())     # rendered methods-section sentence

The backend falls back to 'native' if $ANVIO_RUNNER is not set.
Override it per-executor by passing backend= to AnvioExecutor(...).
"""

import os
import re
import yaml
import json
import shutil
import datetime
import subprocess

from threading import Lock
from collections import namedtuple

import anvio
import anvio.terminal as terminal

from anvio.errors import CommandError, ConfigError


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "A. Murat Eren"
__email__ = "a.murat.eren@gmail.com"


run      = terminal.Run()
progress = terminal.Progress()

BACKENDS = ('native', 'conda', 'docker', 'singularity')

# Resolved at module import from $ANVIO_RUNNER.
# Individual AnvioExecutor instances can override this further.
_env_backend = os.environ.get('ANVIO_RUNNER', 'native').lower()
DEFAULT_BACKEND = _env_backend if _env_backend in BACKENDS else 'native'

RunResult = namedtuple('RunResult', ['stdout', 'stderr', 'returncode', 'cmd'])


####################################################################################################
#
#     REGISTRY ENTRY CLASSES
#
####################################################################################################

class RegistryEntry:
    """Base class for all registry entries (tools, programs, databases).

    Subclasses expose domain-specific attributes. This base class owns
    get_methods_text(), which renders the optional example_text template
    defined in the YAML files.

    Attributes
    ==========
    name : str
        The registry key for this entry.
    citation : str
        The primary citation DOI, or 'PLACEHOLDER' if not yet filled in.
    """

    def __init__(self, name, citation='', example_text=''):
        self.name          = name
        self.citation      = citation
        self._example_text = example_text


    def __repr__(self):
        return f"{self.__class__.__name__}(name='{self.name}', citation='{self.citation}')"


    def get_methods_text(self):
        """Render the example_text template with attributes and citations resolved.

        Template tokens
        ===============
        {self}           : this entry's name
        {self.<attr>}    : any attribute on this entry (citation, description,
                           version, authors, …)
        {<key>.<attr>}   : any attribute of another registry entry, looked up
                           by its key across databases, programs, and tools

        Unknown tokens are left unchanged so partially-populated templates
        degrade gracefully rather than raising.

        Returns
        =======
        str
            The rendered sentence, or an empty string if no example_text is
            defined for this entry.

        Examples
        ========
        >>> prog = registry.get_program('anvi-estimate-metabolism')
        >>> print(prog.get_methods_text())
        Metabolic pathway completeness was estimated using anvi-estimate-metabolism
        (doi:10.1038/s41564-023-01439-2) against the KEGG database (doi:10.1093/nar/gkw1092).

        >>> tool = registry.get_tool('diamond')
        >>> # with a hypothetical example_text on the tool entry:
        >>> # "Sequences were aligned using {self} v{self.version} ({self.citation})."
        >>> print(tool.get_methods_text())
        Sequences were aligned using diamond v2.1.9 (doi:10.1038/s41592-021-01101-x).
        """
        if not self._example_text:
            return ''

        def resolve(match):
            token = match.group(1)
            if token == 'self':
                return self.name
            if token.startswith('self.'):
                attr = token[5:]
                value = getattr(self, attr, None)
                return str(value) if value is not None else f'{{{token}}}'
            if '.' in token:
                key, attr = token.rsplit('.', 1)
                entry = (registry.get_database(key)
                         or registry.get_program(key)
                         or registry.get_tool(key))
                if entry:
                    value = getattr(entry, attr, None)
                    return str(value) if value is not None else f'{{{token}}}'
            return f'{{{token}}}'

        return re.sub(r'\{([^}]+)\}', resolve, self._example_text)


class ProgramEntry(RegistryEntry):
    """Registry entry for an anvi'o program.

    Returned by Registry.get_program().

    Attributes
    ==========
    name : str
        The program name, e.g. 'anvi-estimate-metabolism'.
    authors : str
        Author(s) of the primary paper.
    citation : str
        DOI of that paper.
    uses : list
        Registry keys of the databases this program relies on.
    """

    def __init__(self, name, data):
        super().__init__(name,
                         citation=data.get('citation', ''),
                         example_text=data.get('example_text', ''))
        self.authors = data.get('authors', '')
        self.uses    = data.get('uses', [])


class DatabaseEntry(RegistryEntry):
    """Registry entry for a reference database.

    Returned by Registry.get_database().

    Attributes
    ==========
    name : str
        The database registry key, e.g. 'kegg-data'.
    authors : str
    citation : str
    provider : str
        'third-party' or 'anvio'.
    description : str
        One-sentence description of the database.
    url : str
    """

    def __init__(self, name, data):
        super().__init__(name,
                         citation=data.get('citation', ''),
                         example_text=data.get('example_text', ''))
        self.authors     = data.get('authors', '')
        self.provider    = data.get('provider', '')
        self.description = data.get('description', '')
        self.url         = data.get('url', '')


class ToolEntry(RegistryEntry):
    """Registry entry for a third-party tool.

    Returned by Registry.get_tool(). All {version} placeholders in runner
    strings are resolved at construction time, so runners always contains
    ready-to-use image/spec strings.

    Attributes
    ==========
    name : str
        The tool registry key, e.g. 'diamond'.
    authors : str
    citation : str
    url : str
    version : str
        The pinned version string, e.g. '2.1.9'.
    runners : dict
        Backend → resolved runner string, e.g.
        {'conda': 'bioconda::diamond=2.1.9', 'docker': 'buchfink/diamond:v2.1.9', …}
    """

    def __init__(self, name, data):
        super().__init__(name,
                         citation=data.get('citation', ''),
                         example_text=data.get('example_text', ''))
        self.authors = data.get('authors', '')
        self.url     = data.get('url', '')
        self.version = str(data.get('version', ''))

        # Resolve {version} in runner strings at construction time.
        self.runners = {
            backend: spec.replace('{version}', self.version) if spec and '{version}' in str(spec) else spec
            for backend, spec in data.get('runners', {}).items()
        }


####################################################################################################
#
#     REGISTRY
#
####################################################################################################

class Registry:
    """Loads and merges the three YAML registry files, providing entry lookups.

    The three files under anvio/registry/ are:

      third-party-programs.yaml  — external tools → ToolEntry instances
      anvio-programs.yaml        — anvi'o programs → ProgramEntry instances
      databases.yaml             — reference databases → DatabaseEntry instances

    The instance at module level (``registry``) is not populated until the
    first call to any lookup method, then cached for the lifetime of the
    process. Thread-safe via a double-checked lock.

    Notes
    =====
    - The anvi-* wildcard entry in anvio-programs.yaml is merged into every
      specific program entry as a set of defaults. Fields in a specific entry
      take precedence over the wildcard.
    - All three getter methods return None for unknown keys and never raise.
    """

    _REGISTRY_DIR = os.path.join(os.path.dirname(__file__), 'registry')
    _lock = Lock()


    def __init__(self):
        self._loaded  = False
        self.tools     = {}
        self.programs  = {}
        self.databases = {}


    def _ensure_loaded(self):
        if self._loaded:
            return
        with self._lock:
            if not self._loaded:
                self._load()


    def _load(self):
        def _read(filename):
            path = os.path.join(self._REGISTRY_DIR, filename)
            if not os.path.exists(path):
                raise ConfigError(f"Registry file not found: {path}. "
                                  f"This file should have shipped with anvi'o — "
                                  f"please report this as a bug at https://github.com/merenlab/anvio/issues.")
            with open(path) as f:
                return yaml.safe_load(f) or {}

        raw_tools = _read('third-party-programs.yaml')
        self.tools = {
            name: ToolEntry(name, data)
            for name, data in raw_tools.items()
            if isinstance(data, dict)
        }

        raw_programs = _read('anvio-programs.yaml')
        wildcard = raw_programs.pop('anvi-*', {})
        self.programs = {}
        for name, data in raw_programs.items():
            merged = dict(wildcard)
            if isinstance(data, dict):
                merged.update(data)
            self.programs[name] = ProgramEntry(name, merged)

        raw_databases = _read('databases.yaml')
        self.databases = {
            name: DatabaseEntry(name, data)
            for name, data in raw_databases.items()
            if isinstance(data, dict)
        }

        self._loaded = True


    def get_tool(self, name):
        """Return a ToolEntry for the named third-party tool, or None.

        Parameters
        ==========
        name : str
            Registry key, e.g. 'diamond', 'hmmer', 'prodigal'.

        Returns
        =======
        ToolEntry or None
        """
        self._ensure_loaded()
        return self.tools.get(name)


    def get_runner(self, name, backend):
        """Return the resolved runner string for a tool + backend pair, or None.

        For the 'native' backend this always returns None (no wrapper needed).
        For other backends the returned string has all {version} placeholders
        resolved and is ready to embed directly in a command.

        Parameters
        ==========
        name : str
            Tool registry key, e.g. 'diamond'.
        backend : str
            One of: 'conda', 'docker', 'singularity'.

        Returns
        =======
        str or None
        """
        self._ensure_loaded()
        entry = self.tools.get(name)
        return entry.runners.get(backend) if entry else None


    def get_database(self, name):
        """Return a DatabaseEntry for the named reference database, or None.

        Parameters
        ==========
        name : str
            Database key, e.g. 'kegg-data', 'cog-data', 'gtdb'.

        Returns
        =======
        DatabaseEntry or None
        """
        self._ensure_loaded()
        return self.databases.get(name)


    def get_program(self, name):
        """Return a ProgramEntry for the named anvi'o program, or None.

        The wildcard defaults from the 'anvi-*' entry are already merged in,
        so every returned entry carries at least 'authors' and 'citation'.
        Call get_methods_text() on the result to render the methods-section
        template for this program.

        Parameters
        ==========
        name : str
            Program name, e.g. 'anvi-estimate-metabolism'.

        Returns
        =======
        ProgramEntry or None
        """
        self._ensure_loaded()
        return self.programs.get(name)


# Module-level singleton — shared by all AnvioExecutor instances in a process.
registry = Registry()


####################################################################################################
#
#     ANVIO EXECUTOR
#
####################################################################################################

class AnvioExecutor:
    """Dispatch a third-party tool command to the appropriate execution backend.

    The backend is resolved in this priority order:
      1. Explicit ``backend`` constructor argument.
      2. ANVIO_RUNNER environment variable.
      3. 'native' (run the binary directly from PATH).

    Supported backends
    ==================
    native      : run the binary from PATH (the default, zero overhead).
    docker      : wrap with ``docker run``, mounting work_dir at the same
                  absolute path as on the host (Nextflow-style bind mount).
    singularity : wrap with ``singularity exec``, binding work_dir.
    conda       : wrap with ``conda run`` inside a per-tool conda environment
                  whose name is derived from the registry conda spec.

    Parameters
    ==========
    tool_name : str
        Registry key for the third-party tool, e.g. 'diamond', 'hmmer'.
    backend : str, optional
        Override the execution backend. One of: native, docker, singularity, conda.
    work_dir : str, optional
        Working directory for container-based backends. Defaults to cwd.
        For docker and singularity, this directory is bind-mounted at the
        same absolute path inside the container, so file arguments in cmd
        resolve identically inside and outside the container.
    log : str, optional
        Path to a log file. When provided, the executor writes the command,
        backend, stdout, stderr, and return code to this file after every
        run_command() call — even when the command fails.
    run : terminal.Run, optional
    progress : terminal.Progress, optional

    Notes
    =====
    - Tools not in the registry still work with 'native'. For other backends,
      a missing registry entry raises ConfigError.
    - The conda backend expects a pre-created environment named after the tool
      spec (e.g. anvio-tool-diamond-2.1.9). Use get_conda_env_name() to
      get the expected name, then create it with:
          conda create -n <name> -c bioconda -c conda-forge <package>=<version>
    - For docker and singularity, HASH_PLACEHOLDER entries in the registry must
      be filled in before those backends are usable. Pending hashes are present
      in anvio/registry/third-party-programs.yaml.

    Examples
    ========
    >>> executor = AnvioExecutor('diamond', run=self.run, log=self.run.log_file_path)
    >>> result = executor.run_command(['diamond', 'blastp', '-q', 'in.fa', '-d', 'db'])
    >>> result.returncode
    0

    >>> # Override the backend for a single invocation:
    >>> executor = AnvioExecutor('diamond', backend='docker', work_dir='/data/run1')
    >>> result = executor.run_command(['diamond', 'blastp', ...])
    """

    def __init__(self, tool_name, backend=None, work_dir=None, log=None, run=run, progress=progress):
        self.tool_name = tool_name
        self.run       = run
        self.progress  = progress
        self.work_dir  = os.path.abspath(work_dir or os.getcwd())
        self.log       = log

        self.backend = (backend or DEFAULT_BACKEND).lower()
        if self.backend not in BACKENDS:
            raise ConfigError(f"AnvioExecutor :: unknown backend '{self.backend}'. "
                              f"Valid choices are: {', '.join(BACKENDS)}.")

        self.entry = registry.get_tool(tool_name)
        if self.entry is None and self.backend != 'native' and anvio.DEBUG:
            self.run.warning(f"'{tool_name}' is not in the anvi'o registry. "
                             f"Container-based backends will not be available for this tool.")


    def run_command(self, cmd, raise_on_error=True):
        """Run cmd via the configured backend and return a RunResult.

        Parameters
        ==========
        cmd : list or str
            Command and arguments to execute.
        raise_on_error : bool
            If True (default), raise CommandError when the command exits
            with a non-zero return code.

        Returns
        =======
        RunResult
            Named tuple with fields: stdout, stderr, returncode, cmd.
            ``cmd`` is the fully-resolved argv list that was actually
            executed, including any container wrapper arguments.

        Raises
        ======
        ConfigError
            If the backend is misconfigured or a required helper is absent.
        CommandError
            If raise_on_error is True and the command exits non-zero.
        """
        if isinstance(cmd, str):
            cmd = cmd.split()

        dispatch = {
            'native':      self._run_native,
            'docker':      self._run_docker,
            'singularity': self._run_singularity,
            'conda':       self._run_conda,
        }

        return dispatch[self.backend](list(cmd), raise_on_error)


    def get_conda_env_name(self):
        """Return the expected conda environment name for this tool.

        Returns None if the tool has no conda runner entry.

        Returns
        =======
        str or None
        """
        spec = registry.get_runner(self.tool_name, 'conda')
        return self._conda_env_name(spec) if spec else None


    def _execute(self, final_cmd):
        """Run final_cmd via subprocess and return (stdout, stderr, returncode)."""
        if anvio.DEBUG:
            self.run.info(f'[{self.backend}] cmd', ' '.join(str(x) for x in final_cmd))

        try:
            proc = subprocess.run(
                final_cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
            )
        except OSError as e:
            raise ConfigError(f"AnvioExecutor :: failed to launch '{final_cmd[0]}': {e}")

        return proc.stdout.decode(), proc.stderr.decode(), proc.returncode


    def _finish(self, final_cmd, stdout, stderr, returncode, raise_on_error):
        """Build a RunResult, write the log if requested, and optionally raise on non-zero exit."""
        result = RunResult(stdout=stdout, stderr=stderr,
                           returncode=returncode, cmd=final_cmd)

        if self.log:
            self._write_log(result)

        if raise_on_error and returncode != 0:
            msg = (f"'{self.tool_name}' exited with code {returncode} "
                   f"(backend: {self.backend}, "
                   f"cmd: {' '.join(str(x) for x in final_cmd)})")
            if stderr.strip():
                msg += f"\n{stderr.strip()}"
            raise CommandError(msg)

        return result


    def _write_log(self, result):
        """Write a structured log entry for the completed command to self.log."""
        with open(self.log, 'w') as f:
            f.write(f'# DATE: {datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")}\n')
            f.write(f'# CMD: {" ".join(str(x) for x in result.cmd)}\n')
            f.write(f'# BACKEND: {self.backend}\n')
            f.write(f'# RETURN CODE: {result.returncode}\n')
            if result.stdout:
                f.write(f'\n{result.stdout}')
            if result.stderr:
                f.write(f'\n# STDERR:\n{result.stderr}')


    def _run_native(self, cmd, raise_on_error):
        stdout, stderr, returncode = self._execute(cmd)
        return self._finish(cmd, stdout, stderr, returncode, raise_on_error)


    def _run_docker(self, cmd, raise_on_error):
        image = registry.get_runner(self.tool_name, 'docker')
        if not image:
            raise ConfigError(f"AnvioExecutor :: no Docker image registered for '{self.tool_name}'. "
                              f"Add a 'docker' runner entry to anvio/registry/third-party-programs.yaml.")

        if 'HASH_PLACEHOLDER' in image:
            raise ConfigError(f"AnvioExecutor :: the Docker image for '{self.tool_name}' still contains "
                              f"a HASH_PLACEHOLDER ({image}). Fill in the correct biocontainers build "
                              f"hash in anvio/registry/third-party-programs.yaml before using this backend.")

        if not shutil.which('docker'):
            raise ConfigError("AnvioExecutor :: 'docker' is not in PATH. "
                              "Please install Docker and make sure the daemon is running.")

        # Mount work_dir at the same absolute path so file arguments in cmd
        # resolve identically inside and outside the container (Nextflow-style).
        final_cmd = [
            'docker', 'run', '--rm',
            '-v', f'{self.work_dir}:{self.work_dir}',
            '-w', self.work_dir,
            '-u', f'{os.getuid()}:{os.getgid()}',
            image,
        ] + cmd

        stdout, stderr, returncode = self._execute(final_cmd)
        return self._finish(final_cmd, stdout, stderr, returncode, raise_on_error)


    def _run_singularity(self, cmd, raise_on_error):
        image = registry.get_runner(self.tool_name, 'singularity')
        if not image:
            raise ConfigError(f"AnvioExecutor :: no Singularity image registered for '{self.tool_name}'. "
                              f"Add a 'singularity' runner entry to anvio/registry/third-party-programs.yaml.")

        if 'HASH_PLACEHOLDER' in image:
            raise ConfigError(f"AnvioExecutor :: the Singularity image for '{self.tool_name}' still contains "
                              f"a HASH_PLACEHOLDER ({image}). Fill in the correct biocontainers build "
                              f"hash in anvio/registry/third-party-programs.yaml before using this backend.")

        if not shutil.which('singularity'):
            raise ConfigError("AnvioExecutor :: 'singularity' is not in PATH.")

        final_cmd = [
            'singularity', 'exec',
            '--bind', self.work_dir,
            '--pwd', self.work_dir,
            image,
        ] + cmd

        stdout, stderr, returncode = self._execute(final_cmd)
        return self._finish(final_cmd, stdout, stderr, returncode, raise_on_error)


    def _run_conda(self, cmd, raise_on_error):
        spec = registry.get_runner(self.tool_name, 'conda')
        if not spec:
            raise ConfigError(f"AnvioExecutor :: no conda spec registered for '{self.tool_name}'. "
                              f"Add a 'conda' runner entry to anvio/registry/third-party-programs.yaml.")

        # Prefer mamba/micromamba for speed; fall back to conda.
        conda_exe = (shutil.which('mamba')
                     or shutil.which('micromamba')
                     or shutil.which('conda'))
        if not conda_exe:
            raise ConfigError("AnvioExecutor :: no conda/mamba/micromamba found in PATH. "
                              "Please install one of them first.")

        env_name = self._conda_env_name(spec)

        # Verify the env exists rather than creating one silently — env creation
        # is slow and the user should do it deliberately.
        try:
            check = subprocess.run(
                [conda_exe, 'env', 'list', '--json'],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
            )
            envs = json.loads(check.stdout).get('envs', [])
        except Exception:
            envs = []

        if not any(os.path.basename(e) == env_name for e in envs):
            pkg_spec = spec.split('::')[-1] if '::' in spec else spec
            raise ConfigError(f"AnvioExecutor :: conda environment '{env_name}' does not exist. "
                              f"Create it first and then re-run:\n\n"
                              f"    conda create -n {env_name} -c bioconda -c conda-forge {pkg_spec}\n")

        final_cmd = [conda_exe, 'run', '--no-capture-output', '--name', env_name, '--'] + cmd
        stdout, stderr, returncode = self._execute(final_cmd)
        return self._finish(final_cmd, stdout, stderr, returncode, raise_on_error)


    def _conda_env_name(self, spec):
        """Derive a canonical conda env name from a registry package spec.

        Examples
        ========
        'bioconda::diamond=2.1.9'   → 'anvio-tool-diamond-2.1.9'
        'conda-forge::r-base=4.3.3' → 'anvio-tool-r-base-4.3.3'
        """
        pkg = spec.split('::')[-1] if '::' in spec else spec
        return f'anvio-tool-{pkg.replace("=", "-")}'

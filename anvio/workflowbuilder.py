"""Serve the local anvi'o workflow builder interface."""

import contextlib
import importlib
import importlib.metadata
import importlib.util
import io
import json
import pkgutil
import tempfile

from argparse import Namespace
from datetime import datetime, timezone
from functools import partial
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.parse import urlsplit

import anvio
import anvio.terminal as terminal
import anvio.utils as utils

from anvio.errors import ConfigError


@contextlib.contextmanager
def quiet_anvio_imports():
    """Silence imports while collecting metadata from local CLI modules."""
    with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
        yield


def local_console_scripts():
    """Return local anvi'o console scripts and their import targets."""
    scripts = {}
    entry_points = importlib.metadata.entry_points()
    selected = entry_points.select(group='console_scripts') if hasattr(entry_points, 'select') else entry_points.get('console_scripts', [])
    for entry in selected:
        if entry.name.startswith('anvi-') and entry.value.startswith('anvio.cli.'):
            scripts[entry.name] = entry.value.split(':', 1)[0]

    import anvio.cli as cli

    for module_info in pkgutil.iter_modules(cli.__path__):
        module_name = f'anvio.cli.{module_info.name}'
        scripts.setdefault(f"anvi-{module_info.name.replace('_', '-')}", module_name)

    return dict(sorted(scripts.items()))


def local_program_metadata(command, module_name):
    """Return graph metadata declared by one local anvi'o CLI module."""
    try:
        with quiet_anvio_imports():
            module = importlib.import_module(module_name)
    except Exception as error:
        return {
            'id': command,
            'kind': 'program',
            'title': command,
            'command': command,
            'description': f'Local metadata import failed: {error}',
            'requires': [],
            'provides': [],
            'parameters': [],
        }

    return {
        'id': command,
        'kind': 'program',
        'title': command,
        'command': command,
        'description': getattr(module, '__description__', ''),
        'requires': list(getattr(module, '__requires__', [])),
        'provides': list(getattr(module, '__provides__', [])),
        'parameters': [],
    }


def installed_console_script_paths():
    """Return installed anvi'o console scripts and their module paths."""
    scripts = {}
    entry_points = importlib.metadata.entry_points()
    selected = entry_points.select(group='console_scripts') if hasattr(entry_points, 'select') else entry_points.get('console_scripts', [])
    for entry in selected:
        if not entry.name.startswith('anvi-') or not entry.value.startswith('anvio.cli.'):
            continue

        module_name = entry.value.split(':', 1)[0]
        module_spec = importlib.util.find_spec(module_name)
        if module_spec and module_spec.origin:
            scripts[entry.name] = module_spec.origin

    if not scripts:
        raise ConfigError("The workflow builder could not find installed anvi'o console scripts. "
                          "Please reinstall this package and try again.")

    return scripts


def local_program_network():
    """Use anvi'o's ProgramsNetwork implementation to build the graph."""
    try:
        import anvio.programs as programs

        class InstalledProgramsNetwork(programs.ProgramsNetwork):
            """Build a programs network from installed entry points."""

            def get_anvio_program_names_and_their_paths(self):
                return installed_console_script_paths()

        output_path = Path(tempfile.gettempdir()) / f"anvio-programs-network-{datetime.now(timezone.utc).timestamp()}.json"
        args = Namespace(output_file=str(output_path), program_names_to_focus=None)
        source_manifest = Path(anvio.__file__).resolve().parent.parent / 'pyproject.toml'
        network_class = programs.ProgramsNetwork if source_manifest.exists() else InstalledProgramsNetwork
        with quiet_anvio_imports():
            network_class(args, r=terminal.Run(verbose=False), p=terminal.Progress()).generate()
        try:
            return json.loads(output_path.read_text(encoding='utf-8'))
        finally:
            output_path.unlink(missing_ok=True)
    except Exception as error:
        raise ConfigError("Anvi'o could not build the local programs network required by the workflow builder. "
                          f"The underlying error was: {error}")


def local_artifacts(programs):
    """Return artifacts declared by local CLI modules."""
    artifacts = {}
    for program in programs:
        for artifact_id in program.get('requires', []) + program.get('provides', []):
            artifacts.setdefault(artifact_id, {
                'id': artifact_id,
                'kind': 'data',
                'title': artifact_id,
                'artifactType': 'LOCAL',
                'description': "Artifact declared by the local anvi'o installation.",
                'requires': [],
                'provides': [],
                'parameters': [],
            })
    return sorted(artifacts.values(), key=lambda artifact: artifact['id'])


def local_runtime():
    """Return browser metadata for the anvi'o package serving this interface."""
    programs = [local_program_metadata(command, module_name) for command, module_name in local_console_scripts().items()]
    return {
        'schemaVersion': '0.1.0',
        'generatedAt': datetime.now(timezone.utc).isoformat(),
        'repository': "local anvi'o installation",
        'version': {
            'id': 'local',
            'label': f"local anvi'o {anvio.__version__}",
            'version': anvio.__version__,
            'packagePath': str(Path(anvio.__file__).resolve()),
        },
        'programs': programs,
        'artifacts': local_artifacts(programs),
        'programNetwork': local_program_network(),
        'classes': [],
        'parameterRelations': {
            'label': f"local anvi'o {anvio.__version__}",
            'resolvedRef': 'local-package',
            'programs': {},
        },
    }


class WorkflowBuilderRequestHandler(SimpleHTTPRequestHandler):
    """Serve static workflow builder files and generated local metadata."""

    runtime_javascript = ''


    def log_message(self, format, *args):
        """Keep the terminal focused on workflow-builder lifecycle messages."""
        pass


    def do_GET(self):
        if urlsplit(self.path).path == '/data/anvio-local-runtime.js':
            payload = self.runtime_javascript.encode('utf-8')
            self.send_response(200)
            self.send_header('Content-Type', 'text/javascript; charset=utf-8')
            self.send_header('Content-Length', str(len(payload)))
            self.end_headers()
            self.wfile.write(payload)
            return
        super().do_GET()


class WorkflowBuilder:
    """Launch the local workflow builder web interface."""

    def __init__(self, args, run=terminal.Run()):
        self.args = args
        self.run = run
        self.static_dir = Path(anvio.__file__).resolve().parent / 'data' / 'workflow_builder'

        if not self.static_dir.is_dir():
            raise ConfigError(f"The workflow builder interface directory is missing at '{self.static_dir}'. "
                              f"This anvi'o installation is incomplete. Please reinstall it and try again.")


    def serve(self):
        """Serve the interface until the user interrupts the process."""
        runtime_javascript = f"window.ANVIO_LOCAL_RUNTIME = {json.dumps(local_runtime(), indent=2, sort_keys=True)};\n"
        WorkflowBuilderRequestHandler.runtime_javascript = runtime_javascript
        handler = partial(WorkflowBuilderRequestHandler, directory=str(self.static_dir))
        server = ThreadingHTTPServer((self.args.ip_address, self.args.port_number), handler)
        url = f"http://{self.args.ip_address}:{self.args.port_number}/index.html"

        self.run.info('Workflow builder URL', url)
        self.run.info_single('Press Ctrl+C to stop the workflow builder server.', nl_before=1)
        if not self.args.server_only:
            utils.open_url_in_browser(url, browser_path=self.args.browser_path, run=self.run)

        try:
            server.serve_forever()
        except KeyboardInterrupt:
            self.run.info_single('Workflow builder server stopped.', nl_before=1)
        finally:
            server.server_close()

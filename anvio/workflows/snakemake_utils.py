"""Small utility functions used by anvi'o Snakefiles."""

import os
import re

import anvio.terminal as terminal


# Path to the long-read technology preset map (token -> per-tool presets). See the
# file's own header for the format. Loaded lazily and cached by get_lr_technology_presets().
LR_TECHNOLOGY_PRESETS_YAML = os.path.join(os.path.dirname(__file__), 'lr_technology_presets.yaml')

# Cache so we parse the YAML at most once per process.
_lr_technology_presets_cache = None

# Remembers which (tool, executable) pairs we have already version-probed this process, so the
# advisory version check shells out to '<tool> --version' at most once per process.
_version_probe_done = set()

# Flye's mutually exclusive read-type flags. Exactly one applies to a given long-read assembly;
# used both to validate a config-provided flag and to resolve one from an lr_technology token.
FLYE_READ_TYPE_FLAGS = ["--pacbio-raw", "--pacbio-corr", "--pacbio-hifi",
                        "--nano-raw", "--nano-corr", "--nano-hq"]


def D(debug_message, debug_log_file_path=".SNAKEMAKEDEBUG"):
    """Append a timestamped debug message to a Snakemake debug log."""
    with open(debug_log_file_path, 'a') as output:
            output.write(terminal.get_date() + '\n')
            output.write(str(debug_message) + '\n\n')


def regex_from_ids(ids):
    """Return a wildcard constraint regex from a list of IDs."""
    return r"(?:{})".format("|".join(re.escape(str(i)) for i in ids)) if ids else r"DO_NOT_MATCH"


# Directory holding the per-rule conda environment files anvi'o ships (one <tool>.yaml per
# rule that wraps a third-party program). See conda_envs/README.md for the full list.
CONDA_ENVS_DIR = os.path.join(os.path.dirname(__file__), 'conda_envs')


def get_anvio_conda_yaml_path(tool):
    """Return the absolute path to the conda YAML anvi'o ships for a tool, or None if none exists.

    The file name matches the rule/config key (e.g. 'nanoplot' -> conda_envs/nanoplot.yaml). This
    is resolved at runtime from the installed anvi'o location, so a config using
    'use_anvio_conda_yaml' stays reproducible across machines (no hard-coded repo path).
    """
    path = os.path.join(CONDA_ENVS_DIR, f'{tool}.yaml')
    return path if os.path.exists(path) else None


def get_conda_yaml_path(workflow, tool):
    """Return the absolute conda YAML path for a workflow tool, or None.

    Precedence (mutual exclusivity is enforced by the workflow sanity checks):
      1. 'use_anvio_conda_yaml': True  -> the YAML anvi'o ships for this tool (see above);
      2. 'conda_yaml': <path>          -> that explicit, user-provided path;
      3. otherwise                     -> None (tool comes from $PATH or an existing conda_env).
    """
    if workflow.get_param_value_from_config([tool, 'use_anvio_conda_yaml']) == True:
        # sanity-checked up front, but stay defensive here since this feeds the conda: directive
        path = get_anvio_conda_yaml_path(tool)
        if not path:
            from anvio.errors import ConfigError
            raise ConfigError(f"'{tool}' is set to use the anvi'o-shipped conda env file, but anvi'o "
                              f"does not ship one for this rule (expected '{tool}.yaml' in {CONDA_ENVS_DIR}).")
        return path

    path = workflow.get_param_value_from_config([tool, 'conda_yaml'])
    if not path:
        return None

    return os.path.abspath(os.path.expanduser(path))


def get_conda_env_prefix(workflow, tool):
    """Return the command prefix for a configured existing conda environment."""
    name = workflow.get_param_value_from_config([tool, 'conda_env'])
    return f"conda run -n {name}" if name else ""


def gunzip_file(file_path, log, shell, output_path=None):
    """Uncompress a gzipped file and return the uncompressed path."""
    uncompressed_file_path = output_path or os.path.splitext(file_path)[0]
    shell("gunzip < %s > %s 2>> %s" % (file_path, uncompressed_file_path, log))
    return uncompressed_file_path


def get_lr_technology_presets():
    """Load and cache the long-read technology preset map (lr_technology_presets.yaml).

    Returns the parsed dict with two top-level keys: 'tools' (per-tool tested_versions)
    and 'technologies' (token -> {minimap2, flye} presets). Parsed once per process.
    """
    global _lr_technology_presets_cache

    if _lr_technology_presets_cache is None:
        # imported here rather than at module top to keep the Snakefile import surface light
        import anvio.utils as u

        presets = u.get_yaml_as_dict(LR_TECHNOLOGY_PRESETS_YAML) or {}
        # normalize to always have the two sections so callers can index without guarding
        presets.setdefault('tools', {})
        presets.setdefault('technologies', {})

        # Completeness check (dev-facing): every technology token must define a preset for every
        # tool declared under `tools:`. Token vocabulary is validated against this map elsewhere,
        # but that only checks membership — without this, a token missing one tool's preset would
        # silently fall through to the config-flag path (wrong preset, no error). This map is
        # hand-edited by developers, so we fail loudly here the moment it goes inconsistent.
        expected_tools = set(presets['tools'].keys())
        gaps = {tech: sorted(expected_tools - {t for t, v in (tool_map or {}).items() if v})
                for tech, tool_map in presets['technologies'].items()}
        gaps = {tech: missing for tech, missing in gaps.items() if missing}
        if gaps:
            from anvio.errors import ConfigError
            details = '; '.join(f"'{tech}' is missing: {', '.join(missing)}" for tech, missing in sorted(gaps.items()))
            raise ConfigError(f"The long-read technology preset map ({LR_TECHNOLOGY_PRESETS_YAML}) is "
                              f"incomplete: every technology token must define a preset for each tool "
                              f"listed under 'tools:' ({', '.join(sorted(expected_tools))}). {details}. "
                              f"Please add the missing preset(s) to that file.")

        _lr_technology_presets_cache = presets

    return _lr_technology_presets_cache


def get_lr_technology_map():
    """Return just the token -> per-tool preset mapping from the preset file."""
    return get_lr_technology_presets()['technologies']


def get_valid_lr_technologies():
    """Return the set of accepted lr_technology tokens (the keys of the technology map)."""
    return set(get_lr_technology_map().keys())


def get_lr_preset(technology, tool):
    """Return the preset/flag for a (technology token, tool) pair, or None if not defined.

    `tool` is one of 'minimap2', 'flye'. Returns None when the technology is
    unknown or the tool has no preset for it (callers validate/raise as appropriate).
    """
    return get_lr_technology_map().get(technology, {}).get(tool)


def warn_if_tool_version_untested(tool, executable=None, run=None):
    """Emit a soft WARNING (never an error) if an installed tool's version is untested.

    Best-effort and non-fatal: compares the version reported by `<executable> --version`
    against the `tested_versions` recorded for `tool` in lr_technology_presets.yaml. Any
    failure to locate the program or parse its version is silently ignored — this check
    must never block a run. Mirrors the 'tested_versions' pattern in anvio/drivers/trnscan_se.py.
    """
    import anvio.utils as u

    executable = executable or tool
    run = run or terminal.Run()

    # probe each (tool, executable) at most once per process (belt-and-suspenders: this is already
    # invoked only on real runs, but guarantees no repeated subprocess if called more than once)
    probe_key = (tool, executable)
    if probe_key in _version_probe_done:
        return
    _version_probe_done.add(probe_key)

    tested_versions = get_lr_technology_presets()['tools'].get(tool, {}).get('tested_versions', [])
    if not tested_versions:
        return

    try:
        if not u.is_program_exists(executable, dont_raise=True):
            return

        output, ret_code = u.get_command_output_from_shell('%s --version' % executable)
        if isinstance(output, bytes):
            output = output.decode('utf-8', errors='replace')

        # pull the first version-looking token (e.g. '2.28', '2.9.5', '1.2.0c') out of the output
        match = re.search(r'\d+\.\d+(?:\.\d+)?[a-z]?', output)
        if not match:
            return
        installed_version = match.group(0)

        if installed_version not in tested_versions:
            run.warning(
                f"You have {executable} version '{installed_version}' installed, but anvi'o's "
                f"long-read presets for '{tool}' were only validated against these version(s): "
                f"{', '.join(tested_versions)}. This is just a heads up (not an error): the presets "
                f"may still work perfectly, but if you run into trouble with long-read "
                f"{tool} steps, a version mismatch is a good first thing to check.",
                header=f"UNTESTED {tool.upper()} VERSION", lc="yellow",
            )
    except Exception:
        # version probing is strictly best-effort; never let it break a workflow
        return

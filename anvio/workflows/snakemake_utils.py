"""Small utility functions used by anvi'o Snakefiles."""

import os
import re

import anvio.terminal as terminal


# Path to the long-read technology preset map (token -> per-tool presets). See the
# file's own header for the format. Loaded lazily and cached by get_lr_technology_presets().
LR_TECHNOLOGY_PRESETS_YAML = os.path.join(os.path.dirname(__file__), 'lr_technology_presets.yaml')

# Cache so we parse the YAML at most once per process.
_lr_technology_presets_cache = None


def D(debug_message, debug_log_file_path=".SNAKEMAKEDEBUG"):
    """Append a timestamped debug message to a Snakemake debug log."""
    with open(debug_log_file_path, 'a') as output:
            output.write(terminal.get_date() + '\n')
            output.write(str(debug_message) + '\n\n')


def regex_from_ids(ids):
    """Return a wildcard constraint regex from a list of IDs."""
    return r"(?:{})".format("|".join(re.escape(str(i)) for i in ids)) if ids else r"DO_NOT_MATCH"


def get_conda_yaml_path(workflow, tool):
    """Return the absolute conda YAML path configured for a workflow tool."""
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
    and 'technologies' (token -> {longqc, minimap2, flye} presets). Parsed once per process.
    """
    global _lr_technology_presets_cache

    if _lr_technology_presets_cache is None:
        # imported here rather than at module top to keep the Snakefile import surface light
        import anvio.utils as u

        presets = u.get_yaml_as_dict(LR_TECHNOLOGY_PRESETS_YAML) or {}
        # normalize to always have the two sections so callers can index without guarding
        presets.setdefault('tools', {})
        presets.setdefault('technologies', {})
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

    `tool` is one of 'longqc', 'minimap2', 'flye'. Returns None when the technology is
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

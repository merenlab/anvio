# anvi'o Workflow Builder

The workflow builder is distributed with this anvi'o fork. It runs locally and does not depend on the location of the cloned repository after installation.

## Install

Activate a Python 3.10 anvi'o environment, then install the package from the repository root:

```bash
pip install .
```

## Launch

```bash
anvi-build-workflow
```

The command starts a local HTTP server and opens the interface in your default browser. Keep the terminal open while using the application. Press `Ctrl+C` to stop the server.

To start the server without opening a browser:

```bash
anvi-build-workflow --server-only
```

The terminal prints the local URL to open manually.

## Development Mode

For development, use an editable installation:

```bash
pip install --config-settings editable_mode=compat -e .
```

Editable mode intentionally follows the repository path. Use the standard installation above when you need a standalone command that remains valid after moving the clone.

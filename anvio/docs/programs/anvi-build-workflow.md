Launch a local web interface to design anvi'o workflows.

The workflow builder reads the programs and artifacts available in your local anvi'o installation, then opens an interactive canvas in your browser. You can search for anvi'o programs, place blocks, connect compatible inputs and outputs, edit parameters, validate the resulting workflow, and export it as a Bash script, a Snakemake file, or an image.

The interface is served locally. It does not send your workflow or your parameter values to a remote service.

### Standalone installation

Install this package from any clone or source archive:

{{ codestart }}
pip install .
{{ codestop }}

The installation copies the workflow builder and its metadata into the active Python environment. After installation, `anvi-build-workflow` does not depend on the location of the source directory. You can move or remove the clone without breaking the command.

### Launching the workflow builder

Run the program without parameters:

{{ codestart }}
anvi-build-workflow
{{ codestop }}

The command finds an available port, starts a local HTTP server, and opens the workflow builder in your default browser. Keep the command running while you use the interface. Press `Ctrl+C` when you are finished.

### Starting the server without opening a browser

For a remote session, or when you prefer to open the URL yourself, use `--server-only`:

{{ codestart }}
anvi-build-workflow --server-only
{{ codestop }}

The program prints the URL of the workflow builder. You can also choose the listening address and port:

{{ codestart }}
anvi-build-workflow --server-only \
                    --ip-address localhost \
                    --port-number 8080
{{ codestop }}

### Choosing another browser

To open the interface with a specific browser executable, use `--browser-path`:

{{ codestart }}
anvi-build-workflow --browser-path /path/to/browser
{{ codestop }}

### Notes

The graph is built from the local anvi'o package through anvi'o's own program metadata and program-network implementation. The interface combines that local graph with its embedded workflow-builder metadata. Documentation links target the official help pages at [anvio.org](https://anvio.org/help/).

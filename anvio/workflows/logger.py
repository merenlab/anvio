"""Register the bundled Snakemake logger for ordinary editable installs."""

from anvio.errors import ConfigError


def register_workflow_logger():
    # Import lazily to avoid the plugin -> workflow manifest -> workflows cycle,
    # and keep dependency-list and dry-run paths independent of registration.
    try:
        from snakemake_interface_common.exceptions import InvalidPluginException
        from snakemake_interface_logger_plugins.registry import LoggerPluginRegistry
    except ImportError as e:
        raise ConfigError("Could not import the Snakemake interfaces needed to register the anvi'o workflow logger. "
                          "Please check your anvi'o installation.") from e

    try:
        registry = LoggerPluginRegistry()
        if registry.is_installed('anvio'):
            return

        import snakemake_logger_plugin_anvio

        registry.register_plugin('snakemake_logger_plugin_anvio', snakemake_logger_plugin_anvio)
    # Registry validation uses issubclass, which raises TypeError for non-class handlers.
    except (ImportError, InvalidPluginException, TypeError) as e:
        raise ConfigError("Could not register the anvi'o workflow logger with Snakemake. "
                          "Please check your anvi'o and Snakemake installations.") from e

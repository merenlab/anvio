"""Snakemake 9 logger plugin that records anvi'o workflow job status."""

from logging import LogRecord

from snakemake_interface_logger_plugins.base import LogHandlerBase
from snakemake_interface_logger_plugins.common import LogEvent

from anvio.workflows.scripts.snakemake_log_handler import log_handler


class LogHandler(LogHandlerBase):
    """Adapt Snakemake 9 log records to anvi'o's workflow manifest handler."""

    @property
    def writes_to_stream(self):
        return False

    @property
    def writes_to_file(self):
        return False

    @property
    def has_filter(self):
        return True

    @property
    def has_formatter(self):
        return False

    @property
    def needs_rulegraph(self):
        return False

    def __post_init__(self):
        self.addFilter(self._manifest_event_filter)

    @staticmethod
    def _manifest_event_filter(record):
        event = getattr(record, 'event', None)
        if event is None:
            return 'Complete log' in record.getMessage()

        return event in {
            LogEvent.JOB_INFO,
            LogEvent.JOB_FINISHED,
            LogEvent.JOB_ERROR,
            LogEvent.ERROR,
        }

    def emit(self, record: LogRecord):
        event = getattr(record, 'event', None)
        formatted_message = record.getMessage()
        if event is None and 'Complete log' not in formatted_message:
            return

        message = record.__dict__.copy()
        # ``LogRecord.name`` is the Python logger name (usually "snakemake"),
        # while the legacy callback used the same key for a workflow rule name.
        message.pop('name', None)
        message['level'] = event.value if isinstance(event, LogEvent) else str(event)

        if formatted_message:
            message.setdefault('msg', formatted_message)

        log_handler(message)

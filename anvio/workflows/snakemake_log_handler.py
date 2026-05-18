"""Snakemake log handler for anvi'o workflow manifests."""

import os
import re

from anvio.workflows.manifest import append_manifest_row, update_snakemake_log_path


JOBS = {}
SNAKEMAKE_LOG_PATH = ''


def _as_dict(value):
    if value is None:
        return {}

    if isinstance(value, dict):
        return value

    if isinstance(value, str):
        wildcards = {}
        for part in value.replace(',', ' ').split():
            if '=' in part:
                key, val = part.split('=', 1)
                wildcards[key.strip()] = val.strip()
        return wildcards

    if hasattr(value, '_asdict'):
        return value._asdict()

    if hasattr(value, 'items'):
        return {key: val for key, val in value.items()}

    try:
        return dict(value)
    except (TypeError, ValueError):
        return {}


def _as_list(value):
    if value is None:
        return []

    if isinstance(value, str):
        return [value]

    try:
        return list(value)
    except TypeError:
        return [value]


def _job_id(message):
    for key in ('jobid', 'job_id', 'job'):
        if key in message:
            return str(message[key])

    return None


def _rule_name(message):
    return message.get('name') or message.get('rule') or message.get('rule_name') or ''


def _wildcards(message):
    return _as_dict(message.get('wildcards') or message.get('wildcard_values'))


def _log_path(message):
    logs = _as_list(message.get('log') or message.get('logs') or message.get('logfile') or message.get('logfiles'))
    return ','.join(str(log) for log in logs if log)


def _message_strings(message):
    strings = []

    for value in message.values():
        if isinstance(value, str):
            strings.append(value)
        elif isinstance(value, (list, tuple)):
            strings.extend(str(item) for item in value if isinstance(item, str))

    return strings


def _snakemake_log_path(message):
    for text in _message_strings(message):
        match = re.search(r'Complete log:\s*(\S+)', text)
        if match:
            return match.group(1)

    return ''


def _remember_snakemake_log_path(message):
    global SNAKEMAKE_LOG_PATH

    snakemake_log_path = _snakemake_log_path(message)

    if not snakemake_log_path:
        return

    SNAKEMAKE_LOG_PATH = snakemake_log_path

    manifest_path = os.environ.get('ANVIO_WORKFLOW_MANIFEST_PATH')
    if manifest_path:
        update_snakemake_log_path(manifest_path, SNAKEMAKE_LOG_PATH)


def _entry_from_message(message):
    wildcards = _wildcards(message)

    return {'rule': _rule_name(message),
            'group': wildcards.get('group', ''),
            'read': wildcards.get('readset', ''),
            'log_path': _log_path(message)}


def _remember_job(message):
    job_id = _job_id(message)

    if not job_id:
        return

    JOBS[job_id] = _entry_from_message(message)


def _record_job(message, status):
    manifest_path = os.environ.get('ANVIO_WORKFLOW_MANIFEST_PATH')

    if not manifest_path:
        return

    job_id = _job_id(message)
    entry = JOBS.get(job_id, {}) if job_id else {}
    current_entry = _entry_from_message(message)

    for key, value in current_entry.items():
        if value:
            entry[key] = value

    append_manifest_row(manifest_path,
                        status,
                        entry.get('rule', ''),
                        group=entry.get('group', ''),
                        read=entry.get('read', ''),
                        log_path=entry.get('log_path', ''),
                        snakemake_log_path=SNAKEMAKE_LOG_PATH)


def log_handler(message):
    message = _as_dict(message)
    level = message.get('level')

    _remember_snakemake_log_path(message)

    if level == 'job_info':
        _remember_job(message)
    elif level == 'job_finished':
        _record_job(message, 'succeeded')
    elif level == 'job_error' or (level == 'error' and (_job_id(message) or _rule_name(message))):
        _record_job(message, 'failed')

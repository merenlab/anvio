"""Helpers for workflow status manifests."""

import os
import fcntl


MANIFEST_HEADER = ['status', 'rule', 'group', 'read', 'log_path', 'snakemake_log_path']


def get_workflow_logs_dir(logs_dir, workflow_name):
    logs_parent_dir = os.path.dirname(logs_dir)
    logs_dir_name = os.path.basename(logs_dir)

    if workflow_name in ('pangenomics', 'phylogenomics'):
        workflow_logs_dir_name = workflow_name
    elif logs_dir_name == '00_LOGS':
        workflow_logs_dir_name = workflow_name
    elif logs_dir_name.startswith('00_LOGS-') or logs_dir_name.startswith('00_LOGS_'):
        workflow_logs_dir_name = logs_dir_name[len('00_LOGS') + 1:]
    else:
        workflow_logs_dir_name = workflow_name

    return os.path.join(logs_parent_dir, '00_LOGS', workflow_logs_dir_name)


def initialize_manifest(manifest_path):
    os.makedirs(os.path.dirname(manifest_path), exist_ok=True)

    with open(manifest_path, 'w') as output:
        output.write('\t'.join(MANIFEST_HEADER) + '\n')


def append_manifest_row(manifest_path, status, rule, group='', read='', log_path='', snakemake_log_path=''):
    os.makedirs(os.path.dirname(manifest_path), exist_ok=True)

    row = [status, rule or '', group or '', read or '', log_path or '', snakemake_log_path or '']
    row = [str(value).replace('\t', ' ').replace('\n', ' ') for value in row]

    with open(manifest_path, 'a') as output:
        fcntl.flock(output.fileno(), fcntl.LOCK_EX)
        try:
            output.write('\t'.join(row) + '\n')
            output.flush()
        finally:
            fcntl.flock(output.fileno(), fcntl.LOCK_UN)


def update_snakemake_log_path(manifest_path, snakemake_log_path):
    if not snakemake_log_path or not os.path.exists(manifest_path):
        return

    with open(manifest_path, 'r+') as manifest:
        fcntl.flock(manifest.fileno(), fcntl.LOCK_EX)
        try:
            lines = manifest.readlines()

            if not lines:
                lines = ['\t'.join(MANIFEST_HEADER) + '\n']

            header = lines[0].rstrip('\n').split('\t')
            if header != MANIFEST_HEADER:
                header = MANIFEST_HEADER

            updated_lines = ['\t'.join(header) + '\n']
            for line in lines[1:]:
                row = line.rstrip('\n').split('\t')
                while len(row) < len(MANIFEST_HEADER):
                    row.append('')

                row[5] = snakemake_log_path
                updated_lines.append('\t'.join(row[:len(MANIFEST_HEADER)]) + '\n')

            manifest.seek(0)
            manifest.truncate()
            manifest.writelines(updated_lines)
            manifest.flush()
        finally:
            fcntl.flock(manifest.fileno(), fcntl.LOCK_UN)

"""Helpers for workflow status manifests."""

import os
import fcntl


MANIFEST_HEADER = ['status', 'rule', 'group', 'read', 'log_path']


def initialize_manifest(manifest_path):
    os.makedirs(os.path.dirname(manifest_path), exist_ok=True)

    with open(manifest_path, 'w') as output:
        output.write('\t'.join(MANIFEST_HEADER) + '\n')


def append_manifest_row(manifest_path, status, rule, group='', read='', log_path=''):
    os.makedirs(os.path.dirname(manifest_path), exist_ok=True)

    row = [status, rule or '', group or '', read or '', log_path or '']
    row = [str(value).replace('\t', ' ').replace('\n', ' ') for value in row]

    with open(manifest_path, 'a') as output:
        fcntl.flock(output.fileno(), fcntl.LOCK_EX)
        try:
            output.write('\t'.join(row) + '\n')
            output.flush()
        finally:
            fcntl.flock(output.fileno(), fcntl.LOCK_UN)

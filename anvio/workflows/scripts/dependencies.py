"""Helpers for workflow dependency checks."""

import re
import shlex


SHELL_TOKENS_TO_SKIP = {'',
                        '!',
                        '(',
                        ')',
                        '[',
                        '[[',
                        ']',
                        ']]',
                        'case',
                        'do',
                        'done',
                        'elif',
                        'else',
                        'esac',
                        'fi',
                        'for',
                        'function',
                        'if',
                        'in',
                        'then',
                        'while'}

SHELL_COMMANDS_TO_IGNORE = {'cat',
                            'cd',
                            'cp',
                            'echo',
                            'exit',
                            'gzip',
                            'mkdir',
                            'mv',
                            'printf',
                            'rm',
                            'sed',
                            'test',
                            'touch'}

SHELL_CONTROL_OPERATORS = ('&&', '||', ';')


def _has_open_quote(line):
    in_single_quote = False
    in_double_quote = False
    escaped = False

    for char in line:
        if escaped:
            escaped = False
            continue

        if char == '\\' and not in_single_quote:
            escaped = True
            continue

        if char == "'" and not in_double_quote:
            in_single_quote = not in_single_quote
        elif char == '"' and not in_single_quote:
            in_double_quote = not in_double_quote

    return in_single_quote or in_double_quote


def _logical_shell_lines(shellcmd):
    lines = []
    pending = ''

    for raw_line in shellcmd.splitlines():
        line = raw_line.rstrip()

        if pending:
            line = pending + ' ' + line.strip()

        if line.endswith('\\'):
            pending = line[:-1].rstrip()
            continue

        if _has_open_quote(line):
            pending = line
            continue

        lines.append(line)
        pending = ''

    if pending:
        lines.append(pending)

    return lines


def _strip_inline_comment(line):
    in_single_quote = False
    in_double_quote = False

    for char_index, char in enumerate(line):
        if char == "'" and not in_double_quote:
            in_single_quote = not in_single_quote
        elif char == '"' and not in_single_quote:
            in_double_quote = not in_double_quote
        elif char == '#' and not in_single_quote and not in_double_quote:
            return line[:char_index]

    return line


def _split_shell_fragments(shellcmd):
    fragments = []

    for line in _logical_shell_lines(shellcmd):
        line = _strip_inline_comment(line).strip()

        if not line:
            continue

        for operator in SHELL_CONTROL_OPERATORS:
            line = line.replace(operator, '\n')

        fragments.extend(fragment.strip() for fragment in line.splitlines() if fragment.strip())

    return fragments


def _tokenize_shell_fragment(fragment):
    try:
        return shlex.split(fragment, comments=False, posix=True)
    except ValueError:
        return fragment.split()


def _is_snakemake_placeholder(token):
    return token.startswith('{') and token.endswith('}')


def _is_environment_assignment(token):
    return re.match(r'^[A-Za-z_][A-Za-z0-9_]*=', token) is not None


def _is_shell_redirection(token):
    return token.startswith(('>', '<', '2>', '&>'))


def _normalize_command_token(token):
    return token.strip().strip('(').strip()


def get_programs_from_shell_command(shellcmd):
    """Return external programs from a Snakemake shell command string."""

    programs = []

    for fragment in _split_shell_fragments(shellcmd or ''):
        for token in _tokenize_shell_fragment(fragment):
            token = _normalize_command_token(token)

            if not token:
                continue

            if token.endswith('\\'):
                token = token[:-1]

            if not token:
                continue

            if _is_environment_assignment(token):
                continue

            if _is_snakemake_placeholder(token):
                continue

            if _is_shell_redirection(token):
                continue

            if token in SHELL_TOKENS_TO_SKIP:
                break

            if token in SHELL_COMMANDS_TO_IGNORE:
                break

            if token.startswith('-') or token.startswith('{'):
                break

            programs.append(token)
            break

    return programs


def get_programs_from_snakemake_workflow(snakemake_workflow_object):
    programs = []
    programs_seen = set([])

    for rule in snakemake_workflow_object.rules:
        if not rule.shellcmd:
            continue

        for program in get_programs_from_shell_command(rule.shellcmd):
            if program in programs_seen:
                continue

            programs_seen.add(program)
            programs.append(program)

    return programs

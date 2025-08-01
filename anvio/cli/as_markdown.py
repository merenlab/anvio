#!/usr/bin/env python
# -*- coding: utf-8

import sys

import anvio

from anvio.errors import ConfigError


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['meren']
__provides__ = ['markdown-txt']
__description__ = ("Markdownizides TAB-delmited data with headers in terminal.")


def main():
    try:
        run_program()
    except ConfigError as e:
        print(e)
        sys.exit(-1)


def run_program():
    args = get_args()

    line_counter = 0
    num_fields = None
    field_indices_to_include = []
    field_indices_for_code_columns = []

    if args.max_num_lines_to_show:
        if args.max_num_lines_to_show < 2:
            raise ConfigError("Max number of lines to show can't be less than 2 :/")

    for line in sys.stdin:
        if '\t' not in line:
            raise ConfigError("This program only works with TAB delimited input files :/")

        line_counter += 1

        fields = [f.strip() for f in line.strip('\n').split('\t')]

        if line_counter == 1:
            if args.exclude_columns:
                field_names_to_exclude = [n.strip() for n in args.exclude_columns.split(',')]
                field_indices_to_include = [i for i in range(0, len(fields)) if fields[i] not in field_names_to_exclude]
            else:
                field_indices_to_include = range(0, len(fields))

            if args.code_columns:
                if args.code_columns == 'ALL-COLUMNS':
                    field_indices_for_code_columns = [i for i in range(0, len(fields))]
                else:
                    field_names_for_code_columns = [n.strip() for n in args.code_columns.split(',')]
                    field_indices_for_code_columns = [i for i in range(0, len(fields)) if fields[i] in field_names_for_code_columns]

        if not num_fields:
            num_fields = len(fields)
        else:
            if num_fields != len(fields):
                raise ConfigError(f"The number of TAB-delimited fields do not seem to be equal throughout the lines "
                                  f"of the input file. The first line had {num_fields} TAB-delimited fields, but the "
                                  f"line {line_counter} had {len(fields)} of them :/")

        if line_counter == 1:
            print("|" + "|".join([f"**`{fields[i]}`**" for i in field_indices_to_include]) + "|")
            print("|" + "|".join([":--"] * len(field_indices_to_include)) + "|")
        else:
            output = []

            for i in field_indices_to_include:
                if i in field_indices_for_code_columns:
                    output.append("`" + fields[i].replace('|', '\\|') + "`")
                else:
                    output.append(fields[i].replace('|', '\\|'))

            print("|" + "|".join(output) + "|")

        if line_counter == args.max_num_lines_to_show:
            print("|" + "|".join(["(...)"] * len(field_indices_to_include)) + "|")
            return


def get_args():
    from anvio.argparse import ArgumentParser
    parser = ArgumentParser(description=__description__)
    parser.add_argument('--max-num-lines-to-show', type=int, default=None, help="Maximum number of lines to show.")
    parser.add_argument('--exclude-columns', type=str, default=None, help="Comma-separated list of column names to exclude from the output.")
    parser.add_argument('--code-columns', type=str, default=None, help="Comma-separated list of column names to be wrapped with '`' characters. "
                "If you pass 'ALL-COLUMNS' as a parameter rather than individual ones, all columns will be wrapped with '`' characters because "
                "you're welcome.")
    return parser.get_args(parser)


if __name__ == '__main__':
    main()

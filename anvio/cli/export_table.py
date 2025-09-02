#!/usr/bin/env python
# -*- coding: utf-8
"""Takes a database and a table name as parameters, stores the content as
a TAB-delimited matrix."""

import sys

import anvio
import anvio.db as db
import anvio.utils as utils
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths

from anvio.errors import ConfigError, FilesNPathsError


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['meren']
__description__ = "Export anvi'o database tables as TAB-delimited text files"


def main():
    args = get_args()
    run = terminal.Run()

    A = lambda x: args.__dict__[x] if x in args.__dict__ else None
    args_db_path = A('db_path')
    args_table = A('table')
    args_list = A('list')
    args_fields = A('fields')
    args_output_file = A('output_file')
    args_matrix_format = A('matrix_format')
    args_index = A('index')
    args_columns = A('columns')
    args_values = A('values')

    try:
        if args_index or args_columns or args_values:
            if not args_matrix_format:
                raise ConfigError("The parameters `--index`, `--column`, or `--values` are only relevant if you are also "
                                  "asking your output to be reported in `--matrix-format`.")
            if not (args_index and args_columns and args_values):
                raise ConfigError("If you define any of these parameters, then you must define all of them: "
                                  "`--index`, `--column`, `--values`.")

        filesnpaths.is_file_exists(args_db_path)
        database = db.DB(args_db_path, None, ignore_version=True)
        tables_in_database = database.get_table_names()

        if args_list and not args_table:
            for table in tables_in_database:
                print(table)
            sys.exit()

        if not args_table:
            raise ConfigError("You must specify a table name.")

        if args_table not in tables_in_database:
            raise ConfigError("Table '%s' is not seem to be in this databse :/" % args_table)

        run.info('Database', '"%s" has been initiated with its %d tables.' % (args_db_path, len(tables_in_database)))

        table_columns = database.get_table_structure(args_table)

        if args_list:
            run.info('Table columns', '"%s"' % ', '.join(table_columns))
            sys.exit()

        if args_fields:
            fields_of_interest = [f.strip() for f in args_fields.split(',')]
            table_columns = [f for f in table_columns if f in fields_of_interest]

            if not len(table_columns):
                raise ConfigError("None of the fields you are interested in are in the table header...")

            run.info('Columns to report', '"%s"' % ', '.join(table_columns))

        if args_index and len([v for v in [args_index, args_columns, args_values] if v not in table_columns]):
            raise ConfigError(f"Not all fields you have defined for `--index`, `--column`, `--values` are among "
                              f"the fields of your table. Here are your the field names you have: {', '.join(table_columns)}.")

        table_content = database.get_table_as_dataframe(args_table, columns_of_interest=table_columns)

        run.info('Table', '"%s" has been read with %d entries and %d columns.' % (args_table, len(table_content), len(table_columns)))

        if not args_output_file:
            args_output_file = args_table + '.txt'

        run.info('Matrix format', args_matrix_format, nl_before=1)
        if args_index:
            run.info('Index, columns, values were set', f"{args_index}, {args_columns}, {args_values}.", mc="green")

        if not args_matrix_format:
            utils.store_dataframe_as_TAB_delimited_file(table_content, args_output_file)
        else:
            if len(table_columns) != 3 and not args_index:
                raise ConfigError(f"The `--matrix-format` works automatically with tables that has three fields, but "
                                  f"this one has {len(table_columns)}. So you need to specify which fields should be "
                                  f"used for what when anvi'o tries to pivot this table and report it as a matrix "
                                  f"using the patameters `--index`, `--columns`, and `--values` in your next attempt. "
                                  f"Each of these parameters will have to be a field name, which you can choose from "
                                  f"these beautiful fields in your table: {', '.join(table_columns)}.")

            if not args_index:
                run.warning(f"Since you only used `--matrix-format` with a table that happened to have three fields, anvi'o "
                            f"is automatically using '{table_columns[0]}' as the 'index', '{table_columns[1]}' as the 'columns', and "
                            f"'{table_columns[2]}' as the 'values'. If you don't like those choices, or if the output does not "
                            f"make sense, consider setting these parameters manually.")
                table_as_matrix = table_content.pivot(index=table_columns[0], columns=table_columns[1], values=table_columns[2])
            else:
                table_as_matrix = table_content.pivot(index=args_index, columns=args_columns, values=args_values)

            utils.store_dataframe_as_TAB_delimited_file(table_as_matrix, args.output_file, include_index=True, index_label=table_columns[0])

        run.info('Output', args.output_file)
    except ConfigError as e:
        print(e)
        sys.exit(-1)
    except FilesNPathsError as e:
        print(e)
        sys.exit(-2)


def get_args():
    from anvio.argparse import ArgumentParser
    parser = ArgumentParser(description=__description__)

    parser.add_argument('db_path', metavar = 'DB',
                        help = "Anvi'o database to read from.")

    groupA = parser.add_argument_group('TARGET', "This is where you specify which table and/or fields you wisht o export "
                                                 "from your anvi'o database. If you don't know what tables you have in "
                                                 "your database, you can add the flag `--list` to your command. If your "
                                                 "command includes a table name, then the `--list` will show you the "
                                                 "field names found in that particular table.")
    groupA.add_argument(*anvio.A('table'), **anvio.K('table'))
    groupA.add_argument(*anvio.A('fields'), **anvio.K('fields'))
    groupA.add_argument(*anvio.A('list'), **anvio.K('list'))

    groupB = parser.add_argument_group('REPORTING', "This is where you can ask anvi'o to pivot your table that is in long format "
                                                    "and report it in wide format (i.e., as a matrix), instead. This may not work "
                                                    "every table right off the bat, and you may need to specify the `--index`, `--columns`, "
                                                    "and `--values` variables. Try it only with the `--matrix-format` flag, take a "
                                                    "look at the output, if it doesn't seem to be working for you, improve.")
    groupB.add_argument(*anvio.A('matrix-format'), **anvio.K('matrix-format'))
    groupB.add_argument('--index', type=str, default=None, help="The field name in your table be used as the index of your matrix output.")
    groupB.add_argument('--columns', type=str, default=None, help="The field name in your table be used as the column of your matrix output.")
    groupB.add_argument('--values', type=str, default=None, help="The field name in your table be used as the values of your matrix output.")

    groupC = parser.add_argument_group('OUTPUT', "Set the final resting place of these data.")
    groupC.add_argument(*anvio.A('output-file'), **anvio.K('output-file'))

    return parser.get_args(parser)


if __name__ == '__main__':
    main()

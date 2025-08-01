#!/usr/bin/env python
# -*- coding: utf-8

import sys
from anvio.argparse import ArgumentParser

import anvio
import anvio.utils as utils
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths
import anvio.interactive as interactive
from anvio.bottleroutes import BottleApplication

from anvio.errors import ConfigError, FilesNPathsError, DictIOError


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['ozcan', 'FlorianTrigodet', 'watsonar']
__requires__ = ['contigs-db']
__provides__ = ['contigs-stats', 'interactive', 'svg']
__description__ = "Start the anvi'o interactive interface for viewing or comparing contigs statistics."


def write_stats_to_file(tables, output_file):
    run = terminal.Run()
    filesnpaths.is_output_file_writable(output_file)
    output_table = [['contigs_db'] + tables['header']] + tables['basic_stats'] + tables['hmm'] + tables['scg']

    f = open(output_file, 'w')

    if anvio.AS_MARKDOWN:
        f.write(f"|{'|'.join(map(str, output_table[0]))}|\n")
        f.write(f"|{':--|' + '|'.join([':--:'] * (len(output_table[0][1:])))}|\n")
        for line in output_table[1:]:
            f.write(f"|{'|'.join(map(str, line))}|\n")
    else:
        for line in output_table:
            f.write("%s\n" % "\t".join(map(str, line)))

    f.close()

    run.info(f"{'Output file (as markdown)' if anvio.AS_MARKDOWN else 'Output file'}", output_file)


@terminal.time_program
def main():
    args = get_args()
    run = terminal.Run()

    A = lambda x: args.__dict__[x] if x in args.__dict__ else None

    try:
        utils.is_all_npm_packages_installed()

        if A('report_as_text'):
            if not A('output_file'):
                raise ConfigError("You set the --report-as-text flag but didn't give any --output-file. ")

        d = interactive.ContigsInteractive(args)

        if A('output_file'):
            write_stats_to_file(d.tables, A('output_file'))

        if not A('report_as_text'):
            if args.dry_run:
                run.info_single('Dry run, eh? Fine. Bai!', nl_after=1)
                sys.exit()

            args.mode = 'contigs'
            port_number = utils.get_port_num(args.port_number, args.ip_address, run=run)

            app = BottleApplication(d)
            app.run_application(args.ip_address, port_number)
    except ConfigError as e:
        print(e)
        sys.exit(-1)
    except FilesNPathsError as e:
        print(e)
        sys.exit(-2)
    except DictIOError as e:
        print(e)
        sys.exit(-3)


def get_args():
    parser = ArgumentParser(description=__description__)
    parser.add_argument('input', metavar = 'CONTIG DATABASE(S)', nargs='+',
                        help = "Anvio'o Contig databases to display statistics, \
                                you can give multiple databases by seperating them with space.")
    groupA = parser.add_argument_group('REPORT CONFIGURATION', "Specify what kind of output you want.")
    groupA.add_argument(*anvio.A('report-as-text'), **anvio.K('report-as-text'))
    groupA.add_argument(*anvio.A('output-file'), **anvio.K('output-file'))
    groupB = parser.add_argument_group('SERVER CONFIGURATION', "For power users.")
    groupB.add_argument(*anvio.A('dry-run'), **anvio.K('dry-run'))
    groupB.add_argument(*anvio.A('ip-address'), **anvio.K('ip-address'))
    groupB.add_argument(*anvio.A('port-number'), **anvio.K('port-number'))
    groupB.add_argument(*anvio.A('browser-path'), **anvio.K('browser-path'))
    groupB.add_argument(*anvio.A('server-only'), **anvio.K('server-only'))
    groupB.add_argument(*anvio.A('password-protected'), **anvio.K('password-protected'))

    return parser.get_args(parser)


if __name__ == '__main__':
    main()

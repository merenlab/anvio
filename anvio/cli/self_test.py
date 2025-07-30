#!/usr/bin/env python
# -*- coding: utf-8

import os
import sys
import glob
import shutil
import signal
import itertools
import subprocess

import anvio
import anvio.tests as tests
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths

from anvio.errors import FilesNPathsError, ConfigError


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['meren', 'semiller10', 'ekiefl', 'ivagljiva', 'mschecht']
__requires__ = []
__provides__ = []
__description__ = "A program for anvi'o to test itself"


def __catch_sig(signum, frame):
    """We need to press CTRL+C to kill the server that is run by this script in a\
       subprocess, but then we don't want the script itself to catch it. This is\
       a workaround to avoid that."""
    pass


def main():
    try:
        signal.signal(signal.SIGINT, __catch_sig)
        run_program()
    except FilesNPathsError as e:
        print(e)
        sys.exit(-1)
    except ConfigError as e:
        print(e)
        sys.exit(-2)


def run_program():
    args = get_args()
    run = terminal.Run()

    if args.suite not in tests.known_tests:
        raise ConfigError(f"Well, the test suite '{args.suite}' is not something anvi'o knows about :/ Here is a list "
                          f"of tests available if you would like to try again: {', '.join(list(tests.known_tests.keys()))}.")

    tests_dir_path = os.path.dirname(anvio.tests.__file__)

    test_files_found = [os.path.basename(p) for p in glob.glob(os.path.join(tests_dir_path, '*.sh'))]

    for test in list(itertools.chain(*list(tests.known_tests.values()))):
        if test not in test_files_found:
            raise FilesNPathsError(f"Anvi'o failed to locate the test file for '{test}' ('{tests.known_tests[test]}')... What a "
                                   f"terrible start :( (it was looking for it under the directory '{tests_dir_path}').")

    if args.output_dir:
        output_dir = os.path.abspath(args.output_dir)
        filesnpaths.check_output_directory(output_dir)
        filesnpaths.gen_output_directory(output_dir)
    else:
        output_dir = filesnpaths.get_temp_directory_path()

    os.chdir(tests_dir_path)

    try:
        for test in tests.known_tests[args.suite]:
            command = f'./{test} {output_dir} {"no_interactive" if args.no_interactive else "interactive"} {args.num_threads if args.num_threads > 1 else ""}'
            exitcode = subprocess.call(command, shell=True)
            if exitcode != 0:
                raise ConfigError(f"According to the exit code ('{exitcode}'), anvi'o suspects that something may have gone wrong while "
                                  f"running your tests :/ We hope that the reason is clear to you from the lines above. But if you "
                                  f"don't see anything obvious, and especially if the test ended up running until the end with "
                                  f"reasonable looking final results, you shouldn't worry too much about this error. Life "
                                  f"is short and we all can worry just a bit less.")
    except KeyboardInterrupt:
        pass

    if args.output_dir:
        run.info_single(f"The self-test is done, and all the files anvi'o generated are stored in {args.output_dir}", nl_after=1)
    else:
        shutil.rmtree(output_dir)
        run.info_single("Anvi'o's self-test is done, and the temporary files are all gone.", nl_after=1)


def get_args():
    from anvio.argparse import ArgumentParser
    parser = ArgumentParser(description=__description__)

    available_tests = ', '.join([f"'{test}'" for test in list(tests.known_tests.keys())])
    parser.add_argument('--suite', default='mini', help="A suite of component tests to execute. By default this program will "
                            f"execute the mini test of anvi'o, which will help you to see if your computer and installation is "
                            f"able to perform some of the most basic anvi'o operations, such as generating an anvi'o contigs "
                            f"database, profiling BAM files, or starting an interactive interface. But you are welcome to "
                            f"execute different component tests. Here is a list of what is available to you: "
                            f"{available_tests}")
    parser.add_argument(*anvio.A('output-dir'), **anvio.K('output-dir', {"help": "If you declare an output dir, all your data will be "
                            "stored in there, instead of being stored in a temporary directory to be deleted once the tests are done. "
                            "This is particularly useful if you wish to play with general anvi'o output files"}))
    parser.add_argument(*anvio.A('no-interactive'), **anvio.K('no-interactive'))
    parser.add_argument(*anvio.A('num-threads'), **anvio.K('num-threads'))

    return parser.get_args(parser)


if __name__ == '__main__':
    main()

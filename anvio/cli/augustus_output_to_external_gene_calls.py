#!/usr/bin/env python
# -*- coding: utf-8

import sys

import anvio
import anvio.utils as u
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths

from anvio.errors import ConfigError, FilesNPathsError


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['meren']
__provides__ = ['external-gene-calls']
__requires__ = ['augustus-gene-calls']
__description__ = ("Takes in gene calls by AUGUSTUS v3.3.3, generates an anvi'o external gene calls file. "
                   "It may work well with other versions of AUGUSTUS, too. It is just no one has tested the "
                   "script with different versions of the program")


def main():
    try:
        run_program()
    except ConfigError as e:
        print(e)
        sys.exit(-1)
    except FilesNPathsError as e:
        print(e)
        sys.exit(-2)


def run_program():
    args = get_args()
    run = terminal.Run()
    progress = terminal.Progress()
    pp = terminal.pretty_print

    gene_calls = {}

    filesnpaths.is_file_plain_text(args.input_file)
    filesnpaths.is_output_file_writable(args.output_file)

    progress.new("Parsing the input file")
    progress.update("...")

    L = lambda x: line.startswith(x)
    with open(args.input_file) as input_file:

        # quick test.
        if input_file.readline().find('AUGUSTUS') < 0:
            progress.reset()
            if args.just_do_it:
                run.warning("This file doesn't look like an AUGUSTUS gene calls file, but anvi'o will continue to "
                            "try to parse it as you instructed. If everything works well, please send this "
                            "file to an anvi'o developer so they can update this code to not complain about this.")
            else:
                raise ConfigError("This file doesn't look like an AUGUSTUS gene calls file :/ If you are certain "
                                  "that it is, you can try your chances with this script by adding the flag "
                                  "`--just-do-it` to your command.")

        gene_num = 0
        in_gene = False
        in_protein = False
        aa_sequence = ''

        while 1:
            if gene_num % 100 == 0:
                progress.update("gene %d ..." % gene_num)

            line = input_file.readline()
            if not line:
                progress.reset()
                run.info("Num genes parsed", pp(len(gene_calls)))

                break

            if L('# start gene'):
                in_gene = True

                l = input_file.readline()
                if '\t' not in l:
                    progress.reset()
                    raise ConfigError(f"The following gene entry does not seem to have any TAB characters for anvi'o "
                                      f"to parse fields accurately: \"{l}\"")

                gene_call = l.strip().split('\t')

                if len(gene_call) != 9:
                    progress.reset()
                    raise ConfigError(f"The entry does not seem to have as many fields as anvi'o expects: {gene_call} :/")

                gene_calls[gene_num] = {'contig': gene_call[0],
                                        'start': int(gene_call[3]) - 1,
                                        'stop': int(gene_call[4]),
                                        'direction': 'f' if gene_call[6] == '+' else 'r',
                                        'partial': 0,
                                        'call_type': 1,
                                        'source': 'AUGUSTUS',
                                        'version': 'v3.3.3',
                                        'aa_sequence': ''}

                continue

            if L('# end gene'):
                gene_calls[gene_num]['aa_sequence'] = aa_sequence

                in_gene = False
                in_protein = False
                aa_sequence = ''

                gene_num += 1
                continue

            if L('# protein sequence'):
                if line.strip().endswith(']'):
                    aa_sequence += line.strip()[line.find('[') + 1:-1]
                    in_protein = False
                else:
                    aa_sequence += line.strip()[line.find('[') + 1:]
                    in_protein = True

                continue

            if in_gene and in_protein:
                if line.strip().endswith(']'):
                    aa_sequence += line.strip()[2:-1]
                    in_protein = False
                else:
                    aa_sequence += line[2:].strip()

                continue

    progress.update("Storing the output ...")
    u.store_dict_as_TAB_delimited_file(gene_calls, args.output_file, ['gene_callers_id', 'contig', 'start', 'stop', 'direction', 'partial', 'call_type', 'source', 'version', 'aa_sequence'])
    progress.end()

    run.info("Anvi'o external genes file output", args.output_file, mc='green')

    if not len(gene_calls) and args.just_do_it:
        run.warning("You asked anvi'o to just to id, and she just did it. Anvi'o hopes that you enjoy your zero gene calls :(")


def get_args():
    from anvio.argparse import ArgumentParser
    parser = ArgumentParser(description=__description__)

    parser.add_argument("-i", "--input-file", required=True, help=("Gene calls file from AUGUSTUS (that ends with .gff). Please note that "
                                    "the script is only tested with AUGUSTUS v3.3.3 output (although it may still work with other versions of "
                                    "the program)."))
    parser.add_argument(*anvio.A("output-file"), **anvio.K("output-file"))
    parser.add_argument(*anvio.A("just-do-it"), **anvio.K("just-do-it"))

    return parser.get_args(parser)


if __name__ == '__main__':
    main()

#!/usr/bin/env python
# -*- coding: utf-8
"""Generate phylogenomic tree from aligment file"""

import os
import sys
from anvio.argparse import ArgumentParser

import anvio
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths
from anvio.fastalib import ReadFasta, SequenceSource
from anvio.drivers import driver_modules

from anvio.utils import check_contig_names
from anvio.errors import ConfigError, FilesNPathsError, DictIOError


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['ozcan', 'meren']
__requires__ = ['concatenated-gene-alignment-fasta']
__provides__ = ['phylogeny']
__resources__ = [("View this program in action in the anvi'o phylogenetics workflow", "http://merenlab.org/2017/06/07/phylogenomics/")]
__description__ = "Generate phylogenomic tree from aligment file"


def main():
    args = get_args()
    run = terminal.Run()
    progress = terminal.Progress()

    try:
        program = driver_modules['phylogeny'][args.program if (args.program and args.program in driver_modules) else 'default']

        input_file_path = os.path.abspath(args.fasta_file)
        filesnpaths.is_file_fasta_formatted(input_file_path)

        output_file_path = os.path.abspath(args.output_file)
        filesnpaths.is_output_file_writable(output_file_path)

        # check header validity
        if not args.just_do_it:
            progress.new('Checking FASTA deflines to avoid snafus later')
            progress.update('tick tock ...')
            fasta = SequenceSource(input_file_path)

            while next(fasta):
                if not check_contig_names(fasta.id, dont_raise=True):
                    progress.end()
                    raise ConfigError("At least one of the deflines in your FASTA file does not comply with the 'simple deflines' "
                                      "requirement of anvi'o :/ Anvi'o is very upset for making you do this, but you can simply "
                                      "solve this issue by either (1) removing all the characters after the first space character "
                                      "in the deflines of your sequences in the FASTA file, (2) by simplifying the deflines "
                                      "using `anvi-script-reformat-fasta` with `--simplify-deflines` flag, OR (3) by including the "
                                      "flag `--just-do-it` in your command line, in which case anvi'o will not check your deflines "
                                      "(in which case anvi'o shall not be held accountable for what might go wrong in your "
                                      "downstream analyses of the resulting newick file wiht anvi'o).")
                try:
                    int(fasta.id)
                    is_int = True
                except:
                    is_int = False

                if is_int:
                    progress.end()
                    raise ConfigError("At least one of the contigs in your FASTA file (well, this one to be precise: '%s') looks like "
                                      "a number. For reasons we can't really justify, anvi'o does not like those numeric names, and hereby "
                                      "asks you to make sure every contig name contains at least one alphanumeric character :/ Meanwhile we, "
                                      "the anvi'o developers, are both surprised by and thankful for your endless patience with such eccentric "
                                      "requests. You the real MVP." % fasta.id)

            fasta.close()
            progress.end()

        # We're clear
        run.info("Input aligment file path", input_file_path)
        run.info("Output file path", output_file_path)

        alignments = ReadFasta(input_file_path, quiet=True)
        run.info("Alignment names", ", ".join([i.split()[0] for i in alignments.ids]))

        alignment_lengths = [len(x) for x in alignments.sequences]
        if len(set(alignment_lengths)) == 1:
            run.info("Alignment sequence length", alignment_lengths[0])
        else:
            raise ConfigError("Alignment lengths are not equal in input file.")

        program().run_command(input_file_path, output_file_path)
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

    groupA = parser.add_argument_group('INPUT FILES', "Concatenated aligment files exported using anvi-get-sequences-for-gene-clusters")
    groupA.add_argument(*anvio.A('fasta-file'), **anvio.K('fasta-file', {'required': True}))

    groupB = parser.add_argument_group('OUTPUT FILE', "The output file where the generated newick tree will be stored.")
    groupB.add_argument(*anvio.A('output-file'), **anvio.K('output-file', {'required': True}))

    groupD = parser.add_argument_group('PROGRAM', "The program that will be used for generating tree. Available options: " + ", ".join(list(driver_modules['phylogeny'].keys())))
    groupD.add_argument(*anvio.A('program'), **anvio.K('program'))

    groupF = parser.add_argument_group('OTHERS', "Parameters of convenience")
    groupF.add_argument(*anvio.A('just-do-it'), **anvio.K('just-do-it'))
    return parser.get_args(parser)

if __name__ == '__main__':
    main()

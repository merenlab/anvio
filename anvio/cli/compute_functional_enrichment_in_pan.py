#!/usr/bin/env python
# -*- coding: utf-8
"""A program that computes functional enrichment in a pangenome"""

import sys

import anvio
import anvio.terminal as terminal
import anvio.summarizer as summarizer
import anvio.filesnpaths as filesnpaths

from anvio.errors import ConfigError, FilesNPathsError


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['ivagljiva', 'adw96', 'meren']
__resources__ = [("A description of the enrichment script run by this program can be found in Shaiber et al 2020", "https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-02195-w"),
                 ("An example of pangenome functional enrichment in the context of the Prochlorococcus metapangenome from Delmont and Eren 2018 is included in the pangenomics tutorial", "http://merenlab.org/2016/11/08/pangenomics-v2/")]
__tags__ = ["pangenomics", "functions"]
__requires__ = ['misc-data-layers', 'pan-db', 'genomes-storage-db', 'functions']
__provides__ = ['functional-enrichment-txt']
__description__ = ("A program that computes functional enrichment within a pangenome.")


@terminal.time_program
def main():
    args = get_args()

    # make sure the output files is not going to overwrite anything
    if filesnpaths.is_file_exists(args.output_file, dont_raise=True):
        raise ConfigError(f"There is already a file at '{args.output_file}' :/ Anvi'o has a thing against overwriting existing files. "
                          f"Please remove the exiting file first, or give a different output file name to this program.")

    # make sure we can write to the output file
    filesnpaths.is_output_file_writable(args.output_file)

    if not args.annotation_source and not args.list_annotation_sources:
        raise ConfigError("You must declare a function annotation source with `--annotation-source` parameter to be used for the "
                          "functional enrichment analysis (well, obviously). If you don't know what is available, you can use "
                          "`--list-annotation-sources` flag to find out what is there.")

    try:
        s = summarizer.PanSummarizer(args, lazy_init=True)
        s.functional_enrichment_stats()
    except ConfigError as e:
        print(e)
        sys.exit(-1)
    except FilesNPathsError as e:
        print(e)
        sys.exit(-2)


def get_args():
    from anvio.argparse import ArgumentParser

    parser = ArgumentParser(description=__description__)

    groupA = parser.add_argument_group('CRITICAL INPUTS', "Give anvi'o what it wants.")
    groupA.add_argument(*anvio.A('pan-db'), **anvio.K('pan-db'))
    groupA.add_argument(*anvio.A('genomes-storage'), **anvio.K('genomes-storage'))
    groupA.add_argument(*anvio.A('category-variable'), **anvio.K('category-variable'))

    groupB = parser.add_argument_group('FUNCTION ANNOTATION SOURCE', "Here you tell anvi'o which function annotation source to use.")
    groupB.add_argument(*anvio.A('annotation-source'), **anvio.K('annotation-source'))
    groupB.add_argument(*anvio.A('include-gc-identity-as-function'), **anvio.K('include-gc-identity-as-function'))
    groupB.add_argument(*anvio.A('list-annotation-sources'), **anvio.K('list-annotation-sources'))

    groupC = parser.add_argument_group('OUTPUT OPTIONS', "What comes out the other end. (Please provide at least the output file name.)")
    groupC.add_argument(*anvio.A('output-file'), **anvio.K('output-file', {'required': True}))

    groupD = parser.add_argument_group('OUTPUT OPTIONS FOR FUNCTIONAL ENRICHMENT', "Reporting options that only make sense for input option #1 or #3.")
    groupD.add_argument(*anvio.A('functional-occurrence-table-output'), **anvio.K('functional-occurrence-table-output'))

    groupE = parser.add_argument_group('OPTIONAL THINGIES', "If you want it, here it is, come and get it.")
    groupE.add_argument(*anvio.A('just-do-it'), **anvio.K('just-do-it'))

    return parser.get_args(parser)


if __name__ == '__main__':
    main()

#!/usr/bin/env python
# -*- coding: utf-8
"""A program that computes metabolic enrichment across groups of genomes and metagenomes"""

import sys

import anvio
import anvio.kegg as kegg
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths

from anvio.errors import ConfigError, FilesNPathsError

__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['ivagljiva', 'adw96']
__resources__ = [("A description of the enrichment script run by this program can be found in Shaiber et al 2020", "https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-02195-w"),
                  ("An example of pangenome functional enrichment in the context of the Prochlorococcus metapangenome from Delmont and Eren 2018 is included in the pangenomics tutorial", "http://merenlab.org/2016/11/08/pangenomics-v2/")]
__tags__ = ["metabolism"]
__requires__ = ['kegg-metabolism', 'user-metabolism', 'groups-txt', 'external-genomes', 'internal-genomes', 'functions']
__provides__ = ['functional-enrichment-txt']
__description__ = ("A program that computes metabolic enrichment across groups of genomes and metagenomes")


@terminal.time_program
def main():
    args = get_args()

    try:
        filesnpaths.is_output_file_writable(args.output_file)

        e = kegg.KeggModuleEnrichment(args)
        e.run_enrichment_stats()
    except ConfigError as e:
        print(e)
        sys.exit(-1)
    except FilesNPathsError as e:
        print(e)
        sys.exit(-2)


def get_args():
    from anvio.argparse import ArgumentParser

    parser = ArgumentParser(description=__description__)

    groupB = parser.add_argument_group('CRITICAL INPUTS', "This program requires TWO INPUT FILES: (1) an output file produced by the "
                                       "program `anvi-estimate-metabolism` as input and (2) a groups file (using the `--groups-txt` "
                                       "parameter) to specify which sample is in which group.")
    groupB.add_argument(*anvio.A('modules-txt'), **anvio.K('modules-txt'))
    groupB.add_argument(*anvio.A('groups-txt'), **anvio.K('groups-txt'))

    groupO = parser.add_argument_group('OUTPUT OPTIONS', "What comes out the other end. (Please provide at least the output file name.)")
    groupO.add_argument(*anvio.A('output-file'), **anvio.K('output-file', {'required': True}))

    groupE = parser.add_argument_group('OPTIONAL OPTIONS', "If you want it, here it is, come and get it.")
    groupE.add_argument(*anvio.A('sample-header'), **anvio.K('sample-header'))
    groupE.add_argument(*anvio.A('module-completion-threshold'), **anvio.K('module-completion-threshold',
                                {'help': "This threshold defines the percent completeness score at which we consider a KEGG module to be 'present'"
                                         "in a given sample. That is, if a module's completeness in a sample is less than this value, then we say "
                                         "the module is not present in that sample, and this will affect the module's enrichment score. "
                                         "By extension, if a module's completeness is less than this value in all samples, it will have "
                                         "a very very low enrichment score (ie, it will not be considered enriched at all, because it doesn't occur in "
                                         "any groups). Note that the closer this number is to 0, the more meaningless this whole enrichment analysis is... "
                                         "but hey, this is your show. This threshold CAN be different from the completeness threshold used in `anvi-estimate-metabolism` "
                                         "if you wish. The default threshold is %(default)g."}))
    groupE.add_argument(*anvio.A('use-stepwise-completeness'), **anvio.K('use-stepwise-completeness'))
    groupE.add_argument(*anvio.A('include-samples-missing-from-groups-txt'), **anvio.K('include-samples-missing-from-groups-txt'))
    groupE.add_argument(*anvio.A('qlambda'), **anvio.K('qlambda'))
    groupE.add_argument(*anvio.A('just-do-it'), **anvio.K('just-do-it'))

    return parser.get_args(parser)


if __name__ == '__main__':
    main()

#!/usr/bin/env python
# -*- coding: utf-8

import os
import sys
import anvio
import argparse

import pandas as pd
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths

from anvio.dbinfo import DBInfo as dbi
from anvio.dbops import ContigsDatabase
from anvio.errors import ConfigError, FilesNPathsError
from anvio.kegg import KeggContext, KeggMetabolismEstimator

__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['ivagljiva', 'meren']
__requires__ = ["contigs-db"]
__provides__ = ["metabolic-independence-score"]
__description__ = "Takes a genome as a contigs-db, and tells you whether it can be considered as an organism of high metabolic independence, or not"

# list of metablic modules general completion of which suggests high metabolic independence
# please see docs/programs/anvi-script-classify-hmi-genomes.md for more information. Also,
# ANY CHANGE HERE must be reflected in that help doc, too:
HMI_mods = ['M00049', 'M00050', 'M00007', 'M00140', 'M00005', 'M00083', 'M00120', 'M00854',
            'M00527', 'M00096', 'M00048', 'M00855', 'M00022', 'M00844', 'M00051', 'M00082',
            'M00157', 'M00026', 'M00526', 'M00015', 'M00019', 'M00432', 'M00018', 'M00570',
            'M00126', 'M00115', 'M00028', 'M00924', 'M00122', 'M00125', 'M00023', 'M00631',
            'M00061' ]


def main():
    try:
        run_program()
    except ConfigError as e:
        print(e)
        sys.exit(-1)
    except FilesNPathsError as e:
        print(e)
        sys.exit(-1)


def get_metabolic_independence_estimates(metabolism_file, HMI_module_list, threshold, strategy='pathwise'):
    """Adds up the pathwise completeness scores of each module in the HMI list. Returns the class, and the score."""

    module_df = pd.read_csv(metabolism_file, sep="\t", index_col=0)

    module_not_found = [m for m in HMI_module_list if m not in module_df.index]
    if module_not_found:
        mod_str = ", ".join(module_not_found)
        raise ConfigError(f"Some of the module accession(s) in the provided list were not found in the metabolism "
                          f"output files, which means that they do not exist in the KEGG database you "
                          f"provided. You either have the wrong accession(s), or you need to provide a different "
                          f"database. Here is a list of the missing modules: {mod_str}")

    score_key = f"{strategy}_module_completeness"
    subset = module_df.loc[HMI_module_list][score_key]

    genome_score = subset.sum()

    if genome_score >= threshold:
        return ('High', genome_score)
    else:
        return ('Low', genome_score)


def run_program():
    args = get_args()
    run = terminal.Run()
    P = terminal.pluralize
    pp = terminal.pretty_print

    # learn about the contigs-db and make sure it has KOfams
    dbi(args.contigs_db, expecting='contigs')
    contigs_db = ContigsDatabase(args.contigs_db)
    run.info('Contigs DB', f"{args.contigs_db} ('{contigs_db.meta['project_name']}' with "
                           f"{P('contig', contigs_db.meta['num_contigs'])} and "
                           f"{pp(contigs_db.meta['total_length'])} nucleotides)")

    if 'KEGG_Module' not in contigs_db.meta['gene_function_sources']:
        raise ConfigError("A minor inconvenience: your genome does not seem to be annotated with KEGG. Please run "
                          "the program `anvi-run-kegg-kofams` on it, and try again.")
    contigs_db.disconnect()

    # make sure KEGG data is accessible in this environment
    kegg_ctx = KeggContext(args)

    if not os.path.exists(kegg_ctx.kegg_modules_db_path):
        if args.kegg_data_dir:
            raise ConfigError("Anvi'o is unable to find the MODULES.db in the KEGG data directory you have specified :/")
        else:
            raise ConfigError("Anvi'o is unable to find the MODULES.db in the default KEGG data directory. You either "
                              "need to run `anvi-setup-kegg-data` to install the KEGG data in your environment, or "
                              "you need to use the `--kegg-data-dir` parameter to specify where should anvi'o find it.")

    run.info("KEGG data directory", kegg_ctx.kegg_data_dir)

    global HMI_mods
    # take care of the modules of interest
    if args.module_list:
        filesnpaths.is_file_plain_text(args.module_list)
        HMI_mods = [x.strip() for x in open(args.module_list, 'r').readlines()]
    else:
        HMI_mods = HMI_mods

    run.info("Modules of interest", f"{', '.join(HMI_mods)}")

    # make sure the threshold makes sense
    try:
        float(args.threshold)
        assert(args.threshold > 0)
    except:
        raise ConfigError(f"The threshold value must be a positive number. '{args.threshold}' doesn't look like it.")

    run.info("Threshold for HMI/LMI classification", args.threshold)

    # select which type of completeness score to use
    strategy = 'pathwise'
    if args.use_stepwise_completeness:
        strategy = 'stepwise'
    run.info("Score type", f"{strategy} completeness")

    # keep track of generated files so we know what to clean up at the end
    files_for_cleanup = []

    # prepare for running metabolism estimation
    db_prefix = os.path.basename(args.contigs_db).split('.')[0]
    metabolic_modules_file = f"{db_prefix}_modules.txt"

    filesnpaths.is_output_file_writable(metabolic_modules_file, ok_if_exists=False)

    files_for_cleanup.append(metabolic_modules_file)
    m = KeggMetabolismEstimator(argparse.Namespace(contigs_db=args.contigs_db,
                                                   kegg_data_dir=args.kegg_data_dir,
                                                   include_zeros=True,
                                                   output_file_prefix=db_prefix), run=terminal.Run(verbose=False))
    m.estimate_metabolism()

    if not os.path.exists(metabolic_modules_file):
        raise ConfigError("Anvi'o is unable to find the output from `KeggMetabolismEstimator` :( Something must have gone wrong "
                          "prior to this step, but at this stage of the code we're Jon Snow.")



    genome_metabolic_independence_class, genome_metabolic_independence_score = get_metabolic_independence_estimates(metabolic_modules_file, HMI_mods,
                                                                                                                    args.threshold, strategy=strategy)

    run.warning(None, header='CLASSIFICATION RESULT', lc='green')
    run.info('Metabolic independence', genome_metabolic_independence_class, mc='green')
    run.info('Threshold for classification', args.threshold)
    run.info('Genome score', f"{genome_metabolic_independence_score:.2f}")
    run.info('Completeness metric used', f"{strategy}")


def get_args():
    from anvio.argparse import ArgumentParser

    parser = ArgumentParser(description=__description__)

    groupA = parser.add_argument_group('THE GENOME', "Your genome as a contigs-db.")
    groupA.add_argument(*anvio.A('contigs-db'), **anvio.K('contigs-db'))

    groupE = parser.add_argument_group('OPTIONAL INPUTS', "To go beyond what is default.")
    groupE.add_argument('--module-list', metavar='FILE', help="By default, this script will try to determine whether your genome "
                "represents an HMI or LMI organism based on the sum of the completeness scores of each metabolic "
                "module we have previously identified as predictors of metabolic independence status. You can find the list of "
                "such modules in the documentation of the tool (please see the link at the end of the help menu). But using this "
                "parameter, you can provide a list of KEGG module accession numbers (one accession per line) to be used instead.")
    groupE.add_argument(*anvio.A('kegg-data-dir'), **anvio.K('kegg-data-dir', {'help': "By default, anvi'o will try to use the "
                "default KEGG data on your computer (if you have it setup already). Alternatively, you can use this parameter to "
                "point the path that contains a MODULES.db that holds the HMI module data."}))
    groupE.add_argument('--threshold', metavar='NUM', type=float, help="Genomes with a score over this threshold will be classified "
                "as HMI. Must be a positive number. The default value reflects our findings based on known HMI/LMI genomes. Please see the docs for "
                "more information on this.", default=20)
    groupE.add_argument('--use-stepwise-completeness', default=False, action='store_true', help="Pass this flag if you want stepwise "
                        "completeness scores to be used for the calculation rather than pathwise completeness. Confused about what these are? "
                        "Please see the help page for anvi-estimate-metabolism: https://anvio.org/help/main/programs/anvi-estimate-metabolism")

    groupO = parser.add_argument_group('OPTIONAL THINGIES', "We are not sure if increasing the number of threads will help at this "
                                                            "time. But you can try.")
    groupO.add_argument(*anvio.A('num-threads'), **anvio.K('num-threads'))


    return parser.get_args(parser)


if __name__ == '__main__':
    main()

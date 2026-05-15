#!/usr/bin/env python
# -*- coding: utf-8

import sys
import shutil

from anvio.argparse import ArgumentParser
import anvio
import anvio.dbops as dbops
import anvio.utils as utils
import anvio.terminal as terminal
import anvio.summarizer as summarizer
import anvio.filesnpaths as filesnpaths


from anvio.errors import ConfigError, FilesNPathsError

from anvio.panrep import PanRepresenter


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ["mahmoudyousef98"]
__requires__ = ["pan-db", "genomes-storage-db"]
__provides__ = []
__resources__ = [
    (
        "The role of gene cluster homogeneity described in the Anvi'o pangenomics tutorial",
        "http://merenlab.org/2016/11/08/pangenomics-v2/#inferring-the-homogeneity-of-gene-clusters",
    )
]
__description__ = "Compute homogeneity for gene clusters"


@terminal.time_program
def main():
    run_program()


def run_program():
    args = get_args()
    run = terminal.Run()
    progress = terminal.Progress()

    temp_dir = filesnpaths.get_temp_directory_path()
    args.contigs_fasta = f"{temp_dir}/My_Perfect_FASTA.fasta"
    args.external_gene_calls = f"{temp_dir}/gene_call_ids.txt"
    args.skip_mindful_splitting = True
    args.db_variant = "pangenome_supplemented"
    args.contigs_db = args.output_file

    if args.keep_promoter:
        args.keep_senteny = True
        run.info_single(
            "Since you chose to keep the promoter region that means you also keep the  senteny by definition"
        )

    myPanRep = PanRepresenter(args, temp_dir)
    while myPanRep.all_clusters:
        if myPanRep.first_iteration:
            myPanRep.process_best_genome()

        myPanRep.process_additional_genomes()
        myPanRep.all_clusters -= myPanRep.seen_clusters
    # end while

    myPanRep.assemble_sumplementary_contig()
    myPanRep.write_outputs()
    myPanRep.build_contigs_db()

    if anvio.DEBUG:
        run.warning(
            "The temp directory, %s, is kept. Please don't forget to clean it up "
            "later" % temp_dir,
            header="Debug",
        )
    else:
        run.info_single(
            "Cleaning up the temp directory (you can use `--debug` if you would "
            "like to keep it for testing purposes)",
            nl_before=1,
            nl_after=1,
        )
        shutil.rmtree(temp_dir)


def get_args():
    parser = ArgumentParser(description=__description__)

    groupA = parser.add_argument_group(
        "INPUT FILES", "Input files from the pangenome analysis."
    )
    groupA.add_argument(*anvio.A("external-genomes"), **anvio.K("external-genomes"))
    groupA.add_argument(*anvio.A("pan-db"), **anvio.K("pan-db"))
    groupA.add_argument(*anvio.A("genomes-storage"), **anvio.K("genomes-storage"))

    groupB = parser.add_argument_group(
        "REPORTING",
        "How do you want results to be reported? Anvi'o can produce a TAB-delimited output file for\
        you (for which you would have to provide an output file name). Or the results can be stored\
        in the pan database directly, for which you would have to explicitly ask for it. You can get\
        both as well in case you are a fan of redundancy and poor data analysis practices. Anvi'o\
        does not judge.",
    )
    groupB.add_argument(*anvio.A("output-file"), **anvio.K("output-file"))

    parser.add_argument("--kmer-size", default=4)
    parser.add_argument("--split-length", default=20000)
    parser.add_argument("--best-genome")  # check if it exists
    parser.add_argument("--min-num-contigs")  # check if it's actually an int
    parser.add_argument("--keep-senteny", action="store_true")
    parser.add_argument("--keep-promoter", action="store_true")

    return parser.get_args(parser)


if __name__ == "__main__":
    main()

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
    project_name = args.project_name
    temp_dir = filesnpaths.get_temp_directory_path()
    args.contigs_fasta = f"{temp_dir}/{project_name}.fasta"
    args.external_gene_calls = f"{temp_dir}/{project_name}_gene_calls.txt"
    args.db_variant = "pan-genome"

    if args.keep_promoter:
        args.keep_senteny = True
        run.info_single( "Since you chose to keep the promoter region that means you also keep the  senteny by definition")


    try:
        int(args.max_num_contigs)
    except ValueError as e:
        raise ConfigError("max-num-contigs should be a strictly positive integer")

    try:
        int(args.gap_size)
    except ValueError as e:
        raise ConfigError("gap-size should be a strictly positive integer")

    try:
        float(args.gap_size)
    except ValueError as e:
        raise ConfigError("alpha should be between 0 and 1 inclusive")

    try:
        float(args.gap_size)
    except ValueError as e:
        raise ConfigError("alpha should be between 0 and 1 inclusive")

#FilesNPathsError

    myPanRep = PanRepresenter(args, temp_dir)
    # print(f" This is it {myPanRep.external_genomes.__dict__.keys()}")
    # print(f" This is contigs_dbs_found {myPanRep.external_genomes.contigs_dbs_found}")
    # print(f" This is External_genome_names {myPanRep.external_genomes}")
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

    groupA = parser.add_argument_group("INPUT FILES", "Input files from the pangenome analysis.")
    groupA.add_argument(*anvio.A("external-genomes"), **anvio.K("external-genomes", {'required': True}))
    groupA.add_argument(*anvio.A("genomes-storage"), **anvio.K("genomes-storage", {'required': True}))
    groupA.add_argument(*anvio.A("pan-db"), **anvio.K("pan-db"))

    groupB = parser.add_argument_group("MY STUFF", "All these option are related the contigs_db that will be outputed.")
    groupB.add_argument("--gap-size", metavar= 'INT', default=20)
    groupB.add_argument("--alpha", metavar= 'FLOAT', default=0.8)
    groupB.add_argument("--representative", metavar= 'GENOME-NAME')  # check if it exists
    groupB.add_argument("--max-num-contigs", metavar= 'INT', default = 99999)  # check if it's actually an int
    groupB.add_argument("--keep-synteny", action="store_true")
    groupB.add_argument("--keep-promoter", action="store_true")

    groupB.add_argument(
        "--test-flag",
        metavar="FILE_PATH",
        help="The Help text goes HERE",
    )

    groupC = parser.add_argument_group(
        "OUTPUT-STUFF",
        "All these option are related the contigs_db that will be outputed.",
    )
    groupC.add_argument(*anvio.A("output-file"), required=True, **anvio.K("output-file"))
    groupC.add_argument(*anvio.A("project-name"), **anvio.K("project-name"))
    groupC.add_argument(*anvio.A("description"), **anvio.K("description"))
    groupC.add_argument(*anvio.A("kmer-size"), **anvio.K("kmer-size"))
    groupC.add_argument(*anvio.A("split-length"), **anvio.K("split-length"))
    groupC.add_argument(*anvio.A("skip-mindful-splitting"), **anvio.K("skip-mindful-splitting"))

    return parser.get_args(parser)


if __name__ == "__main__":
    main()

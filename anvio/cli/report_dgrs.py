#!/usr/bin/env python
# -*- coding: utf-8
import sys
import argparse

import anvio
import anvio.dgrs as dgrs
import anvio.terminal as terminal

# debugging + profiling
import io
import cProfile, pstats

from anvio.errors import ConfigError, FilesNPathsError
from anvio.terminal import time_program, Run

profiler = cProfile.Profile()

with terminal.SuppressAllOutput():
    import anvio.data.hmm as hmm_data
    available_hmm_sources = list(hmm_data.sources.keys())



__author__ = "Developers of anvi'o (see AUTHORS.txt)"
__copyright__ = "Copyleft 2015-2024, the Meren Lab (http://merenlab.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ["katysloz"]
__requires__ = ["fasta", "contigs-db", 'profile-db']
__provides__ = []
__description__ = "A program designed to find Diversity-generating retroelements based on high single nucleotide variants and  nucleotide similarity"

run = Run()

@time_program

def main():
    args = get_args()
    try:
        if args.hmm_usage:
            args.hmm_usage = [p.strip() for p in args.hmm_usage.split(',') if p.strip()]
        D = dgrs.DGR_Finder(args)
        D.process(args)
    except ConfigError as e:
        print(e)
        sys.exit(-1)
    except FilesNPathsError as e:
        print(e)
        sys.exit(-1)


#call function here for filtering hits of blast from root
def float_range(mini,maxi):
    """Return function handle of an argument type function for
        ArgumentParser checking a float range: mini <= arg <= maxi
        mini - minimum acceptable argument
        maxi - maximum acceptable argument"""

    # Define the function with default arguments
    def float_range_checker(arg):
        """New Type function for argparse - a float within predefined range."""

        try:
            f = float(arg)
        except ValueError:
            raise argparse.ArgumentTypeError("must be a decimal number")
        if f < mini or f > maxi:
            raise argparse.ArgumentTypeError("must be in range [" + str(mini) + "-" + str(maxi)+"]")
        return f

    # Return function handle to checking function
    return float_range_checker

my_float_range= float_range(0.5,1.0)


def get_args():

    available_hmm_sources_pretty = '; '.join([f"'{s}' (type: {hmm_data.sources[s]['kind']})" for s in available_hmm_sources])

    from anvio.argparse import ArgumentParser

    parser = ArgumentParser(description=__description__)

    groupX = parser.add_argument_group('PARAMETERS OF CONVENIENCE', "Because life is already very hard as it is.")
    groupX.add_argument(*anvio.A('just-do-it'), **anvio.K('just-do-it'))
    groupX.add_argument(*anvio.A('num-threads'), **anvio.K('num-threads'))
    groupX.add_argument(*anvio.A('verbose'), **anvio.K('verbose'))

    groupA = parser.add_argument_group('INPUT DATA', "Contigs.db AND a Profile.db preferably a merged profile.db. The tool searches for DGRs based on areas of high SNV density")
    groupA.add_argument(*anvio.A('contigs-db'), **anvio.K('contigs-db', {'required': True}))
    groupA.add_argument(*anvio.A('profile-db'), **anvio.K('profile-db', {'required': True})) #ADD THAT HAS TO BE MERGED_DB, unless have

    groupE = parser.add_argument_group('CONTIGS AND PROFILE DB INPUT ARGUMENTS', "Options for using the Contigs.db and Profile.db input for this program")
    groupE.add_argument("-I","--hmm-usage", required = False, help="The name of the HMM run with your Contigs.db. If not provided, "
                        "anvi'o will use the 'Reverse_Transcriptase' HMM by default (type: 6 clades of DGR Retroelements from "
                        "doi.org/10.1038/s41467-021-23402-7 including other known reverse transcriptases). You can provide a "
                        "comma-separated list of names for multiple profiles (but in that case don't put a space between each "
                        f"profile name). As a reminder here is the list of anvi'o installed profiles available to you: {available_hmm_sources_pretty}.",
                        type=str, default = None)
    groupE.add_argument(*anvio.A('gene-caller'), **anvio.K('gene-caller', {'help': "The gene caller to show gene calls if you are using a contigs.db. This is used to tell the program that you want to find the genes that your Variable Regions occurred in."}))

    groupG = parser.add_argument_group('COLLECTIONS MODE', "Options for collections mode. This is when you want to restrict your search to an individual collection, where anvi'o will look in bins for DGRs, importantly working using contigs not splits.")
    groupG.add_argument("--collections-mode", help="The flag to use that searches through specified splits that are in the collection that is specified with the '-C / --collection' flag.", action = "store_true")
    groupG.add_argument(*anvio.A('collection-name'), **anvio.K('collection-name', {'help':"The name of the singular collection that you want to search for DGRs in"}))

    groupJ = parser.add_argument_group('SNV CLUSTER LOCATOR PARAMETERS', "Options for adjusting how to locate clusters of SNVs")
    groupJ.add_argument("-d","--departure-from-reference-percentage", help="Minimum departure from reference to consider a SNV. Default is 0.1", type=float, default=0.1)
    groupJ.add_argument("--minimum-snv-density", help="The minimum percentage of SNVs over a SNV cluster needed for it to be valid, think about it as in the number of SNVs over the length of a possible VR, Default = 0.1", type=float, default=0.1)
    groupJ.add_argument("-w","--snv-window-size", help="Size of sliding window in bp for SNV density calculation. Default = 50", type=int, default=50, metavar="INT")
    groupJ.add_argument("--snv-window-step", help="Step size in bp for sliding window movement. Default = 10", type=int, default=10, metavar="INT")
    groupJ.add_argument("--variable-buffer-length", help="Length of bp added to your high SNV density 'window'. Default = 35", type=int, default=35)

    groupB = parser.add_argument_group('BLASTN ARGUMENTS', "BLASTn parameters for potential Template and Variable Region search")
    groupB.add_argument("--word-size", help="BLASTn word size parameter. Default = 8", type=str, default=8, metavar="INT")

    groupC = parser.add_argument_group('LOCATING VR OPTIONS', "Options for the fine tuning this program's location of variable regions'")
    groupC.add_argument("--skip-Ns", help="Skip 'N' bases when searching for mismatches", action = 'store_true', default=True)
    groupC.add_argument("--skip-dashes", help="Skip '-' bases when searching for mismatches", action = 'store_true',default=True)
    groupC.add_argument("--discovery-mode", help="By default, anvi'o uses SNVs occurring in the first and second codon position of ORF to identify DGRs. "
                        "This constraint allows for a fast search and more reliable results. If you feel daring, you can use this flag and let anvi'o use "
                        "ANY SNVs to identify regions of interest for the VR/TR search.", action = "store_true", default=False)

    groupD = parser.add_argument_group('FILTER BLAST RESULTS FOR TR/VRS', "Parameters for refining how stringent your search for template and variable regions is")
    groupD.add_argument("--repeat-threshold", help="This is how repeated a VR or TR sequence can be after being masked by [pytantan](https://doi.org/10.1093/nar/gkq1212). The default is 50%%.", type=float, default=0.5)
    groupD.add_argument("--num-imperfect-tandem-repeats", help="Number of imperfect tandem repeats from the BLAST search in both the TR and VR regions, this is computed using [pytrf](https://github.com/lmdu/pytrf). The default is 10", type=int, default=10)
    groupD.add_argument("--repeat-motif-coverage", help="Defines the coverage threshold for imperfect repeat motifs in the VR sequence. This is calculated as (motif length × number of repeats) / total VR sequence length. Imperfect tandem repeats from the BLAST search are computed using [pytrf](https://github.com/lmdu/pytrf). The default is 0.8", default=0.8)
    groupD.add_argument("--initial-mismatch-bias-threshold", help="Minimum fraction of mismatches to dominant base for initial BLAST hit filtering (before trimming). This is the first stage of a two-stage filter. Default = 0.6", type=float, default=0.6, metavar="FLOAT")
    groupD.add_argument("--trimmed-mismatch-bias-threshold", help="Minimum fraction of mismatches to dominant base required in the trimmed optimal window. This is the second stage that finds the best region within each BLAST hit. Default = 0.95", type=float, default=0.95, metavar="FLOAT")
    groupD.add_argument("--minimum-vr-length", help="Minimum length in bp for a valid VR/TR region after trimming. Default = 50", type=int, default=50, metavar="INT")
    groupD.add_argument("-n","--number-of-mismatches", help="Minimum number of mismatches required in the TR/VR alignment. Default = 7", type=int, default=7, metavar="INT")
    groupD.add_argument("--allow-any-base", help="By default, anvi'o only reports DGRs with A-base mutagenesis (the canonical mechanism). "
                        "Use this flag to also report DGRs with other dominant bases (C, G, or T). This is useful for exploratory "
                        "analysis but may include false positives.", default=False, action="store_true")
    groupD.add_argument("--min-mismatching-base-types-vr", help="The minimum number of mismatching base 'types' in the variable region, to ensure the variable region has multiple bases present. Has to be an integer between 1 and 4. (NB if 1 then has no variety in the mismatching bases). Default = 2", type=int, default=2, metavar="INT")
    groupD.add_argument("--min-base-types-tr", help="The minimum number of base 'types' in the template region, to ensure the variable region has multiple bases present. Has to be an integer between 1 and 4. (NB if 1 then has no variety in the mismatching bases). Default = 2", type=int, default=2, metavar="INT")
    groupD.add_argument("--min-base-types-vr", help="The minimum number of base 'types' in the variable region, to ensure the variable region has multiple bases present. Has to be an integer between 1 and 4. (NB if 1 then has no variety in the mismatching bases). Default = 2", type=int, default=2, metavar="INT")
    groupD.add_argument("--snv-matching-proportion", help="A blanket proportion of SNVs allowed in matching bases of the TR and VR. If there are more than 30 SNVs in the VR then 30%% of these can be in matching positions, if there are less than 30 SNVs then 25%% of these can be in matching positions. This is a conservative approach based on short-read metagenomes, that can be changed here.", type=float, metavar="FLOAT")
    groupD.add_argument("--snv-codon-position", help="The maximum percentage of SNVs that are allowed in the 3rd codon position of the variable region.", type=float, default=0.33, metavar="FLOAT")

    groupF = parser.add_argument_group('OUTPUT DIRECTORY', "Where to put all the output files.")
    groupF.add_argument(*anvio.A('output-dir'), **anvio.K('output-dir'))
    groupF.add_argument("--parameter-output", help="Add this flag if you want to output the parameters and their values you have input into a tsv file", action="store_true")
    groupF.add_argument(*anvio.A('overwrite-output-destinations'), **anvio.K('overwrite-output-destinations'))

    groupH = parser.add_argument_group('REPORTING GENOMIC CONTEXT AROUND DGRS',
                    "Once the DGRs are computed, anvi'o can go back to contigs where they are found, and for each DGR site report the surrounding genes and their functions as flat "
                    "text files in the output directory so you can actually take quick look at them.")
    groupH.add_argument('--skip-recovering-genomic-context', default=False, action="store_true", help="Of course you can skip this step because why should anyone have nice things.")
    groupH.add_argument('--num-genes-to-consider-in-context', default=3, type=int, help="With this parameter you can adjust the number of genes anvi'o should consider to characterize the genomic context "
                    "surrounding DGRs. If you set nothing here, anvi'o will go 3 genes upstream and downstream from each DGRs VR. If there are not enough genes in the contig given the position of either of the VRs, anvio' will report whatever it can.", metavar="INT")

    groupI = parser.add_argument_group('COMPUTING VARIABLE REGION VARIABILITY PROFILING', "What is the "
                    "proportion of variable regions represented across samples? A multi million question "
                    "that anvi'o WILL answer for you IF  you are using short-reads and you have `r1` and `r2` columns in a `samples-txt` "
                    "file that points to raw FASTQ reads you have used to generate your merged PROFILE.db. "
                    "Truly mind boggling stuff.")
    groupI.add_argument(*anvio.A('samples-txt'), **anvio.K('samples-txt'))
    groupI.add_argument('--initial-variable-primer-length', default=12, type=int, help="This is the initial section of the variable region primer "
                    "that actually comes before the variable region itself. This is to ensure that the whole primer is specific enough for the "
                    "primer search to find only the variable regions within the short reads. Default = 12 bp.", metavar="INT")
    groupI.add_argument('--whole-primer-length', default = 65,type=int, help="Calculating variable region ratios "
                    "requires anvi'o to 'design' an in silico primer, based on the variable region "
                    "and the upstream/downstream genomic context to search for short reads in raw sequencing data "
                    "to find the ratio of DGR activity per sample. This variable controls the length of a sequence "
                    "used to build this primer to search for short reads. The longer it is, the more specific "
                    "the primers will be to survey short reads. But if it is too long, "
                    "then there will not be enough reads to 'keep' depending on the minimum "
                    "short read length of the sequencing. The default value will probably good for the vast "
                    "majority of cases. Default = 65", metavar="INT")
    groupI.add_argument('--skip-compute-DGR-variability-profiling', default=False, action="store_true", help="This flag "
                    "will skip computing DGR variability profiling, which is an extremely costly step. It may be a good idea "
                    "to first run the workflow without computing activity, take a look at the output files to make sure "
                    "everything looks alright, and then run the workflow without this flag.")

    groupK = parser.add_argument_group('CALCULATE ACTIVITY FROM KNOWN DGRS?', "This input option is ONLY relevant to those who have "
                    "already run the entire workflow and have their dgrs output file (i.e., 'DGRs_found.tsv'). This input option will "
                    "use the existing dgrs reported in the input file and recalculate the dgr activity output. It is "
                    "a great way to calculate dgrs or refine them from a smaller set of genomes / metagenomes, and scale "
                    "up the characterizations of their activity to thousands of metagenomes quickly (*cough* cheater *cough*). "
                    "Please note that if you use this flag, most other options EXCEPT those that are listed below 'COMPUTING VARIABLE "
                    "REGION VARIABILITY PROFILING' will be irrelevant AND the flag for discovery-mode needs activating if the DGRs were found with this flag "
                    "or else we will only compute primers based on 1st and 2nd codon positions, so you can skip all and take a look at the parameters there. "
                    "Remember you need a samples-txt to compute DGR activity in those samples.")
    groupK.add_argument(*anvio.A('pre-computed-dgrs'), **anvio.K('pre-computed-dgrs'))

    return parser.get_args(parser)


if __name__ == '__main__':
    main()
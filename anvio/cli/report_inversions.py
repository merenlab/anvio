#!/usr/bin/env python
# -*- coding: utf-8

import sys

import anvio
import anvio.inversions as inversions

from anvio.errors import ConfigError, FilesNPathsError
from anvio.terminal import time_program, Run

__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['FlorianTrigodet', 'meren', 'ekiefl']
__requires__ = ["bams-and-profiles-txt"]
__provides__ = ["inversions-txt"]
__description__ = "Reports inversions"

run = Run()

@time_program
def main():
    args = get_args()

    try:
        I = inversions.Inversions(args)
        I.process()
    except ConfigError as e:
        print(e)
        sys.exit(-1)
    except FilesNPathsError as e:
        print(e)
        sys.exit(-1)


def get_args():
    from anvio.argparse import ArgumentParser

    parser = ArgumentParser(description=__description__)

    groupA = parser.add_argument_group('INPUT DATA', "Essentially a BAMs and profiles file and nothing more.")
    groupA.add_argument(*anvio.A('bams-and-profiles'), **anvio.K('bams-and-profiles'))

    groupB = parser.add_argument_group('CALCULATE ACTIVITY FROM KNOWN INVERSIONS?', "This input option is ONLY relevant to those who have "
                    "already run the entire workflow and have their consensus inversions (i.e., 'CONSENSUS-INVERSIONS.txt') or "
                    "sample-specific inversions (i.e., 'INVERSIONS-IN-[SAMPLE-NAME].txt') calculated. This input option will "
                    "use the existing inversions reported in the input file and recalculate the inversion activity output. It is "
                    "a great way to calculate inversions or refine them from a smaller set of genomes / metagenomes, and scale "
                    "up the characterizations of their activity to thousands of metagenomes quickly (*cough* cheater *cough*). "
                    "Please note that if you use this flag, most other options EXCEPT those that are listed below 'KEY ALGORITHMIC "
                    "COMPONENT 04: COMPUTING INVERSION ACTIVITY' will be irrelevant, so you can skip all and take a look at the "
                    "parameters there.")
    groupB.add_argument(*anvio.A('pre-computed-inversions'), **anvio.K('pre-computed-inversions'))

    groupC = parser.add_argument_group('KEY ALGORITHMIC COMPONENT 01: IDENTIFYING REGIONS OF INTEREST', "How should anvi'o identify regions of interest "
                    "based on REV/REV and FWD/FWD paired-end reads? Defaults will be good for most cases.")
    groupC.add_argument('--min-coverage-to-define-stretches', default=10, type=int, help="Value to break up contigs into 'stretches' of "
                    "high-coverage regions of FWD/FWD and REV/REV reads. The lower the value, the more noise. This acts as a low-pass "
                    "filter if it helps you imagine how it works.", metavar="INT")
    groupC.add_argument('--min-stretch-length', default=50, type=int, help="These are not the stretches you are looking for (unless they "
                    "are longer than this, obv).", metavar="INT")
    groupC.add_argument('--min-distance-between-independent-stretches', default=2000, type=int, help="Our 'low pass' filter may break a "
                    "single stretch of reasonable coverage of FWD/FWD and REV/REV reads into multiple pieces. To recover from that, we "
                    "wish to merge the fragmented ones if they are closer to one another than this value.", metavar="INT")
    groupC.add_argument('--num-nts-to-pad-a-stretch', default=100, type=int, help="Some leeway towards upstream and downstream context "
                    "that is essential to not miss key information due to coverage variation that may influence the beginnings and ends "
                    "of final stretches.", metavar="INT")

    groupD = parser.add_argument_group('KEY ALGORITHMIC COMPONENT 02: FINDING PALINDROMES', "Some essential parameters to find palindromes in sequence "
                    "stretches anvi'o identified in the previous step")
    groupD.add_argument(*anvio.A('min-palindrome-length'), **anvio.K('min-palindrome-length'))
    groupD.add_argument(*anvio.A('max-num-mismatches'), **anvio.K('max-num-mismatches'))
    groupD.add_argument(*anvio.A('min-distance'), **anvio.K('min-distance'))
    groupD.add_argument(*anvio.A('min-mismatch-distance-to-first-base'), **anvio.K('min-mismatch-distance-to-first-base'))
    groupD.add_argument(*anvio.A('palindrome-search-algorithm'), **anvio.K('palindrome-search-algorithm'))

    groupE = parser.add_argument_group('KEY ALGORITHMIC COMPONENT 03: CONFIRMING INVERSIONS', "Which palindromes are inversions? A one million dollar "
                    "question that is quite difficult to get right (but as anvi'o does get it right frequently).")
    groupE.add_argument('--check-all-palindromes', default=False, action="store_true", help="In a given 'stretch' "
                    "anvi'o has identified as a region of interest due to the coverage of FWD/FWD and REV/REV reads "
                    "(see KEY ALGORITHMIC COMPONENT 01), there may be multiple palindromic sequences. By default, "
                    "anvi'o will stop its search as soon as it finds evidence among short reads in the BAM file "
                    "that suggests a given palindrome represents an active inversion, WITHOUT testing any of the "
                    "remaining palindromes in the same region. While this strategy is extremely efficient as it "
                    "will not go through all reads over and over again for each palindrome, it will miss if there "
                    "are multiple inversions in a single stretch, which is extremely unlikely. But using this this "
                    "you can instruct anvi'o to test ALL palindromes in a stretch. We understand -- we are "
                    "unreasonable people, too.")
    groupE.add_argument('--process-only-inverted-reads', default=False, action="store_true", help="At one point, "
                    "anvi'o will have all the regions of interest in contigs that include palindromes that look "
                    "promising. At that point, it will access to short reads in the BAM file to determine which "
                    "palindromes in fact represent active inversions by searching for unique constructs that "
                    "can only occur when a genomic region did funny things. One option is to search for such "
                    "reads that are evidence of inversion activity among paired-end reads that are in FWD/FWD "
                    "or REV/REV orientation. Which is defined by the fetch-filter 'inversions'. However, this "
                    "may be too limiting as there may be paired-end reads that are too downstream or too upstream "
                    "to the region of interest, and thus mapping FWD/REV or REV/FWD orientation just like every "
                    "other non-rebellious paired end read, YET including one of the constructs. Thus, searching "
                    "all reads may in fact help the identification of more inversions (especially those that are "
                    "covered less), with only an added disadvantage of compute time, which should be negligible "
                    "in almost all instances (since anvi'o at this stage is only focusing on very specific "
                    "regions of genomes). But this parameter is here in case you insist on only using inverted "
                    "paired-end reads and assert your authority. You do you and turn on the flag, you rebellious "
                    "scientist who will likely miss a lot of additional inversions like a boss.")

    groupF = parser.add_argument_group('KEY ALGORITHMIC COMPONENT 04: COMPUTING INVERSION ACTIVITY', "What is the "
                    "proportion of invertible repeat orientation across samples? A two million dollars question "
                    "that anvi'o WILL answer for you IF you have `r1` and `r2` columns in your `bams-and-profiles-txt` "
                    "file that points to raw FASTQ reads you have used to generate your BAM files for each sample. "
                    "Truly crazy stuff.")
    groupF.add_argument(*anvio.A('num-threads'), **anvio.K('num-threads'))
    groupF.add_argument('--oligo-primer-base-length', default=12, type=int, help="Calculating inversion ratios "
                    "require anvi'o to 'design' an in silico primer based on palindromes associated with inversions "
                    "and the upstream/downstream genomic context to search for short reads in raw sequencing data"
                    "to find the ratio of inversion activity per sample. This variable is to control how much of "
                    "the palindrome for a given inversion should be used to build a primer to search for short "
                    "reads. The longer it is, the more specific the primers will be to survey short reads. But if "
                    "if it is too long, then there will not be enough reads to 'keep' depending on the minimum "
                    "short read length of the sequencing. The default value will probably good for the vast "
                    "majority of cases.", metavar="INT")
    groupF.add_argument('--skip-compute-inversion-activity', default=False, action="store_true", help="This flag "
                    "will help you skip computing inversion activity, which is an extremely costly step, even if "
                    "your `bams-and-profiles-txt` contains `r1` and `r2` entries. It may be a good idea to first "
                    "run the workflow without computing activity, take a look at the CONSENSUS file to make sure "
                    "everything looks alright, and then run the workflow without this flag.")
    groupF.add_argument('--end-primer-search-after-x-hits', default=None, type=int, help="For very large "
                    "datasets, primer search can take a very long time. By setting a small integer here, you can ask "
                    "anvi'o to stop searching primers after a few hits. Once the total number of primer hits reach to "
                    "number for a given sample, anvi'o will stop searching further and continue with the next sample. "
                    "This flag is only good for testing, since it will prematurely end the search without testing "
                    "all primers")
    groupF.add_argument('--min-frequency-to-report', default=1, type=int, help="By default, anvi'o will only report "
                    "primers supported by more than one read from the `r1` and `r2` entry in your `bams-and-profiles-txt`. "
                    "The reason for this filtering is sequencing errors, which  can create very low frequency "
                    "entries in the activity report. Use this flag to change that minimum threshold.")

    groupG = parser.add_argument_group('KEY ALGORITHMIC COMPONENT 05: REPORTING GENOMIC CONTEXT AROUND INVERSIONS',
                    "Once the consensus inversions are computed, anvi'o can go back to contigs on which they are "
                    "found, and for each inversion site report the surrounding genes and their functions as flat "
                    "text files in the output directory so you can actually take quick look at them.")
    groupG.add_argument('--skip-recovering-genomic-context', default=False, action="store_true", help="Of course "
                    "you can skip this step because why should anyone have nice things.")
    groupG.add_argument('--num-genes-to-consider-in-context', default=3, type=int, help="With this parameter you "
                    "can adjust the number of genes anvi'o should consider to characterize the genomic context "
                    "surrounding inversions. If you set nothing here, anvi'o will go 3 genes upstream from the "
                    "first palindrome of the inversion and 3 genes downstream from the second. If there are not "
                    "enough genes in the contig given the position of either of the palindromes, anvio' will "
                    "report whatever it can.", metavar="INT")
    groupG.add_argument(*anvio.A('gene-caller'), **anvio.K('gene-caller', {'help': "The gene caller to show gene calls "
                    "from."}))

    groupM = parser.add_argument_group('KEY ALGORITHMIC COMPONENT 06: FINDING MOTIFS IN INVERTED REPEATS',
                    "Site-specific inversions are carried out by site-specific recombinases which recognize "
                    "imperfect palindromic sequences on each end of an inversion. They appear as a single pair "
                    "of inverted repeats, the one anvi'o report in the consensus output. We can use DNA motif search "
                    "to connect inversions that are potentially inverted by the same site-specific recombinase.")
    groupM.add_argument('--skip-search-for-motifs', default=False, action="store_true", help="Skip the search for "
                    "conserved motifs in inverted repeats.")
    groupM.add_argument('--num-of-motifs', default=None, type=int, help="By default anvi'o will search for as many motifs as "
                    "inversions. It can be time consuming if you have a lot of inversions. You can use that flag to "
                    "search for less motifs.", metavar="INT")

    groupH = parser.add_argument_group('TARGETED SEARCH & PROCESSING', "By default anvi'o process every region of every "
                    "contig it is given. Using the following three parameters, you can limit your search to a particular "
                    "contig, start position, or end position. You can use ANY COMBINATION of these three. For instance, "
                    "if you focus on a particular contig only, all regions identified in it by FWD/FWD or REV/REV reads. "
                    "You can limit the search space on that contig using target region start and end parameters.")
    groupH.add_argument(*anvio.A('verbose'), **anvio.K('verbose'))
    groupH.add_argument(*anvio.A('target-contig'), **anvio.K('target-contig'))
    groupH.add_argument(*anvio.A('target-region-start'), **anvio.K('target-region-start'))
    groupH.add_argument(*anvio.A('target-region-end'), **anvio.K('target-region-end'))

    groupH = parser.add_argument_group('OUTPUT DIRECTORY', "Where to put all the output files.")
    groupH.add_argument(*anvio.A('output-dir'), **anvio.K('output-dir'))

    return parser.get_args(parser)


if __name__ == '__main__':
    main()

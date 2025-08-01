#!/usr/bin/env python
# -*- coding: utf-8

import sys
import anvio.profiler as profiler

import anvio

from anvio.terminal import time_program
from anvio.argparse import ArgumentParser
from anvio.errors import ConfigError, FilesNPathsError


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['meren', 'ekiefl', 'ozcan']
__tags__ = ["metagenomics", "profile_db", "contigs_db", "bam", "variability", "clustering"]
__resources__ = [("The usage of the profiler in metagenomics workflow", "http://merenlab.org/2016/06/22/anvio-tutorial-v2/#anvi-profile")]
__requires__ = ['bam-file', 'contigs-db']
__provides__ = ['single-profile-db', 'misc-data-items-order', 'variability-profile']
__description__ = ("The flagship anvi'o program to profile a BAM file. Running this program on a BAM file "
                   "will quantify coverages per nucleotide position in read recruitment results and will "
                   "average coverage and detection data per contig. It will also calculate single-nucleotide, "
                   "single-codon, and single-amino acid variants, as well as structural variants, such as "
                   "insertion and deletions, to eventually stores all data into a single anvi'o profile database. "
                   "For very large projects, this program can demand a lot of time, memory, and storage "
                   "resources. If all you want is to learn coverages of your nutleotides, genes, contigs, or your "
                   "bins collections from BAM files very rapidly, and/or you do not need anvi'o single profile "
                   "databases for your project, please see other anvi'o programs that profile BAM files, "
                   "`anvi-script-get-coverage-from-bam` and `anvi-profile-blitz`")
__resources__ = [("Another description as part of the metagenomic workflow", "http://merenlab.org/2016/06/22/anvio-tutorial-v2/#anvi-profile")]


def get_args():
    parser = ArgumentParser(__description__)

    #############################################################################################################################
    groupI = parser.add_argument_group('INPUT FILES', 'There are two possible inputs for anvio profiler. You must\
                                                  to declare either of these two.')
    #############################################################################################################################
    groupI.add_argument('-i', '--input-file', metavar = 'INPUT_BAM', default = None,
                        help = 'Sorted and indexed BAM file to analyze. Takes a long time depending on the\
                                length of the file and parameters used for profiling.')
    groupI.add_argument(*anvio.A('contigs-db'), **anvio.K('contigs-db', {'required': False}))
    groupI.add_argument(*anvio.A('blank-profile'), **anvio.K('blank-profile'))

    #############################################################################################################################
    groupS = parser.add_argument_group("PROJECT DETAILS", "Set your sample name for this profile. If you don't decalere ."
                                        "anything here, anvi'o will come up with a sample name based on the BAM file name, and "
                                        "you will not be able to change it later.")
    #############################################################################################################################
    groupS.add_argument(*anvio.A('sample-name'), **anvio.K('sample-name'))
    groupS.add_argument(*anvio.A('description'), **anvio.K('description'))

    #############################################################################################################################
    groupQ = parser.add_argument_group("CONTIG FILTERS", "Of all reference contigs that are used for read recruitment, you can "
                                        "focus only those that match certain lenght criteria or a subset you defined explicilty.")
    #############################################################################################################################
    groupQ.add_argument(*anvio.A('min-contig-length'), **anvio.K('min-contig-length'))
    groupQ.add_argument(*anvio.A('max-contig-length'), **anvio.K('max-contig-length'))
    groupQ.add_argument(*anvio.A('contigs-of-interest'), **anvio.K('contigs-of-interest'))
    groupQ.add_argument(*anvio.A('list-contigs'), **anvio.K('list-contigs'))

    #############################################################################################################################
    groupJ = parser.add_argument_group("SHORT READ FILTERS", "Choose which reads to work (or not to work) with, like a pro. "
                                        "These parameters will define which short reads from a BAM file will be considered "
                                        "during anvi'o profiling.")
    #############################################################################################################################
    groupJ.add_argument(*anvio.A('min-percent-identity'), **anvio.K('min-percent-identity', {'help': 
                            'Ignore any reads with a percent identity to the reference less '
                            'than this number, e.g. 95. If not provided, all reads in the BAM '
                            'file will be used (and things will run faster).', 'default': None}))
    groupJ.add_argument(*anvio.A('fetch-filter'), **anvio.K('fetch-filter'))

    #############################################################################################################################
    groupK = parser.add_argument_group("POPULATION GENETICS", "All things related to downstream population genetics analyses. "
                                        "Anvi'o is able to profile single-nucleotide variants and single-codon variants that "
                                        "occur in reads stored in your BAM file. Using these parameters you can have more "
                                        "control over how it is done.")
    #############################################################################################################################
    groupK.add_argument(*anvio.A('report-variability-full'), **anvio.K('report-variability-full'))
    groupK.add_argument(*anvio.A('min-coverage-for-variability'), **anvio.K('min-coverage-for-variability'))
    groupK.add_argument(*anvio.A('skip-SNV-profiling'), **anvio.K('skip-SNV-profiling'))
    groupK.add_argument(*anvio.A('skip-INDEL-profiling'), **anvio.K('skip-INDEL-profiling'))
    groupK.add_argument(*anvio.A('profile-SCVs'), **anvio.K('profile-SCVs'))
    groupK.add_argument(*anvio.A('skip-edges'), **anvio.K('skip-edges'))

    #############################################################################################################################
    groupX = parser.add_argument_group('HIERARCHICAL CLUSTERING', "Do you want your splits to be clustered? Yes? No?\
                                        Maybe? Remember: By default, anvi-profile will not perform hierarchical clustering\
                                        on your splits; but if you use `--blank` flag, it will try. You can skip that by\
                                        using the `--skip-hierarchical-clustering` flag.")
    #############################################################################################################################
    groupX.add_argument(*anvio.A('cluster-contigs'), **anvio.K('cluster-contigs'))
    groupX.add_argument(*anvio.A('skip-hierarchical-clustering'), **anvio.K('skip-hierarchical-clustering'))
    groupX.add_argument(*anvio.A('distance'), **anvio.K('distance', {'help':
                            'The distance metric for the hierarchical clustering. Only relevant if you are\
                             using `--cluster-contigs` flag. The default is "%(default)s".'}))
    groupX.add_argument(*anvio.A('linkage'), **anvio.K('linkage', {'help':
                            'The linkage method for the hierarchical clustering. Just like the distance metric\
                             this is only relevant if you are using it with `--cluster-contigs` flag. The\
                             default is "%(default)s".'}))

    #############################################################################################################################
    groupM = parser.add_argument_group('OUTPUT CONTROL', 'Assert your dominance over output dsetinations.')
    #############################################################################################################################
    groupM.add_argument(*anvio.A('output-dir'), **anvio.K('output-dir'))
    groupM.add_argument(*anvio.A('overwrite-output-destinations'), **anvio.K('overwrite-output-destinations'))

    #############################################################################################################################
    groupZ = parser.add_argument_group('PERFORMANCE', 'Performance settings for profiler')
    #############################################################################################################################
    groupZ.add_argument(*anvio.A('num-threads'), **anvio.K('num-threads'))
    groupZ.add_argument(*anvio.A('queue-size'), **anvio.K('queue-size'))
    groupZ.add_argument(*anvio.A('write-buffer-size-per-thread'), **anvio.K('write-buffer-size-per-thread'))
    groupZ.add_argument('--force-multi', action='store_true',
                        help="This is not useful to non-developers. It forces the multi-process "
                             "routine even when 1 thread is chosen.")

    return parser.get_args(parser)


@time_program
def main():
    args = get_args()

    try:
        profiler.BAMProfiler(args)._run()
    except ConfigError as e:
        print(e)
        sys.exit(-1)
    except FilesNPathsError as e:
        print(e)
        sys.exit(-1)

if __name__ == '__main__':
    main()

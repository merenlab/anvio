#!/usr/bin/env python
# -*- coding: utf-8
"""A DIAMOND and MCL-based pangenome workflow"""

import sys

import anvio
import anvio.panops as panops
import anvio.terminal as terminal
import anvio.constants as constants

from anvio.errors import ConfigError, FilesNPathsError, HDF5Error


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['meren']
__requires__ = ['genomes-storage-db', 'gene-clusters-txt']
__provides__ = ['pan-db', 'misc-data-items-order', 'gene-clusters']
__description__ = ("An anvi'o program to compute a pangenome from an anvi'o genome storage")
__resources__ = [("A tutorial on pangenomics", "http://merenlab.org/2016/11/08/pangenomics-v2/"),]


run = terminal.Run()
progress = terminal.Progress()


def main():
    args = get_args()

    if not args.gene_clusters_txt:
        run.warning('If you publish results from this workflow, please do not forget to cite DIAMOND '
                    '(doi:10.1038/nmeth.3176), unless you use it with --use-ncbi-blast flag, and MCL '
                    '(http://micans.org/mcl/ and doi:10.1007/978-1-61779-361-5_15)', lc = 'yellow')

    try:
        pan = panops.Pangenome(args, run, progress)
        pan.process()
    except ConfigError as e:
        print(e)
        sys.exit(-1)
    except FilesNPathsError as e:
        print(e)
        sys.exit(-2)
    except HDF5Error as e:
        print(e)
        sys.exit(-2)


def get_args():
    # some footwork to help with the help menus:
    if panops.additional_param_sets_for_sequence_search['diamond']:
        DIAMOND = f"the following parameters: \"{panops.additional_param_sets_for_sequence_search['diamond']}\""
    else:
        DIAMOND = "nothing as the default additional parameters for DIAMOND is empty."

    if panops.additional_param_sets_for_sequence_search['ncbi_blast']:
        NCBI_BLAST = f"the following parameters: \"{panops.additional_param_sets_for_sequence_search['ncbi_blast']}\""
    else:
        NCBI_BLAST = "nothing, as the default additional parameters for the NCBI BLAST is empty."


    from anvio.argparse import ArgumentParser

    parser = ArgumentParser(description=__description__)

    groupA = parser.add_argument_group('INPUT DATA: GENOMES & GENES', "You first tell anvi'o where your genes are (which is stored in a fancy 'genomes storage' "
                                "file you have ALREADY generated using the program `anvi-genomes-storage` (if this means nothing to you, see the anvi'o "
                                "pangenomics tutorial)). They you tell anvi'o if you wish to limit your analysis to a smaller set of genomes than those "
                                "anvi'o will find in that genomes-storage-db, or whether or not it should exclude partial gene calls.")
    groupA.add_argument(*anvio.A('genomes-storage'), **anvio.K('genomes-storage', {'required': True}))
    groupA.add_argument(*anvio.A('genomes-names'), **anvio.K('genomes-names'))
    groupA.add_argument('--exclude-partial-gene-calls', default = False, action = 'store_true', help = "By default, anvi'o includes all partial\
                                gene calls from the analysis, which, in some cases, may inflate the number of gene clusters identified and\
                                introduce extra heterogeneity within those gene clusters. Using this flag, you can request anvi'o to exclude\
                                partial gene calls from the analysis (whether a gene call is partial or not is an information that comes directly\
                                from the gene caller used to identify genes during the generation of the contigs database).")


    groupB = parser.add_argument_group('DEFAULT GENE CLUSTER FORMATION: SEQUENCE SEARCH', "The first step in this workflow is to perform reciprocal sequence search "
                                "within the gene pool of input genomes. If you don't change anything in this section, things will work fine, and anvi'o will "
                                "use DIAMOND by default to accomplish the search with default parameters anvi'o developers set for you. But please go "
                                "through each parameter anyway to have an idea regarding what is available to you to play with.")
    groupB.add_argument('--use-ncbi-blast', default = False, action = 'store_true', help = "This program uses DIAMOND by default (and you should, too), "
                                "but if you'd like to use the good ol' blastp by the NCBI instead, it can't stop you.")
    groupB.add_argument('--additional-params-for-seq-search', type=str, metavar = "CMD LINE PARAMS", help = "OK. This is very important. While "
                                f"anvi'o has some defaults for whichever approach you choose to use for your sequence search, you can assume full control "
                                f"over what is passed to the search program. Put anything you wish anvi'o to send your search program in double quotes, and "
                                f"they will be passed to the program. If you don't use this parameter, in addition to the additional parameters anvi'o will "
                                f"use to call your search algorithm of preference, anvi'o will pass to DIAMOND {DIAMOND}, and to the NCBI blast {NCBI_BLAST}. "
                                f"If you use this parameter, it will completely overwrite what you see above. This means, if you are about to use DIAMOND "
                                f"and would like to enable sensitive mode for DIAMOND along with the current anvi'o default additional parameter for it, then "
                                f"you should set this parameter like this manually: --additional-params-for-seq-search "
                                f"\"{panops.additional_param_sets_for_sequence_search['diamond']} --sensitive\". DO NOT EVER FORGET THE DOUBLE QUOTES, unless "
                                f"you hate your computer and want to see it melting beforey our eyes.")


    groupB = parser.add_argument_group('DEFAULT GENE CLUSTER FORMATION: RESOLVING NETWORK', "The second step in this workflow is to use the sequence search results "
                                "to determine gene clusters. Which is a straightforward step thanks to the MCL algorithm. The following parameters are there "
                                "to tweak this step, even though default choices will work for the vast majority of cases.")
    groupB.add_argument('--minbit', type = float, default = 0.5, metavar = "MINBIT", help = "The minimum minbit value. The minbit heuristic \
                                provides a mean to set a to eliminate weak matches between two amino acid sequences. We learned it from ITEP \
                                (Benedict MN et al, doi:10.1186/1471-2164-15-8), which is a comprehensive analysis workflow for pangenomes, \
                                and decided to use it in the anvi'o pangenomic workflow, as well. Briefly, If you have two amino acid sequences,\
                                'A' and 'B', the minbit is defined as 'BITSCORE(A, B) / MIN(BITSCORE(A, A), BITSCORE(B, B))'. So the minbit score\
                                between two sequences goes to 1 if they are very similar over the entire length of the 'shorter' amino acid sequence,\
                                and goes to 0 if (1) they match over a very short stretch compared even to the length of the shorter amino acid sequence\
                                or (2) the match between sequence identity is low. The default is %(default)g.")
    groupB.add_argument('--mcl-inflation', type = float, default = 2.0, metavar = "INFLATION", help = "MCL inflation parameter, that defines\
                                the sensitivity of the algorithm during the identification of the gene clusters. More information on this\
                                parameter and it's effect on cluster granularity is here: (http://micans.org/mcl/man/mclfaq.html#faq7.2).\
                                The default is %(default)g.")
    groupB.add_argument('--min-percent-identity', type = float, default = 0.0, metavar = "PERCENT", help = "Minimum percent identity\
                                between the two amino acid sequences for them to have an edge for MCL analysis. This value will be used\
                                to filter hits from Diamond search results. Because percent identity is not a predictor of a good match (since\
                                it does not communicate many other important factors such as the alignment length between the two sequences and\
                                its proportion to the entire length of those involved), we suggest you rely on 'minbit' parameter. But you know\
                                what? Maybe you shouldn't listen to anyone, and experiment on your own! The default is %(default)g percent.")
    groupB.add_argument('--min-occurrence', type = int, default = 1, metavar = 'NUM_OCCURRENCE', help = "Do you not want singletons? You don't?\
                                Well, this parameter will help you get rid of them (along with doubletons, if you want). Anvi'o will remove\
                                gene clusters that occur less than the number you set using this parameter from the analysis. The default\
                                is %(default)d, which means everything will be kept. If you want to remove singletons, set it to 2, if you want to\
                                remove doubletons as well, set it to 3, and so on.")

    groupC = parser.add_argument_group('ALTERNATIVE GENE CLUSTER FORMATION: USER-PROVIDED GENE CLUSTERS', "As an alternative approach, you can provide "
                                "anvi'o with your own gene cluster affiliations. In this case, anvi'o would not use the default way of computing gene "
                                "clusters de novo, but take your word for which genes go together.")
    groupC.add_argument(*anvio.A('gene-clusters-txt'), **anvio.K('gene-clusters-txt'))

    groupD = parser.add_argument_group("GENE CLUSTER CHARACTERIZATION", "Parameters that mostly impact how anvi'o characterizes your gene clusters once "
                                "they are identified. Such as aligning within gene cluster sequences, computing homogeneity indices for them, etc. "
                                "Important stuff.")
    groupD.add_argument('--skip-alignments', default = False, action = 'store_true', help = "By default, anvi'o attempts to align amino acid\
                                sequences in each gene cluster using multiple sequence alignment via muscle. You can use this flag to skip\
                                that step and be upset later.")
    groupD.add_argument(*anvio.A('align-with'), **anvio.K('align-with'))
    groupD.add_argument('--skip-homogeneity', default=False, action='store_true', dest='skip_homogeneity',help="By default, anvi'o attempts to calculate homogeneity\
                                values for every gene cluster, given that they are aligned. You can use this flag to have anvi'o skip\
                                homogeneity calculations. Anvi'o will ignore this flag if you decide to skip alignments")
    groupD.add_argument('--quick-homogeneity', default=False, action='store_true', dest='quick_homogeneity',help="By default, anvi'o will use a homogeneity\
                                algorithm that checks for horizontal and vertical geometric homogeneity (along with functional). With this\
                                flag, you can tell anvi'o to skip horizontal geometric homogeneity calculations. It will be less accurate but quicker.\
                                Anvi'o will ignore this flag if you skip homogeneity calculations or alignments all together.")

    groupE = parser.add_argument_group("OTHERS", "Sweet parameters of convenience.")
    groupE.add_argument(*anvio.A('project-name'), **anvio.K('project-name'))
    groupE.add_argument(*anvio.A('description'), **anvio.K('description'))
    groupE.add_argument(*anvio.A('output-dir'), **anvio.K('output-dir', {'metavar':'PAN_DB_DIR'}))
    groupE.add_argument(*anvio.A('overwrite-output-destinations'), **anvio.K('overwrite-output-destinations'))
    groupE.add_argument(*anvio.A('num-threads'), **anvio.K('num-threads'))

    groupF = parser.add_argument_group("ORGANIZING GENE CLUSTERs", "These are stuff that will change the clustering dendrogram of your gene clusters.")
    groupF.add_argument(*anvio.A('skip-hierarchical-clustering'), **anvio.K('skip-hierarchical-clustering', {'help': "Anvi'o attempts\
                                to generate a hierarchical clustering of your gene clusters once it identifies them so you can use\
                                `anvi-display-pan` to play with it. But if you want to skip this step, this is your flag."}))
    groupF.add_argument(*anvio.A('enforce-hierarchical-clustering'), **anvio.K('enforce-hierarchical-clustering', {'help': "If you\
                                want anvi'o to try to generate a hierarchical clustering of your gene clusters even if the number of gene clusters exceeds\
                                its suggested limit for hierarchical clustering, you can use this flag to enforce it. Are you are a\
                                rebel of some sorts? Or did computers make you upset? Express your anger towards machine using this\
                                flag."}))
    groupF.add_argument(*anvio.A('distance'), **anvio.K('distance', {'default': None, 'help':
                      'The distance metric for the clustering of gene clusters. If you do not use this flag,\
                       the default distance metric will be used for each clustering configuration\
                       which is "%s".' % constants.distance_metric_default}))
    groupF.add_argument(*anvio.A('linkage'), **anvio.K('linkage', {'default': None, 'help':
                      'The same story with the `--distance`, except, the system default for this one\
                       is %s.' % constants.linkage_method_default}))

    return parser.get_args(parser)


if __name__ == '__main__':
    main()

#!/usr/bin/env python
# -*- coding: utf-8
"""A script to remove replicated genomes from a list of internal and external genome databases or fasta files"""

import sys
from anvio.argparse import ArgumentParser

import anvio
import anvio.terminal as terminal

from anvio.errors import ConfigError, FilesNPathsError
from anvio.genomesimilarity import Dereplicate

__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['ekiefl', 'mahmoudyousef98']
__requires__ = ['external-genomes', 'internal-genomes', 'fasta', 'genome-similarity']
__provides__ = ['fasta']
__description__ = ("Identify redundant (highly similar) genomes")


@terminal.time_program
def main():
    args = get_args()

    try:
        derep = Dereplicate(args)
        derep.process()
        derep.report()
        derep.clean()
    except ConfigError as e:
        print(e)
        sys.exit(-1)
    except FilesNPathsError as e:
        print(e)
        sys.exit(-1)


def get_args():
    parser = ArgumentParser(description=__description__)

    groupA = parser.add_argument_group('INPUT OPTIONS', "Tell anvi'o what you want.")
    groupA.add_argument(*anvio.A('internal-genomes'), **anvio.K('internal-genomes'))
    groupA.add_argument(*anvio.A('external-genomes'), **anvio.K('external-genomes'))
    groupA.add_argument(*anvio.A('fasta-text-file'), **anvio.K('fasta-text-file'))

    groupD = parser.add_argument_group('IMPORT RESULTS', "Alternatively, if you have previous ANI or mash similarity\
                                        computations on your genomes, you can import the result directory here to use. Please note\
                                        that file names must remain unchanged for anvi'o to find them")
    groupD.add_argument('--ani-dir', type=str, metavar='PATH', help="You can import the directory created by `anvi-compute-genome-similarity`\
                        if `--program` parameter was set to `fastANI` or `pyANI` and use it for dereplication")
    groupD.add_argument('--mash-dir', type=str, metavar='PATH', help="You can import the directory created by `anvi-compute-genome-similarity`\
                        if `--program` parameter was set to `sourmash` and use it for dereplication")

    groupB = parser.add_argument_group('OUTPUT OPTIONS', "Tell anvi'o where to store your results.")
    groupB.add_argument(*anvio.A('output-dir'), **anvio.K('output-dir', {'required': True}))
    groupB.add_argument('--skip-fasta-report', action='store_true', help='By default, if any sequence source is\
                        provided, FASTA files of non-redundant genomes are reported. With this flag, no FASTA files\
                        are reported.')
    groupB.add_argument('--report-all', action='store_true', help='By default, only FASTA files of non-redundant\
                        genomes are reported, i.e. single representatives from each cluster. With this flag, all\
                        genome FASTAS will be reported.')

    groupC = parser.add_argument_group('Program', "Tell anvi'o which similarity program to run.")
    groupC.add_argument('--program', type=str, default=None, help="Tell anvi'o which program to run to process genome similarity.\
                        For ANI, you can either use pyANI or fastANI. If accuracy is paramount (for example, distinguishing things less\
                        than 1 percent different), or for dealing with genomes < 80 percent similar,\
                        pyANI is what we recommend. However, fastANI is much faster. If you for some reason want to use mash\
                        similarity, you can use sourmash, but its really not intended for genome comparisons.",
                        choices=['pyANI','fastANI','sourmash'])

    group_FASTANI = parser.add_argument_group('fastANI Settings', "Tell anvi'o to tell fastANI what settings to set.\
                                                                   Only if `--program` is set to `fastANI`")
    group_FASTANI.add_argument('--fastani-kmer-size', type=int, default=16, help="Choose a kmer. The default is %(default)s.")
    group_FASTANI.add_argument('--fragment-length', type=int, default=3000, help="Choose a fragment length. The default is %(default)s.")
    group_FASTANI.add_argument('--min-fraction', type=int, default=0.25, help="Minimum fraction of alignment to be shared between genome pairs \
                        to calculate ANI. If reference and query genome size differ, smaller one among the two is considered. The default is %(default)s.")

    groupE = parser.add_argument_group('pyANI Settings', "Tell anvi'o to tell pyANI what method you wish to use and what settings to set.\
                        Only if `--program` is set to `pyANI`")
    groupE.add_argument('--method', default='ANIb', type=str, help="Method for pyANI. The default is %(default)s.\
                         You must have the necessary binary in path for whichever method you choose. According to\
                         the pyANI help for v0.2.7 at https://github.com/widdowquinn/pyani, the method 'ANIm' uses\
                         MUMmer (NUCmer) to align the input sequences. 'ANIb' uses BLASTN+ to align 1020nt fragments\
                         of the input sequences. 'ANIblastall': uses the legacy BLASTN to align 1020nt fragments\
                         Finally, 'TETRA': calculates tetranucleotide frequencies of each input sequence",\
                         choices=['ANIm', 'ANIb', 'ANIblastall', 'TETRA'])
    groupE.add_argument(*anvio.A('min-alignment-fraction'), **anvio.K('min-alignment-fraction', params_dict={'default':0.25}))
    groupE.add_argument(*anvio.A('significant-alignment-length'), **anvio.K('significant-alignment-length'))
    groupE.add_argument(*anvio.A('use-full-percent-identity'), **anvio.K('use-full-percent-identity'))
    groupE.add_argument(*anvio.A('min-full-percent-identity'), **anvio.K('min-full-percent-identity'))

    groupF = parser.add_argument_group('sourmash settings', "Tell anvi'o to run sourmash with specific settings. Only\
                        if `--program` is set to `sourmash`")
    groupF.add_argument('--kmer-size', type=int, default=13, metavar='INT', help="Set the k-mer size for mash\
                        similarity checks. The default is %(default)s.")
    groupF.add_argument('--scale', type=int, default=1000, metavar='INT', help='Set the compression ratio for\
                        fasta signature file computations. The default is 1000. Smaller ratios decrease sensitivity,\
                        while larger ratios will lead to large fasta signatures.')

    groupG = parser.add_argument_group('Dereplication Parameters', "Some parameters to guide your dereplication")
    groupG.add_argument('--similarity-threshold', type=float, required=True, help="If two genomes have a similarity\
                        greater than or equal to this threshold, they will belong to the same cluster.\
                        Since measures of 'similarity' depend strongly on what method is used for calculation, and\
                        since the threshold at which two genomes should be considered 'similar enough' to be considered\
                        redundant will depend on the application, anvi'o refuses to provide a default parameter.\
                        If you're using pyANI, maybe 0.90 is what you're after. If you're using sourmash, maybe 0.25\
                        is what you're after. Or maybe not? Anvi'o is feeling nervous about this decision.")
    groupG.add_argument('--cluster-method', type=str, default='simple_greedy', help="Currently, genomes are clustered\
                        based on a simple greedy algorithm. Let's say your similarity threshold is 0.90. If genome A\
                        is 0.95 similar to B, and B is 0.95 similar to C, and C is 0.95 similar to D, then {A,B,C,D}\
                        will form a cluster. This is *even though* D may share a similarity to A of merely 0.80,\
                        which is below similarity threshold. You want better alternatives? Contact the developers with\
                        your ideas.",
                        choices=['simple_greedy'])
    groupG.add_argument('--representative-method', default='centrality', type=str, help="After genomes are grouped into\
                        redundancy clusters, you can define how anvi'o picks the representative genome from the\
                        cluster. 'Qscore' computes the genome with the highest completion and lowest redundancy as\
                        the representative. 'length' returns the longest genome. 'centrality' returns the genome with\
                        the highest average similarity to everything in the cluster, i.e. the most central. The\
                        default is %(default)s", choices=['Qscore', 'length', 'centrality'])

    groupH = parser.add_argument_group('OTHER IMPORTANT STUFF', "Yes. You're almost done.")
    groupH.add_argument(*anvio.A('num-threads'), **anvio.K('num-threads'))
    groupH.add_argument(*anvio.A('just-do-it'), **anvio.K('just-do-it'))
    groupH.add_argument(*anvio.A('skip-checking-genome-hashes'), **anvio.K('skip-checking-genome-hashes'))
    groupH.add_argument(*anvio.A('log-file'), **anvio.K('log-file'))

    return parser.get_args(parser)


if __name__ == '__main__':
    main()

#!/usr/bin/env python

import sys

import anvio
import anvio.terminal as terminal

from anvio.argparse import ArgumentParser
from anvio.genomereorientation import GenomeReorienter
from anvio.errors import ConfigError, FilesNPathsError

__copyright__ = "Copyleft 2015-2025, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['meren', 'ahenoch']
__requires__ = ['fasta-txt']
__provides__ = ['fasta']
__description__ = ("Reorient circular genomes and scaffold fragmented genomes in a fasta-txt so their coordinates "
                   "match a chosen reference genome using minimap2 alignment and seqkit manipulation.")


@terminal.time_program
def main():
    args = get_args()

    try:
        reorienter = GenomeReorienter(args)
        reorienter.process()

    except ConfigError as e:
        print(e)
        sys.exit(-1)
    except FilesNPathsError as e:
        print(e)
        sys.exit(-2)


def get_args():
    parser = ArgumentParser(description=__description__)

    groupA = parser.add_argument_group('INPUT FILES')
    groupA.add_argument('-f', '--fasta-txt', required=True,
                        help="Two-column TAB-delimited file with genome name and FASTA file path.")
    groupA.add_argument(*anvio.A('output-dir'), **anvio.K('output-dir', {'required': True}))

    groupRef = parser.add_argument_group('SELECTION & PROCESSING OF REFERENCE',
                                         "IMPORTANT: the reference MUST be a single contig, since orienting and ordering "
                                         "contigs requires a single continuous coordinate system to work against, and "
                                         "coordinates from different contigs of a fragmented reference do not form one "
                                         "axis. The program will STOP WITH AN ERROR otherwise. Please note that this is a "
                                         "separate matter from circularity: the reference additionally needs to be a "
                                         "COMPLETE, CIRCULAR sequence if anvi'o is going to ROTATE it, which happens when "
                                         "the reference is auto-selected, or when --use-dnaa-for-reference-orientation is "
                                         "used, since rotating a linear sequence or a fragment of a chromosome is "
                                         "biologically meaningless. If your reference is a genomic locus rather than a "
                                         "complete chromosome, name it explicitly with --reference, which never rotates "
                                         "anything. If you do not set --reference, the program will auto-select a reference "
                                         "by choosing the genome with the fewest contigs (ties are broken by the largest "
                                         "total length), and it will refuse to continue if even that genome turns out to "
                                         "have more than one contig. The only way to work with a fragmented reference is "
                                         "--use-auto-reference-as-is, which skips all steps of rotation over the reference, "
                                         "and uses it as is.")
    groupRef.add_argument('--reference', required=False,
                          help="Name of the entry in fasta-txt to use as the reference orientation. It must be a single "
                               "contig, but it does not have to be a complete circular genome: a single contig covering a "
                               "genomic locus of interest works just as well, since anvi'o never rotates a reference you "
                               "name explicitly here. If omitted, auto-selection applies.")
    groupRef.add_argument('--use-auto-reference-as-is', action='store_true',
                          help="When anvi'o selects the reference genome automatically (i.e., when `--reference` parameter is not "
                               "used) it first chooses the genome with the fewest contigs and then longest among those that have "
                               "the same number of contigs, and THEN, it rotates this auto-picked reference to a reasonable start "
                               "point by taking into consideration all the genomes. Using this flag, you can ask anvi'o to not "
                               "tinker with the auto-picked reference orientation. This is most useful when the collection of "
                               "genomes (or contigs) all start with an evolutionarily meaningful positions, and all you want to "
                               "do is to reverse complement those that need it so all genomes (or contigs) have the same "
                               "orientation. This is ALSO the only supported way to move forward when none of the genomes in "
                               "your fasta-txt file is a single contig -- but please be aware of what you are giving up: since "
                               "coordinates that come from different contigs of a fragmented reference do not form a single "
                               "continuous axis, the contig ordering, the scaffolding output, and the reported start positions "
                               "will not be trustworthy. This flag is obviously not compatible with `--reference` and "
                               "`--use-dnaa-for-reference-orientation`.")
    groupRef.add_argument('--use-dnaa-for-reference-orientation', action='store_true',
                          help="Use DnaA gene location to orient the reference genome. The program will identify the DnaA "
                               "gene using an HMM profile (Bac_DnaA_C from Pfam), and rotate the reference to start near "
                               "the DnaA gene, which typically marks the origin of replication in bacterial genomes. Since "
                               "this rotates the reference, it requires the reference genome to be a single contig. This "
                               "option is useful for bacterial genomes but may not work well for plasmids or viral genomes "
                               "without DnaA, or fragments from genomes such as specific genomic loci you are interested in "
                               "(for such cases you should certainly consider naming your reference explicitly with "
                               "`--reference`, which never rotates it, or using `--use-auto-reference-as-is`).")

    groupScaffold = parser.add_argument_group('SCAFFOLDING OPTIONS',
                                              "These options control how multi-contig (fragmented) genomes are ordered "
                                              "and oriented relative to the reference. Single-contig genomes are always "
                                              "processed using the circular reorientation workflow.")
    groupScaffold.add_argument('--scaffold-fragmented', action='store_true',
                               help="Insert N characters between ordered contigs to represent missing sequence. "
                                    "The number of Ns inserted equals the gap size on the reference (distance between "
                                    "consecutive contig alignments minus any overlap). Without this flag, contigs are "
                                    "simply concatenated in order without gap padding.")
    groupScaffold.add_argument('--min-contig-length', type=int, default=1000,
                               help="Minimum contig length (in bp) to include in scaffolding. Contigs shorter than "
                                    "this threshold are excluded from alignment and will not appear in output. This "
                                    "threshold also serves as the smallest fragment anvi'o is willing to create when it "
                                    "cuts a contig that runs past the ends of the reference, which keeps the handful of "
                                    "nucleotides `minimap2` routinely soft-clips off the end of a divergent alignment "
                                    "from being promoted into contigs of their own. Default: %(default)s bp.")
    groupScaffold.add_argument('--keep-query-contigs-intact', action='store_true',
                               help="Never cut or rotate a query contig, even when keeping it in one piece misrepresents how it "
                                    "aligns to the reference. By default anvi'o keeps your contigs intact whenever it "
                                    "can, but it will cut one into fragments when the contig runs past the beginning or "
                                    "the end of the reference. This is becasue such a contig will have no single position "
                                    "on the reference to be placed at, since one part of it belongs to the very start of "
                                    "the reference coordinates and another part of it to the very end (which is exactly what happens "
                                    "when your reference and your query genome were circularized at different "
                                    "positions), or because one of its ends hangs off the reference with too little "
                                    "reference left over for that sequence to sit on. Anvi'o cuts such contigs, aligns each "
                                    "fragment on its own merit, and tells you about it. Which means the output FASTA file for a "
                                    "genome may have MORE contigs than the input FASTA did! If your downstream analyses cannot "
                                    "tolerate that, or you are doing an analysis that requires your contigs to be kept intact, "
                                    "then you can use this flag if you promise that you will keep in mind the tradeoff: your "
                                    "contigs will survive in one piece, and in exchange parts of them will not actually be "
                                    "aligned to the reference even though the output file will look like they are. This flag "
                                    "will ALSO suppress the rotation of circularly permuted contigs (contigs an assembler cut "
                                    "out of a cycle in its assembly graph at an arbitrary point, which anvi'o will recognize "
                                    "because their two ends are neighbours on the reference and will rotates back into co-linearity). "
                                    "The cost of preventing this step is that a large part of such a contig will sit in the "
                                    "wrong place. Long story short, anvi'o will report every contig as they were in the query genome "
                                    "apart from reverse-complementing them when necessary, and re-ordering them to match the order "
                                    "of the reference.")

    groupViz = parser.add_argument_group('VISUALIZATION',
                                         "By default, the program will generate synteny ribbon plots showing "
                                         "alignment patterns between the reference and each query genome before "
                                         "and after reorientation/scaffolding. These visualizations help assess "
                                         "the quality of the reorientation process.")
    groupViz.add_argument('--skip-visualizing-alignments', action='store_true',
                          help="Skip generating alignment visualizations. This can speed up processing when "
                               "you only need the reoriented FASTA files without visual inspection.")
    groupViz.add_argument('--plot-width', type=int,
                          help="Width of alignment plots in characters. By default, plots automatically scale "
                               "to use the full terminal width. Use this to override with a specific width. "
                               "Note: widths below 100 characters may not display properly.")
    groupViz.add_argument('--plot-height', type=int, default=20,
                          help="Height of alignment plots in characters. Default: %(default)s")

    groupB = parser.add_argument_group('ADVANCED')
    groupB.add_argument('--threads', type=int, default=1,
                        help="Number of threads for minimap2.")
    groupB.add_argument('--log-file-path',
                        help="Write a detailed log to this file (otherwise only concise reporting is printed).")

    return parser.get_args(parser)


if __name__ == '__main__':
    main()

#!/usr/bin/env python
"""A foldseek and MCL-based structure-informed pangenome workflow"""

import sys

import anvio
import anvio.panops as panops
import anvio.terminal as terminal

from anvio.errors import ConfigError, FilesNPathsError, HDF5Error


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['metehaansever', 'FlorianTrigodet', 'meren']
__requires__ = ['pan-db', 'genomes-storage-db']
__can_use__ = ['structure-db', 'structures-txt']
__provides__ = ['pan-db']
__description__ = ("An anvi'o program to merge the gene clusters of a pangenome into protein structure-informed gene "
                   "clusters, based on the similarity of their predicted structures")
__resources__ = []


run = terminal.Run()
progress = terminal.Progress()


def main():
    args = get_args()

    run.warning('If you publish results from this workflow, please do not forget to cite foldseek '
                '(doi:10.1038/s41587-023-01773-0) and MCL (http://micans.org/mcl/ and '
                'doi:10.1007/978-1-61779-361-5_15)', lc='yellow')

    try:
        pan = panops.StructurePangenome(args, run, progress)
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
    from anvio.argparse import ArgumentParser

    parser = ArgumentParser(description=__description__)

    groupA = parser.add_argument_group('INPUT PANGENOME', "The conventional pangenome whose gene clusters will be merged into "
                                "protein structure-informed gene clusters (PSGCs), and the genomes storage it was computed from. "
                                "The genomes are always the ones of this pangenome.")
    groupA.add_argument(*anvio.A('pan-db'), **anvio.K('pan-db'))
    groupA.add_argument(*anvio.A('genomes-storage'), **anvio.K('genomes-storage', {'required': True}))

    groupB = parser.add_argument_group('STRUCTURES', "Where the structures come from. Each gene cluster is represented by a "
                                "single structure, and you must pick exactly one of these sources. A structure-db generated with "
                                "`anvi-gen-structure-database --pan-db --select-representative` is the way to go: its protein "
                                "for a given gene cluster is the representative of that gene cluster. If a structure-db holds "
                                "more than one protein with a structure for a gene cluster, anvi'o uses the representative of "
                                "that gene cluster in your pangenome, and stops if it is not among them.")
    groupB.add_argument(*anvio.A('structure-db'), **anvio.K('structure-db', {'required': False}))
    groupB.add_argument(*anvio.A('structures-txt'), **anvio.K('structures-txt'))
    groupB.add_argument(*anvio.A('foldseek-search-results'), **anvio.K('foldseek-search-results'))

    groupC = parser.add_argument_group('STRUCTURE-INFORMED GENE CLUSTER FORMATION', "Foldseek compares the structures of all gene "
                                "clusters with each other, and MCL groups gene clusters with similar structures into PSGCs.")
    groupC.add_argument(*anvio.A('min-tm-score'), **anvio.K('min-tm-score'))
    groupC.add_argument(*anvio.A('structure-mcl-inflation'), **anvio.K('structure-mcl-inflation'))
    groupC.add_argument(*anvio.A('keep-gcs-without-structure'), **anvio.K('keep-gcs-without-structure'))

    groupD = parser.add_argument_group('OUTPUT DATA', "If you skip this section, anvi'o will put the new pan-db next to the input "
                                "pan-db, named after the project name of the input pangenome.")
    groupD.add_argument(*anvio.A('output-file'), **anvio.K('output-file', {'metavar': 'STRUCTURE_PAN_DB',
                                'help': "File path for the new structure-informed pan-db. It must end with '-STRUCTURE-PAN.db'."}))
    groupD.add_argument(*anvio.A('project-name'), **anvio.K('project-name'))
    groupD.add_argument(*anvio.A('description'), **anvio.K('description'))
    groupD.add_argument(*anvio.A('intermediate-data-dir'), **anvio.K('intermediate-data-dir'))
    groupD.add_argument(*anvio.A('overwrite-output-destinations'), **anvio.K('overwrite-output-destinations'))

    groupE = parser.add_argument_group("PSGC CHARACTERIZATION", "Aligning the sequences within each PSGC and computing homogeneity "
                                "indices for them.")
    groupE.add_argument('--skip-alignments', default = False, action = 'store_true', help = "By default, anvi'o aligns the amino "
                                "acid sequences in each PSGC. You can use this flag to skip that step.")
    groupE.add_argument(*anvio.A('align-with'), **anvio.K('align-with'))
    groupE.add_argument('--skip-homogeneity', default=False, action='store_true', dest='skip_homogeneity', help="Skip the "
                                "computation of homogeneity indices for PSGCs. Anvi'o ignores this flag if you skip alignments.")
    groupE.add_argument('--quick-homogeneity', default=False, action='store_true', dest='quick_homogeneity', help="Skip the "
                                "horizontal geometric homogeneity calculations. Less accurate but quicker. Anvi'o ignores this flag "
                                "if you skip homogeneity calculations or alignments.")

    groupF = parser.add_argument_group("ORGANIZING PSGCs", "These change the clustering dendrogram of your PSGCs.")
    groupF.add_argument(*anvio.A('skip-hierarchical-clustering'), **anvio.K('skip-hierarchical-clustering'))
    groupF.add_argument(*anvio.A('enforce-hierarchical-clustering'), **anvio.K('enforce-hierarchical-clustering'))
    groupF.add_argument(*anvio.A('distance'), **anvio.K('distance', {'default': None}))
    groupF.add_argument(*anvio.A('linkage'), **anvio.K('linkage', {'default': None}))

    groupG = parser.add_argument_group("OTHERS", "Sweet parameters of convenience.")
    groupG.add_argument(*anvio.A('compare-with-input-pan'), **anvio.K('compare-with-input-pan'))
    groupG.add_argument(*anvio.A('num-threads'), **anvio.K('num-threads'))

    return parser.get_args(parser)


if __name__ == '__main__':
    main()

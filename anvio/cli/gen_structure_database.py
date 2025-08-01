#!/usr/bin/env python
# -*- coding: utf-8

import sys
import anvio
from anvio.argparse import ArgumentParser

import anvio.terminal as terminal
import anvio.structureops as structops

from anvio.errors import ConfigError, FilesNPathsError, ModellerError


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['ekiefl']
__resources__ = [("A conceptual tutorial on the structural biology capabilities of anvio",
                  "http://merenlab.org/2018/09/04/structural-biology-with-anvio/"),
                 ("A practical tutorial section in the infant gut tutorial",
                  "http://merenlab.org/tutorials/infant-gut/#chapter-vii-linking-genomic-heterogeneity-to-protein-structures")]
__requires__ = ['contigs-db', 'pdb-db']
__provides__ = ['structure-db']
__description__ = ("Creates a database of protein structures. Predict protein structures using "
                   "template-based homology modelling of genes in your contigs database, or import "
                   "pre-computed PDB structures you already have.")


@terminal.time_program
def main():
    args = get_args()
    args.structure_db = args.output_db_path

    try:
        structops.StructureSuperclass(args, create=True)._run()
    except ConfigError as e:
        print(e)
        sys.exit(1)
    except FilesNPathsError as e:
        print(e)
        sys.exit(2)
    except ModellerError as e:
        print(e)
        sys.exit(3)


def get_args():
    parser = ArgumentParser(description=__description__)

    groupD = parser.add_argument_group('DATABASES', 'Declaring relevant anvi\'o databases. First things first.')
    groupD.add_argument(*anvio.A('contigs-db'), **anvio.K('contigs-db'))
    groupD.add_argument("--pdb-db", type=str, default=None, help = \
                        """By default, this program accesses the structure files it needs from an
                        internal anvi'o database that can be set up with anvi-setup-pdb-database. If
                        a required structure is not in this database, it will instead be
                        downloaded from the RCSB PDB server. This parameter exists only if a) you
                        created a database and b) it exists in a custom location. In this case,
                        please provide that path here. Otherwise we vibing.""")

    groupG = parser.add_argument_group('GENES', 'Specifying which genes you want to be modelled.')
    groupG.add_argument(*anvio.A('genes-of-interest'), **anvio.K('genes-of-interest'))
    groupG.add_argument(*anvio.A('gene-caller-ids'), **anvio.K('gene-caller-ids'))

    groupO = parser.add_argument_group('OUTPUT', 'Output file and output style.')
    groupO.add_argument(*anvio.A('output-db-path'), **anvio.K('output-db-path'))
    groupO.add_argument(*anvio.A('dump-dir'), **anvio.K('dump-dir'))

    groupI = parser.add_argument_group('NON-MODELLER OPTIONS', 'Alternatives to using MODELLER')
    groupI.add_argument(*anvio.A('external-structures'), **anvio.K('external-structures'))

    groupM = parser.add_argument_group('MODELLER PARAMS', 'Parameters for MODELLER\'s homology modeling.')
    groupM.add_argument("--num-models", "-N", type=int, default=1, help = \
                        """This parameter determines the number of predicted structures that are
                        solved for a given protein.  The original atomic positions for each model
                        are perturbed by an amount defined by --deviation, which leads to
                        differences between each model. Therefore, whichever of the N models is
                        chosen to be the "best" model is more likely to be accurate when
                        --num-models is high, since more of the solution space is searched. It
                        should be kept in mind that the largest determinant of a model's accuracy is
                        determined by the protein templates used, so no need to go overboard with an
                        excessively large --num-models. The default is %(default)d.""")
    groupM.add_argument("--deviation", "-d", type=float, default=4.0, help = \
                        """Deviation (angstroms)""")
    groupM.add_argument("--modeller-database", "-D", type=str, default="pdb_95", help = \
                        """Which database do you want to search the structures of? Default is
                        "%(default)s". If you have your own database it must have either the extension
                        .bin or .pir. If you don't have a database or don't know what this
                        means, don't worry, we will both inform you and take care of you.""")
    groupM.add_argument("--scoring-method", "-b", type=str, default="DOPE_score", help = \
                        """How should the best model be decided? The metric used could be any of
                        GA341_score, DOPE_score, and molpdf. GA341 is an absolute measure,
                        where a good model will have a score near 1.0, whereas anything below 0.6
                        can be considered bad. DOPE and molpdf scores are relative energy measures,
                        where lower scores are better. DOPE has been generally shown to be a better
                        distinguisher between good and bad models than molpdf. The default is %(default)s.
                        To learn more see the MODELLER tutorial:
                        https://salilab.org/modeller/tutorial/basic.html.""")
    groupM.add_argument("--very-fast", action = 'store_true', help = \
                        """If provided, a very fast optimization is done for each model at the cost
                        of accuracy. It is recommended to use a --num-models of 1, since the
                        optimization is so crude that all models will likely converge to the same
                        solution.""")
    groupM.add_argument("--percent-cutoff", "-p", type=float, default=30, help = \
                        """If a protein in the database has a percent identity to the gene of
                        interest that is less than this parameter, then
                        it is not considered as a template. The default is %(default)f.""")
    groupM.add_argument("--alignment-fraction-cutoff", "-a", type=float, default=0.80, help = \
                        """If a protein in the database aligns to a fraction of the gene of interest
                        that is less than this parameter, the template is not considered. For example,
                        if --alignment-cutoff is set to 0.90, and the fraction of the gene of interest
                        that is covered by a potential template is 0.80 in their alignment, the template
                        does not align to enough of the gene of interest to be considered. The default
                        is %(default)f.""")
    groupM.add_argument("--max-number-templates", "-t", type=int, default=5, help = \
                        """Generally speaking it is best to use as many templates as possible given
                        that they have high proper percent identity to the gene of interest. Taken
                        from https://salilab.org/modeller/methenz/andras/node4.html: 'The use of
                        several templates generally increases the model accuracy. One strength of
                        MODELLER is that it can combine information from multiple template
                        structures, in two ways. First, multiple template structures may be aligned
                        with different domains of the target, with little overlap between them, in
                        which case the modeling procedure can construct a homology-based model of
                        the whole target sequence. Second, the template structures may be aligned
                        with the same part of the target, in which case the modeling procedure is
                        likely to automatically build the model on the locally best template
                        [43,44]. In general, it is frequently beneficial to include in the modeling
                        process all the templates that differ substantially from each other, if they
                        share approximately the same overall similarity to the target sequence.' The
                        default is %(default)d.""")

    groupE = parser.add_argument_group('EXTRA', 'Everything else.')
    groupE.add_argument("--skip-DSSP", action = "store_true", help = \
                        """Dictionary of Secondary Structure of Proteins (DSSP) is a program that takes
                        as its input a protein structure file and outputs predicted secondary
                        structure (alpha helix, beta strand, etc.), measures of solvent
                        accessibility, and hydrogen bonds for each residue in the protein. If for
                        some reason you don't want this, provide this flag.""")
    groupE.add_argument("--modeller-executable", type=str, help = \
                        """The MODELLER program to use. For example, `mod9.19`. Anvi'o will try and find
                        it if not provided""")
    groupE.add_argument("--offline-mode", action = "store_true", help = \
                        """Anvi'o first tries to obtain template structures from a database (see
                        --pdb-db for details). If the requested template does not exist in the
                        database, its structure will be downloaded from the RCSB PDB server.
                        However, if you don't have access to internet, or hate the RCSB PDB,
                        provide this flag so that all operations of this program remain offline. If
                        the template structure is not in the database, then no template structure
                        for you.""")
    groupE.add_argument(*anvio.A('num-threads'), **anvio.K('num-threads'))

    params = {
        'default': 25,
        'help': anvio.K('write-buffer-size-per-thread')['help'] + \
                ' If --num-threads is 1, this parameter is ignored because the DB is written to after each gene'
    }
    groupE.add_argument(*anvio.A('write-buffer-size-per-thread'), **anvio.K('write-buffer-size-per-thread', params))

    return parser.get_args(parser)


if __name__ == '__main__':
    main()

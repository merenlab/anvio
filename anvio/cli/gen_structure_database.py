#!/usr/bin/env python

import sys
import argparse
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
__requires__ = ['contigs-db']
__can_use__ = ['pdb-db', 'genes-of-interest-txt']
__provides__ = ['structure-db']
__description__ = ("Creates a database of protein structures. Predict protein structures for genes in "
                   "your contigs database using either template-based homology modelling (MODELLER) or "
                   "AlphaFold2 (ColabFold), or import pre-computed PDB structures you already have.")


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
    groupD.add_argument("--pdb-db", type=str, default=None, help =
                        """By default, this program accesses the structure files it needs from an
                        internal anvi'o database that can be set up with anvi-setup-pdb-database. If
                        a required structure is not in this database, it will instead be
                        downloaded from the RCSB PDB server. This parameter exists only if a) you
                        created a database and b) it exists in a custom location. In this case,
                        please provide that path here. Otherwise we vibing.""")

    groupG = parser.add_argument_group('GENES', 'Specifying which genes you want structures for.')
    groupG.add_argument(*anvio.A('genes-of-interest'), **anvio.K('genes-of-interest'))
    groupG.add_argument(*anvio.A('gene-caller-ids'), **anvio.K('gene-caller-ids'))

    groupO = parser.add_argument_group('OUTPUT', 'Output file and output style.')
    groupO.add_argument(*anvio.A('output-db-path'), **anvio.K('output-db-path'))
    groupO.add_argument(*anvio.A('dump-dir'), **anvio.K('dump-dir'))

    groupN = parser.add_argument_group('ENGINE', 'Which structure-prediction engine should anvi\'o use?')
    groupN.add_argument("--engine", type=str, default='modeller', choices=['modeller', 'colabfold'], help =
                        """The structure-prediction engine to use. 'modeller' (the default) predicts
                        structures with template-based homology modelling. 'colabfold' predicts structures
                        with AlphaFold2 via ColabFold. Each engine has its own set of parameters (see the
                        MODELLER PARAMS and COLABFOLD PARAMS groups), while a few parameters apply to both
                        (see SHARED PREDICTION PARAMS). This has no effect if you import structures with
                        --external-structures instead of predicting them.""")

    groupS = parser.add_argument_group('SHARED PREDICTION PARAMS', 'Parameters that apply to whichever '
                                       'prediction engine you selected with --engine.')
    groupS.add_argument("--num-models", "-N", type=int, default=1, help =
                        """The number of structures ("models") to compute per protein; the best one is
                        kept. For the MODELLER engine, each of the N models perturbs the initial atomic
                        positions (by an amount set with --deviation) and is optimized independently, so a
                        larger N searches more of the solution space (though template quality matters most,
                        so there is no need to go overboard). For the ColabFold engine, this is the number
                        of AlphaFold2 models (1-5) to run, and the highest-confidence one is kept. In both
                        cases, fewer models is faster but explores less. The default is %(default)d.""")
    groupS.add_argument("--skip-DSSP", action = "store_true", help =
                        """Dictionary of Secondary Structure of Proteins (DSSP) is a program that takes
                        as its input a protein structure file and outputs predicted secondary
                        structure (alpha helix, beta strand, etc.), measures of solvent
                        accessibility, and hydrogen bonds for each residue in the protein. If for
                        some reason you don't want this, provide this flag.""")

    groupM = parser.add_argument_group('MODELLER PARAMS', 'Parameters for MODELLER\'s homology modeling. '
                                       'These only apply when `--engine modeller` (the default) is set.')
    groupM.add_argument("--deviation", "-d", type=float, default=4.0, help =
                        """Deviation (angstroms)""")
    groupM.add_argument("--modeller-database", "-D", type=str, default="pdb_95", help =
                        """Which database do you want to search the structures of? Default is
                        "%(default)s". If you have your own database it must have either the extension
                        .bin or .pir. If you don't have a database or don't know what this
                        means, don't worry, we will both inform you and take care of you.""")
    groupM.add_argument("--scoring-method", "-b", type=str, default="DOPE_score", help =
                        """How should the best model be decided? The metric used could be any of
                        GA341_score, DOPE_score, and molpdf. GA341 is an absolute measure,
                        where a good model will have a score near 1.0, whereas anything below 0.6
                        can be considered bad. DOPE and molpdf scores are relative energy measures,
                        where lower scores are better. DOPE has been generally shown to be a better
                        distinguisher between good and bad models than molpdf. The default is %(default)s.
                        To learn more see the MODELLER tutorial:
                        https://salilab.org/modeller/tutorial/basic.html.""")
    groupM.add_argument("--very-fast", action = 'store_true', help =
                        """If provided, a very fast optimization is done for each model at the cost
                        of accuracy. It is recommended to use a --num-models of 1, since the
                        optimization is so crude that all models will likely converge to the same
                        solution.""")
    groupM.add_argument("--percent-cutoff", "-p", type=float, default=30, help =
                        """If a protein in the database has a percent identity to the gene of
                        interest that is less than this parameter, then
                        it is not considered as a template. The default is %(default)f.""")
    groupM.add_argument("--alignment-fraction-cutoff", "-a", type=float, default=0.80, help =
                        """If a protein in the database aligns to a fraction of the gene of interest
                        that is less than this parameter, the template is not considered. For example,
                        if --alignment-cutoff is set to 0.90, and the fraction of the gene of interest
                        that is covered by a potential template is 0.80 in their alignment, the template
                        does not align to enough of the gene of interest to be considered. The default
                        is %(default)f.""")
    groupM.add_argument("--max-number-templates", "-t", type=int, default=5, help =
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
    groupM.add_argument("--modeller-executable", type=str, help =
                        """The MODELLER program to use. For example, `mod9.19`. Anvi'o will try and find
                        it if not provided""")
    groupM.add_argument("--offline-mode", action = "store_true", help =
                        """Anvi'o first tries to obtain template structures from a database (see
                        --pdb-db for details). If the requested template does not exist in the
                        database, its structure will be downloaded from the RCSB PDB server.
                        However, if you don't have access to internet, or hate the RCSB PDB,
                        provide this flag so that all operations of this program remain offline. If
                        the template structure is not in the database, then no template structure
                        for you.""")

    groupC = parser.add_argument_group('COLABFOLD PARAMS', 'Parameters for structure prediction with ColabFold '
                                       '(AlphaFold2). These only apply when `--engine colabfold` is set. ColabFold '
                                       'is not assumed to be on your $PATH: provide the name of the conda environment '
                                       'it is installed in with --colabfold-conda-env.')
    groupC.add_argument("--colabfold-conda-env", type=str, default=None, metavar='ENV_NAME', help =
                        """The name of the conda environment in which ColabFold (colabfold_batch /
                        colabfold_search) is installed. Anvi'o will run every ColabFold command via
                        `conda run -n ENV_NAME`. If you leave this blank, anvi'o will assume the ColabFold
                        programs are directly available (e.g. on your $PATH).""")
    groupC.add_argument("--colabfold-msa-server", action='store_true', help =
                        """Use the public MMseqs2 MSA server hosted by the ColabFold team to generate the
                        multiple sequence alignment (MSA) step online. This requires an internet connection
                        and is fine for a handful of sequences. Mutually exclusive with --colabfold-db. You
                        must choose one of the two.""")
    groupC.add_argument("--colabfold-db", type=str, default=None, metavar='DB_DIR', help =
                        """Generate the MSA step locally instead of using the public server, by searching a
                        local ColabFold database in this directory (the one you set up with ColabFold's
                        `setup_databases.sh`). Recommended when you have many sequences. Mutually exclusive
                        with --colabfold-msa-server. You must choose one of the two.""")
    groupC.add_argument("--num-recycle", type=int, default=None, metavar='INT', help =
                        """Number of prediction recycles for ColabFold. Increasing this can improve prediction
                        quality at the cost of runtime. If not set, ColabFold's default is used.""")
    groupC.add_argument("--amber", action='store_true', help =
                        """Relax the best predicted structure with OpenMM/Amber. This can improve the quality of
                        side-chains at the cost of longer runtime.""")
    groupC.add_argument("--colabfold-additional-parameters", type=str, default=None, metavar='CMD LINE PARAMS', help =
                        """Additional parameters to pass verbatim to `colabfold_batch`, for anything not exposed
                        as a dedicated flag above. Because the value itself begins with dashes, you must attach
                        it with an equals sign so it is not mistaken for anvi'o's own flags, e.g.
                        --colabfold-additional-parameters="--num-seeds 2 --use-dropout". Anvi'o does not validate
                        these, so use with care.""")
    groupC.add_argument("--only-msa", action='store_true', help =
                        """Run only ColabFold's multiple sequence alignment (MSA) step and stop, writing the MSAs
                        and a small checkpoint manifest into --dump-dir. This creates a resumable checkpoint so
                        you can run the CPU-heavy MSA step and the GPU-heavy prediction step separately (e.g. on
                        different cluster nodes). This only works with a local ColabFold database (--colabfold-db):
                        the public MSA server (--colabfold-msa-server) generates the MSA and predicts the structure
                        in a single step that cannot be split. No structure database is produced in this mode; you
                        finish the job later with --only-predict. Requires --dump-dir.""")
    groupC.add_argument("--only-predict", action='store_true', help =
                        """Skip the MSA step and predict structures from an MSA checkpoint that an earlier
                        --only-msa run wrote to --dump-dir, then build the structure database. You must provide the
                        SAME --contigs-db and genes of interest as the --only-msa run: anvi'o verifies that the
                        sequences match the checkpoint before predicting, and refuses to continue if they do not.
                        Requires --dump-dir (pointing at the --only-msa output) and -o. Mutually exclusive with
                        --only-msa.""")

    groupImport = parser.add_argument_group('IMPORT STRUCTURES', 'Instead of predicting structures, import '
                                            'pre-computed ones. This bypasses the --engine choice entirely.')
    groupImport.add_argument(*anvio.A('external-structures'), **anvio.K('external-structures'))

    groupE = parser.add_argument_group('EXTRA', 'Everything else.')
    groupE.add_argument(*anvio.A('num-threads'), **anvio.K('num-threads'))

    params = {
        'default': 25,
        'help': anvio.K('write-buffer-size')['help'] +
                ' If --num-threads is 1, this parameter is ignored because the DB is written to after each gene'
    }
    groupE.add_argument(*anvio.A('write-buffer-size'), **anvio.K('write-buffer-size', params))
    groupE.add_argument(*anvio.A('write-buffer-size-per-thread'), **anvio.K('write-buffer-size-per-thread', {'help': argparse.SUPPRESS}))

    return parser.get_args(parser)


if __name__ == '__main__':
    main()

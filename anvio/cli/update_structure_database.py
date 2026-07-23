#!/usr/bin/env python

import sys
import argparse
import anvio
from anvio.argparse import ArgumentParser

import anvio.structureops as structops

from anvio.errors import ConfigError, FilesNPathsError, ModellerError


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['ekiefl']
__requires__ = ["contigs-db", "structure-db"]
__can_use__ = ["genes-of-interest-txt"]
__description__ = ("Add or re-run genes from an already existing structure database. All settings used "
                   "to generate your database will be used in this program")


def main():
    args = get_args()

    try:
        structops.StructureSuperclass(args, create=False)._run()
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
    groupG = parser.add_argument_group('GENES', 'Specify which genes you want to be modelled. If a gene already '
                                                'exists in the DB, it will be overwritten if --overwrite is set. '
                                                'Otherwise, an error will be raised.')
    groupO = parser.add_argument_group('OUTPUT', 'Output file and output style.')
    groupM = parser.add_argument_group('MODELLER PARAMS', 'Parameters for MODELLER\'s homology modeling.')
    groupC = parser.add_argument_group('COLABFOLD PARAMS', 'Parameters for structure prediction with ColabFold '
                                       '(AlphaFold2). These only apply when the structure database you are updating '
                                       'was created with the ColabFold engine. The scientific parameters (number of '
                                       'models, recycles, Amber relaxation, etc.) are read back from the database so '
                                       'the update stays consistent with how it was created; the parameters below are '
                                       'machine-specific and so must be provided again here.')
    groupE = parser.add_argument_group('EXTRA', 'Everything else.')

    groupD.add_argument(*anvio.A('contigs-db'), **anvio.K('contigs-db'))
    groupD.add_argument(*anvio.A('structure-db'), **anvio.K('structure-db'))
    groupG.add_argument(*anvio.A('genes-of-interest'), **anvio.K('genes-of-interest'))
    groupG.add_argument(*anvio.A('gene-caller-ids'), **anvio.K('gene-caller-ids'))
    groupG.add_argument(*anvio.A('external-structures'), **anvio.K('external-structures'))
    groupO.add_argument(*anvio.A('dump-dir'), **anvio.K('dump-dir'))
    groupM.add_argument('--list-modeller-params', action='store_true', help='Since you are updating an existing DB, modeller params are set in '
                                                                            'place. You can have this program list them by providing this flag')

    groupC.add_argument("--colabfold-conda-env", type=str, default=None, metavar='ENV_NAME', help =
                        """The name of the conda environment in which ColabFold (colabfold_batch /
                        colabfold_search) is installed. Anvi'o will run every ColabFold command via
                        `conda run -n ENV_NAME`. If you leave this blank, anvi'o will assume the ColabFold
                        programs are directly available (e.g. on your $PATH).""")
    groupC.add_argument("--colabfold-msa-server", action='store_true', help =
                        """Use the public MMseqs2 MSA server hosted by the ColabFold team to generate the
                        multiple sequence alignment (MSA) step online. Provide this if the database you are
                        updating was created this way. Mutually exclusive with --colabfold-db.""")
    groupC.add_argument("--colabfold-db", type=str, default=None, metavar='DB_DIR', help =
                        """Generate the MSA step locally by searching a local ColabFold database in this
                        directory (the one you set up with ColabFold's `setup_databases.sh`). Provide this
                        if the database you are updating was created this way. This path is machine-specific,
                        so it is not stored in the structure database and must be given again here. Mutually
                        exclusive with --colabfold-msa-server.""")
    groupC.add_argument("--only-msa", action='store_true', help =
                        """Run only ColabFold's multiple sequence alignment (MSA) step for the genes you are
                        adding and stop, writing the MSAs and a small checkpoint manifest into --dump-dir. This
                        creates a resumable checkpoint so you can run the CPU-heavy MSA step and the GPU-heavy
                        prediction step separately (e.g. on different cluster nodes). This only works with a
                        local ColabFold database (--colabfold-db). No genes are added to the database in this
                        mode; you finish the job later with --only-predict. Requires --dump-dir.""")
    groupC.add_argument("--only-predict", action='store_true', help =
                        """Skip the MSA step and predict structures from an MSA checkpoint that an earlier
                        --only-msa run wrote to --dump-dir, then add the resulting structures to the database.
                        You must provide the SAME --contigs-db and genes of interest as the --only-msa run.
                        Requires --dump-dir (pointing at the --only-msa output). Mutually exclusive with
                        --only-msa.""")

    groupE.add_argument("--rerun-genes", action='store_true', help =
                        """Supply if you would like to rerun structural modelling for your genes of
                        interest if they are already present in your DB""")
    groupE.add_argument("--modeller-executable", type=str, help =
                        """The MODELLER program to use. For example, `mod9.19`. Anvi'o will try and find
                        it if not provided.""")

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

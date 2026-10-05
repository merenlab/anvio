#!/usr/bin/env python
"""A script to export run ANI on every contig in a FASTA file."""

import os
import sys
import shutil
import math
from anvio.argparse import ArgumentParser

import anvio
import anvio.utils as utils
import anvio.fastalib as fastalib
import anvio.terminal as terminal
import anvio.clustering as clustering
import anvio.filesnpaths as filesnpaths

from anvio.drivers import pyani
from anvio.drivers.pyani_methods import validate_ani_method

from anvio.errors import ConfigError, FilesNPathsError


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['meren']
__requires__ = ['fasta']
__provides__ = ['genome-similarity']
__description__ = "Run ANI between contigs in a single FASTA file"


def main():
    try:
        run_program()
    except ConfigError as e:
        print(e)
        sys.exit(-1)
    except FilesNPathsError as e:
        print(e)
        sys.exit(-1)


def run_program():
    args = get_args()
    run = terminal.Run()

    validate_ani_method(args.method)
    ani_backend = getattr(args, 'ani_backend', 'pyani-plus')

    if ani_backend not in ('legacy', 'pyani-plus'):
        raise ConfigError(f"Unknown ANI backend '{ani_backend}'. Choose 'legacy' or 'pyani-plus'.")
    if getattr(args, 'pyani_plus_program', None) is not None and ani_backend == 'legacy':
        raise ConfigError("--pyani-plus-program can only be used with the pyANI-plus backend.")

    if ani_backend == 'legacy':
        program = pyani.PyANI(args)
    elif ani_backend == 'pyani-plus':
        from anvio.drivers import pyani_plus
        program = pyani_plus.PyANIPlus(args)
    filesnpaths.is_file_fasta_formatted(args.fasta_file)
    filesnpaths.check_output_directory(args.output_dir, ok_if_exists=False)

    temp_dir = filesnpaths.get_temp_directory_path()

    fasta = fastalib.SequenceSource(args.fasta_file)
    num_contigs = 0
    while next(fasta):
        if not utils.check_contig_names(fasta.id, dont_raise=True):
            raise ConfigError("At least one of the deflines in your FASTA File does not comply with the 'simple deflines' "
                              "requirement of anvi'o. You can either use the script `anvi-script-reformat-fasta` to take "
                              "care of this issue, or do it manually, but you must know that anvi'o is very upset for "
                              "making you deal with this. Here is more info on simple deflines %s" %
                                   ('http://merenlab.org/2016/06/22/anvio-tutorial-v2/#take-a-look-at-your-fasta-file'))
        num_contigs += 1
    fasta.reset()

    run.info("Num contigs found for ANI", num_contigs)
    run.info("Temporary output directory", temp_dir)
    run.info("Output directory", os.path.abspath(args.output_dir))

    while next(fasta):
        with open(os.path.join(temp_dir, fasta.id + '.fa'), 'w') as f:
            f.write('>%s\n%s\n' % (fasta.id, fasta.seq))

    results = program.run_command(temp_dir)

    clusterings = {}
    omitted_trees = []
    for report_name in results:
        matrix = results[report_name]
        if ani_backend == 'pyani-plus' and any(
                value is None or not math.isfinite(float(value))
                for row in matrix.values() for value in row.values()):
            omitted_trees.append(report_name)
            continue
        clusterings[report_name] = clustering.get_newick_tree_data_for_dict(matrix,
            linkage=args.linkage, distance=args.distance)

    os.mkdir(args.output_dir)
    for report_name in results:
        output_path_for_report = os.path.join(args.output_dir, args.method + '_' + report_name)

        utils.store_dict_as_TAB_delimited_file(results[report_name], output_path_for_report + '.txt')
        if report_name in clusterings:
            with open(output_path_for_report + '.newick', 'w') as f:
                f.write(clusterings[report_name])

    if omitted_trees:
        run.warning("pyANI-plus reports undefined comparisons in %s. Anvi'o preserved those values in the matrix and "
                    "omitted their Newick trees rather than substitute a distance." % ', '.join(omitted_trees))

    if anvio.DEBUG:
        run.warning("The temp directory, %s, is kept. Please don't forget to clean it up "
                    "later" % temp_dir, header="Debug")
    else:
        run.info_single('Cleaning up the temp directory (you can use `--debug` if you would '
                        'like to keep it for testing purposes)', nl_before=1, nl_after=1)
        shutil.rmtree(temp_dir)


def get_args():
    parser = ArgumentParser(description=__description__)

    parser.add_argument(*anvio.A('fasta-file'), **anvio.K('fasta-file', {'required': True}))
    parser.add_argument(*anvio.A('output-dir'), **anvio.K('output-dir', {'required': True }))
    parser.add_argument(*anvio.A('pan-db'), **anvio.K('pan-db', {'required': False}))
    parser.add_argument(*anvio.A('num-threads'), **anvio.K('num-threads', {
        'help': "CPU allocation for pyANI-plus child processes on Linux via taskset. Default is 1 or ANVIO_THREADS; this controls CPU affinity, not every OS thread. Other backends keep their existing thread behavior."
    }))
    parser.add_argument(*anvio.A('log-file'), **anvio.K('log-file'))
    parser.add_argument('--method', default='ANIb', type=str, help="ANI method (default: %(default)s). Supported methods are ANIb and ANIm. \
                         ANIblastall and TETRA are retired in the Python 3.13 port.")
    parser.add_argument(*anvio.A('ani-backend'), **anvio.K('ani-backend', {'default': 'pyani-plus'}))
    parser.add_argument(*anvio.A('pyani-plus-program'), **anvio.K('pyani-plus-program'))
    parser.add_argument(*anvio.A('distance'), **anvio.K('distance', {'help': 'The distance metric for the hierarchical \
                         clustering. The default is "%(default)s".'}))
    parser.add_argument(*anvio.A('linkage'), **anvio.K('linkage', {'help': 'The linkage method for the hierarchical \
                         clustering. The default is "%(default)s".'}))
    parser.add_argument(*anvio.A('just-do-it'), **anvio.K('just-do-it'))

    return parser.get_args(parser)


if __name__ == '__main__':
    main()

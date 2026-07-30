#!/usr/bin/env python
# -*- coding: utf-8

import sys
import os
import tempfile
import shutil
from concurrent.futures import ProcessPoolExecutor, as_completed
from anvio.argparse import ArgumentParser

import anvio
import anvio.utils as utils
import anvio.dbops as dbops
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths

from anvio.errors import ConfigError, FilesNPathsError


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['meren']
__provides__ = ["coverages-txt", "detection-txt"]
__requires__ = ["profile-db", "contigs-db", "genes-of-interest-txt"]
__description__ = ("Export gene coverage and detection data for all genes associated with "
                   "contigs described in a profile database")


def process_gene_chunk(args_dict, gene_chunk, chunk_index, tmp_dir):
    """
    Worker function: initializes its own ProfileSuperclass for a subset of genes
    and writes coverages + detection to temporary files in tmp_dir.
    Returns the paths to the two temp files produced.
    """
    import types
    import anvio.dbops as dbops

    args = types.SimpleNamespace(**args_dict)

    contigs_db = dbops.ContigsDatabase(args.contigs_db)
    where_clause = f"gene_callers_id in ({','.join(str(g) for g in gene_chunk)})"
    splits_for_chunk = set(
        contigs_db.db.get_single_column_from_table(
            "genes_in_splits", "split", unique=True, where_clause=where_clause
        )
    )
    contigs_db.disconnect()

    args.split_names_of_interest = splits_for_chunk

    profile_db = dbops.ProfileSuperclass(args)

    cov_path = os.path.join(tmp_dir, f'chunk-{chunk_index}-GENE-COVERAGES.txt')
    det_path = os.path.join(tmp_dir, f'chunk-{chunk_index}-GENE-DETECTION.txt')

    cov_file = open(cov_path, 'w')
    det_file = open(det_path, 'w')

    samples = profile_db.p_meta['samples']

    def write_to_file():
        for gene_callers_id in profile_db.gene_level_coverage_stats_dict:
            line_cov = [gene_callers_id]
            line_det = [gene_callers_id]
            for sample_name in samples:
                line_cov.append(profile_db.gene_level_coverage_stats_dict[gene_callers_id][sample_name]['mean_coverage'])
                line_det.append(profile_db.gene_level_coverage_stats_dict[gene_callers_id][sample_name]['detection'])
            cov_file.write('\t'.join(map(str, line_cov)) + '\n')
            det_file.write('\t'.join(map(str, line_det)) + '\n')

        profile_db.gene_level_coverage_stats_dict = {}
        profile_db.split_coverage_values_per_nt_dict = {}

    profile_db.init_gene_level_coverage_stats_dicts(
        callback=write_to_file,
        callback_interval=1000,
        gene_caller_ids_of_interest={int(g) for g in gene_chunk}
    )

    cov_file.close()
    det_file.close()

    return cov_path, det_path, samples


def chunk_list_fixed_size(lst, chunk_size):
    """Yield successive fixed-size chunks from lst."""
    lst = list(lst)
    for start in range(0, len(lst), chunk_size):
        yield lst[start:start + chunk_size]


def main():
    args = get_args()
    run = terminal.Run()
    p = terminal.pluralize

    try:
        if args.gene_caller_id and args.genes_of_interest:
            raise ConfigError("You either should define a single gene caller id, or a list of "
                              "gene callers of interest. Not both")

        # ------------------------------------------------------------------ #
        # Resolve genes of interest
        # ------------------------------------------------------------------ #
        if args.genes_of_interest:
            genes_of_interest = set([
                g.strip()
                for g in utils.get_column_data_from_TAB_delim_file(
                    args.genes_of_interest, column_indices=[0], expected_number_of_fields=1
                )[0] if g
            ])
            user_supplied_genes = True
        elif args.gene_caller_id:
            genes_of_interest = set([str(args.gene_caller_id).strip()])
            user_supplied_genes = True
        else:
            genes_of_interest = None
            user_supplied_genes = False

        # ------------------------------------------------------------------ #
        # If parallelism is requested and no gene filter was given, fetch all
        # gene caller IDs from the contigs DB so we can chunk them.
        # ------------------------------------------------------------------ #
        if not genes_of_interest and args.num_threads > 1:
            run.info_single("No gene filter provided — fetching all gene caller IDs from the "
                            "contigs db to enable parallelism.", mc='green')
            contigs_db = dbops.ContigsDatabase(args.contigs_db)
            genes_of_interest = set(
                str(g) for g in contigs_db.db.get_single_column_from_table(
                    "genes_in_contigs", "gene_callers_id", unique=True
                )
            )
            contigs_db.disconnect()

        # ------------------------------------------------------------------ #
        # Validate splits of interest only when the user supplied a gene filter.
        # When we populated genes_of_interest ourselves from the full contigs DB
        # every gene is guaranteed to have a split, so this check is unnecessary.
        # ------------------------------------------------------------------ #
        if user_supplied_genes and genes_of_interest:
            contigs_db = dbops.ContigsDatabase(args.contigs_db)
            where_clause = f"gene_callers_id in ({','.join(genes_of_interest)})"
            splits_of_interest = set(
                contigs_db.db.get_single_column_from_table(
                    "genes_in_splits", "split", unique=True, where_clause=where_clause
                )
            )
            contigs_db.disconnect()

            if not splits_of_interest:
                raise ConfigError(
                    f"You have provided {p('gene caller id', len(genes_of_interest))} of interest "
                    f"but there were no split names in the contigs matching them :/"
                )

            run.warning(
                f"You have provided {p('gene caller id', len(genes_of_interest))} which matched "
                f"to {p('split', len(splits_of_interest))}. Very good job.", lc="green"
            )

            split_names_in_profile_db = utils.get_all_item_names_from_the_database(args.profile_db)
            splits_only_in_contigs_db = splits_of_interest.difference(split_names_in_profile_db)

            if splits_only_in_contigs_db:
                run.warning(
                    f"Hear this: {p('split name', len(splits_only_in_contigs_db))} that matched "
                    f"to the gene calls you were interested in occurred only in the contigs db "
                    f"but not in the profile database. Anvi'o will remove them from the list of "
                    f"splits of interest."
                )
                splits_of_interest -= splits_only_in_contigs_db

            if not splits_of_interest:
                raise ConfigError(
                    "No split names were left to work with :( Genes you are interested in are "
                    "not occurring on any contig that survived the `anvi-profile` step :/ Sorry!"
                )

            args.split_names_of_interest = splits_of_interest

        # ------------------------------------------------------------------ #
        # Validate output paths up front
        # ------------------------------------------------------------------ #
        gene_coverages_path = args.output_file_prefix + '-GENE-COVERAGES.txt'
        gene_detection_path = args.output_file_prefix + '-GENE-DETECTION.txt'
        filesnpaths.is_output_file_writable(gene_coverages_path)
        filesnpaths.is_output_file_writable(gene_detection_path)

        # ------------------------------------------------------------------ #
        # Dispatch: parallel if num_threads > 1, otherwise original single path
        # ------------------------------------------------------------------ #
        if args.num_threads > 1:
            gene_chunks = list(chunk_list_fixed_size(genes_of_interest, args.chunk_size))
            num_chunks = len(gene_chunks)

            run.info_single(
                f"Processing {len(genes_of_interest)} gene(s) across {num_chunks} chunk(s) of "
                f"up to {args.chunk_size} gene(s) each, using up to {args.num_threads} "
                f"parallel worker(s).",
                mc='green'
            )

            _run_parallel(args, gene_chunks, args.num_threads,
                          gene_coverages_path, gene_detection_path, run)
        else:
            _run_single(args, genes_of_interest, gene_coverages_path, gene_detection_path)

        run.info('Gene coverages', gene_coverages_path)
        run.info('Gene detection', gene_detection_path)

    except ConfigError as e:
        print(e)
        sys.exit(-1)
    except FilesNPathsError as e:
        print(e)
        sys.exit(-1)


def _run_single(args, genes_of_interest, gene_coverages_path, gene_detection_path):
    """Original single-process export path."""
    profile_db = dbops.ProfileSuperclass(args)

    cov_file = open(gene_coverages_path, 'w')
    det_file = open(gene_detection_path, 'w')

    header = ['gene_callers_id'] + profile_db.p_meta['samples']
    cov_file.write('\t'.join(header) + '\n')
    det_file.write('\t'.join(header) + '\n')

    def write_to_file():
        for gene_callers_id in profile_db.gene_level_coverage_stats_dict:
            line_cov = [gene_callers_id]
            line_det = [gene_callers_id]
            for sample_name in profile_db.p_meta['samples']:
                line_cov.append(profile_db.gene_level_coverage_stats_dict[gene_callers_id][sample_name]['mean_coverage'])
                line_det.append(profile_db.gene_level_coverage_stats_dict[gene_callers_id][sample_name]['detection'])
            cov_file.write('\t'.join(map(str, line_cov)) + '\n')
            det_file.write('\t'.join(map(str, line_det)) + '\n')

        profile_db.gene_level_coverage_stats_dict = {}
        profile_db.split_coverage_values_per_nt_dict = {}

    if genes_of_interest:
        profile_db.init_gene_level_coverage_stats_dicts(
            callback=write_to_file,
            callback_interval=1000,
            gene_caller_ids_of_interest={int(g) for g in genes_of_interest}
        )
    else:
        profile_db.init_gene_level_coverage_stats_dicts(
            callback=write_to_file,
            callback_interval=1000
        )

    cov_file.close()
    det_file.close()


def _run_parallel(args, gene_chunks, num_workers, gene_coverages_path, gene_detection_path, run):
    """
    Submits fixed-size gene chunks to a process pool capped at num_workers.
    At most num_workers * chunk_size genes are resident in memory at any time.
    """
    num_chunks = len(gene_chunks)
    args_dict = vars(args).copy()
    args_dict.pop('split_names_of_interest', None)

    tmp_dir = tempfile.mkdtemp(prefix='anvi-export-gene-coverage-', dir=os.getcwd())
    samples = None

    try:
        chunk_cov_paths = {}
        chunk_det_paths = {}

        with ProcessPoolExecutor(max_workers=num_workers) as executor:
            futures = {
                executor.submit(process_gene_chunk, args_dict, chunk, idx, tmp_dir): idx
                for idx, chunk in enumerate(gene_chunks)
            }

            for future in as_completed(futures):
                idx = futures[future]
                cov_path, det_path, chunk_samples = future.result()
                chunk_cov_paths[idx] = cov_path
                chunk_det_paths[idx] = det_path
                if samples is None:
                    samples = chunk_samples
                run.info_single(f"Chunk {idx + 1}/{num_chunks} complete.", mc='cyan')

        # Write header once, then stream each chunk file into the final output
        header_line = '\t'.join(['gene_callers_id'] + samples) + '\n'

        for final_path, chunk_paths in [
            (gene_coverages_path, chunk_cov_paths),
            (gene_detection_path, chunk_det_paths),
        ]:
            with open(final_path, 'w') as out_fh:
                out_fh.write(header_line)
                for idx in sorted(chunk_paths):
                    with open(chunk_paths[idx]) as in_fh:
                        shutil.copyfileobj(in_fh, out_fh)

    finally:
        shutil.rmtree(tmp_dir, ignore_errors=True)


def get_args():
    parser = ArgumentParser(description=__description__)

    groupA = parser.add_argument_group('DATABASES', "Anvi'o databases to read from")
    groupA.add_argument(*anvio.A('profile-db'), **anvio.K('profile-db'))
    groupA.add_argument(*anvio.A('contigs-db'), **anvio.K('contigs-db'))

    groupB = parser.add_argument_group('OUTPUT', "Define a prefix for your output files")
    groupB.add_argument(*anvio.A('output-file-prefix'), **anvio.K('output-file-prefix', {'required': True}))

    groupC = parser.add_argument_group('GENES', "Gene calls you want to work with. Without these "
                                       "parameters anvi'o will report everything it finds in the "
                                       "profile database.")
    groupC.add_argument(*anvio.A('gene-caller-id'), **anvio.K('gene-caller-id'))
    groupC.add_argument(*anvio.A('genes-of-interest'), **anvio.K('genes-of-interest'))

    groupD = parser.add_argument_group('PERFORMANCE', "Parallelism and memory options")
    groupD.add_argument('--num-threads', metavar='INT', type=int, default=1,
                        help="Maximum number of worker processes to run concurrently. "
                             "When greater than 1, all genes (filtered or not) are processed "
                             "in parallel chunks (default: %(default)s).")
    groupD.add_argument('--chunk-size', metavar='INT', type=int, default=500,
                        help="Number of genes per chunk. Smaller values reduce peak memory per "
                             "worker at the cost of more ProfileSuperclass init overhead. "
                             "The total number of chunks will be ceil(genes / chunk-size), "
                             "processed in a rolling pool of --num-threads workers "
                             "(default: %(default)s).")

    return parser.get_args(parser)


if __name__ == '__main__':
    main()


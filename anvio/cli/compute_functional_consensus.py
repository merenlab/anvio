#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Compute functional consensus across gene clusters."""

import sys

import anvio
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths

from anvio.errors import ConfigError, FilesNPathsError

__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['cvanni']
__tags__ = ['pangenomics', 'functions', 'clustering']
__requires__ = ['gene-clusters-txt', 'contigs-db', 'genomes-storage-db']
__provides__ = ['functions-txt']
__description__ = (
    "Evaluate functional annotation consensus within gene clusters from any "
    "clustering algorithm, compute within-source Jaccard-based consistency and "
    "cross-source coherence scores, classify clusters (PURE / HIGH_CONSENSUS / "
    "MIXED / LOW_EVIDENCE), and optionally propagate dominant annotations to "
    "unannotated genes."
)


@terminal.time_program
def main():
    args = get_args()
    try:
        run_program(args)
    except ConfigError as e:
        print(e)
        sys.exit(-1)
    except FilesNPathsError as e:
        print(e)
        sys.exit(-1)


def run_program(args):
    from anvio.clusterannotation import GeneClusterFunctionalConsensus

    run = terminal.Run()

    # ------------------------------------------------------------------ #
    # Resolve gene-cluster source                                          #
    # ------------------------------------------------------------------ #
    if args.pan_db and args.gene_clusters_txt:
        raise ConfigError(
            "Please use either --pan-db or --gene-clusters-txt, not both."
        )
    if not args.pan_db and not args.gene_clusters_txt:
        raise ConfigError(
            "You must supply gene clusters via --pan-db or --gene-clusters-txt."
        )

    # Load gene clusters
    if args.pan_db:
        gene_clusters = _gene_clusters_from_pan_db(args.pan_db)
        run.info('Gene clusters source', f"Pan database: {args.pan_db}")
    else:
        filesnpaths.is_file_exists(args.gene_clusters_txt)
        gene_clusters = args.gene_clusters_txt
        run.info('Gene clusters source', f"File: {args.gene_clusters_txt}")

    # ------------------------------------------------------------------ #
    # Resolve annotation source                                            #
    # ------------------------------------------------------------------ #
    n_annotation_sources = sum([
        bool(args.contigs_db),
        bool(args.genomes_storage),
        bool(args.external_genomes),
    ])
    if n_annotation_sources == 0:
        raise ConfigError(
            "You must provide an annotation source: --contigs-db (single genome), "
            "--genomes-storage, or --external-genomes."
        )
    if n_annotation_sources > 1:
        raise ConfigError(
            "Please provide only one of: --contigs-db, --genomes-storage, "
            "--external-genomes."
        )

    if args.contigs_db:
        annotations = args.contigs_db
        run.info('Annotation source', f"Contigs database: {args.contigs_db}")
    elif args.genomes_storage:
        annotations = args.genomes_storage
        run.info('Annotation source', f"Genomes storage: {args.genomes_storage}")
    else:
        annotations = _parse_external_genomes(args.external_genomes)
        run.info('Annotation source', f"External genomes: {args.external_genomes} "
                                      f"({len(annotations)} genomes)")

    # ------------------------------------------------------------------ #
    # Validate output                                                      #
    # ------------------------------------------------------------------ #
    filesnpaths.is_output_file_writable(args.output_file_prefix + '-consensus.txt')

    # ------------------------------------------------------------------ #
    # Annotation sources filter                                            #
    # ------------------------------------------------------------------ #
    annotation_sources = None
    if args.annotation_sources:
        annotation_sources = [s.strip() for s in args.annotation_sources.split(',')]

    # ------------------------------------------------------------------ #
    # Run                                                                  #
    # ------------------------------------------------------------------ #
    c = GeneClusterFunctionalConsensus(
        gene_clusters=gene_clusters,
        annotations=annotations,
        annotation_sources=annotation_sources,
        min_coverage=args.min_coverage,
        min_jacc_sc=args.min_jacc_sc,
        min_cross_source_coherence=args.min_cross_source_coherence,
        propagate=not args.skip_propagation,
        pfam_data_dir=getattr(args, 'pfam_data_dir', None),
        cogs_data_dir=getattr(args, 'cog_data_dir', None),
    )

    consensus_df, propagated_df = c.compute()

    # ------------------------------------------------------------------ #
    # Write outputs                                                        #
    # ------------------------------------------------------------------ #
    consensus_path = args.output_file_prefix + '-consensus.txt'
    consensus_df.to_csv(consensus_path, sep='\t', index=False)
    run.info('Consensus table', consensus_path, mc='green')

    if not propagated_df.empty:
        propagated_path = args.output_file_prefix + '-propagated.txt'
        propagated_df.to_csv(propagated_path, sep='\t', index=False)
        run.info('Propagated annotations', propagated_path, mc='green')


# ---------------------------------------------------------------------- #
# Helpers                                                                  #
# ---------------------------------------------------------------------- #

def _gene_clusters_from_pan_db(pan_db_path: str):
    """Load gene clusters from an anvi'o pan database into a DataFrame."""
    import pandas as pd
    import anvio.db as anvio_db
    import anvio.tables as t
    from anvio.filesnpaths import is_file_exists

    is_file_exists(pan_db_path)
    _db = anvio_db.DB(pan_db_path, None, ignore_version=True)
    df = _db.get_table_as_dataframe(
        t.pan_gene_clusters_table_name,
        columns_of_interest=['gene_caller_id', 'gene_cluster_id', 'genome_name'],
    )
    _db.disconnect()

    df = df.rename(columns={
        'gene_caller_id': 'gene_callers_id',
        'gene_cluster_id': 'cluster_id',
    })
    return df


def _parse_external_genomes(external_genomes_path: str) -> list:
    """Parse an anvi'o external-genomes file into a list of (name, path) tuples."""
    import anvio.utils as utils
    from anvio.filesnpaths import is_file_exists

    is_file_exists(external_genomes_path)
    d = utils.get_TAB_delimited_file_as_dictionary(external_genomes_path)
    pairs = []
    for name, row in d.items():
        db_path = row.get('contigs_db_path')
        if not db_path:
            from anvio.errors import ConfigError
            raise ConfigError(
                f"External genomes file '{external_genomes_path}' is missing "
                f"a 'contigs_db_path' column."
            )
        is_file_exists(db_path)
        pairs.append((name, db_path))
    return pairs


# ---------------------------------------------------------------------- #
# Argument parser                                                          #
# ---------------------------------------------------------------------- #

def get_args():
    from anvio.argparse import ArgumentParser

    parser = ArgumentParser(description=__description__)

    # ---- Gene clusters -------------------------------------------------
    gc = parser.add_argument_group(
        'GENE CLUSTER INPUT',
        "Where do the gene clusters come from?  Provide exactly one."
    )
    gc.add_argument(
        *anvio.A('pan-db'), **anvio.K('pan-db', {'required': False}),
    )
    gc.add_argument(
        *anvio.A('gene-clusters-txt'),
        **anvio.K('gene-clusters-txt', {'required': False}),
    )

    # ---- Annotations ---------------------------------------------------
    ann = parser.add_argument_group(
        'ANNOTATION INPUT',
        "Where do the functional annotations come from?  Provide exactly one."
    )
    ann.add_argument(
        *anvio.A('contigs-db'), **anvio.K('contigs-db', {'required': False}),
    )
    ann.add_argument(
        *anvio.A('genomes-storage'), **anvio.K('genomes-storage', {'required': False}),
    )
    ann.add_argument(
        *anvio.A('external-genomes'), **anvio.K('external-genomes', {'required': False}),
    )

    # ---- Source filter -------------------------------------------------
    src = parser.add_argument_group(
        'ANNOTATION SOURCES',
        "Which functional annotation sources to evaluate."
    )
    src.add_argument(
        *anvio.A('annotation-sources'),
        **anvio.K('annotation-sources'),
    )

    # ---- Thresholds ----------------------------------------------------
    thr = parser.add_argument_group(
        'THRESHOLDS',
        "Quality thresholds.  Defaults are sensible starting points."
    )
    thr.add_argument(
        '--min-coverage',
        type=float,
        default=0.5,
        metavar='FLOAT',
        help="Minimum fraction of genes in a cluster that must be annotated "
             "for a source to be evaluated.  Default: 0.5",
    )
    thr.add_argument(
        '--min-jacc-sc',
        type=float,
        default=0.6,
        metavar='FLOAT',
        help="Minimum coverage-scaled Jaccard score to classify a cluster as "
             "HIGH_CONSENSUS (rather than MIXED).  Default: 0.6",
    )
    thr.add_argument(
        '--min-cross-source-coherence',
        type=float,
        default=0.5,
        metavar='FLOAT',
        help="Minimum cross-source keyword similarity before a PURE or "
             "HIGH_CONSENSUS cluster is downgraded to MIXED.  Default: 0.5",
    )

    # ---- Output --------------------------------------------------------
    out = parser.add_argument_group('OUTPUT')
    out.add_argument(
        '-o', '--output-file-prefix',
        required=True,
        metavar='PREFIX',
        help="Output file prefix.  The program writes two files: "
             "PREFIX-consensus.txt and PREFIX-propagated.txt",
    )
    out.add_argument(
        '--skip-propagation',
        default=False,
        action='store_true',
        help="Do not produce the propagated annotations output file.",
    )

    # ---- Optional data dirs --------------------------------------------
    opt = parser.add_argument_group(
        'OPTIONAL DATA DIRECTORIES',
        "Override the default paths to Pfam and COG data used for the "
        "functionally_coherent sub-flag."
    )
    opt.add_argument(
        '--pfam-data-dir',
        default=None,
        metavar='PATH',
        help="Path to the Pfam data directory (must contain Pfam-A.clans.tsv).",
    )
    opt.add_argument(
        '--cog-data-dir',
        default=None,
        metavar='PATH',
        help="Path to the COG data directory.",
    )

    return parser.get_args(parser)


if __name__ == '__main__':
    main()

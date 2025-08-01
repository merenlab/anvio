#!/usr/bin/env python
# -*- coding: utf-8

import sys
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
__requires__ = ["profile-db", "contigs-db"]
__description__ = ("Export gene coverage and detection data for all genes associated with "
                   "contigs described in a profile database")


def main():
    args = get_args()
    run = terminal.Run()
    p = terminal.pluralize

    try:
        if args.gene_caller_id and args.genes_of_interest:
            raise ConfigError("You either should define a single gene caller id, or a list of "
                              "gene callers of interest. Not both")

        if args.genes_of_interest:
            genes_of_interest = set([g.strip() for g in utils.get_column_data_from_TAB_delim_file(args.genes_of_interest, column_indices=[0], expected_number_of_fields=1)[0] if g])
        elif args.gene_caller_id:
            genes_of_interest = set([str(args.gene_caller_id).strip()])
        else:
            genes_of_interest = None

        if not genes_of_interest:
            splits_of_interest = None
        else:
            # so there are some gene calls of interest. which means, we will first want to identify
            # which split names they are associated with so we can pass that information to the profile
            # super through `args.split_names_of_interest` variable to limit what is initialized.
            contigs_db = dbops.ContigsDatabase(args.contigs_db)
            where_clause = f"gene_callers_id in ({','.join(genes_of_interest)})"
            splits_of_interest = set(contigs_db.db.get_single_column_from_table("genes_in_splits", "split", unique=True, where_clause=where_clause))

            if not len(splits_of_interest):
                raise ConfigError(f"You have provided {p('gene caller id', len(genes_of_interest))} of interest but there were no split names in the contigs "
                                  f"matching them :/")
            else:
                run.warning(f"You have provided {p('gene caller id', len(genes_of_interest))} which matched to {p('split', len(splits_of_interest))}. "
                            f"Very good job.", lc="green")

            # now we have a list of split names of interest, but since not every split occurs in the
            # profile database, we need to make sure those that don't occur in the profile databae are
            # removed before the profile db is initialized.
            split_names_in_profile_db = utils.get_all_item_names_from_the_database(args.profile_db)
            splits_only_in_contigs_db = splits_of_interest.difference(split_names_in_profile_db)

            if splits_only_in_contigs_db:
                run.warning(f"Hear this: {p('split name', len(splits_only_in_contigs_db))} that matched to the gene calls you were interested in "
                            f"occurred only in the contigs db but not in the profile database. Which means, during the profiling, those contigs "
                            f"were not considered (most likely becasue they were shorter than the minimum contig length cutoff). Anvi'o will "
                            f"remove them from the list of splits of interest and you will see what happens downstream.")

                # remove split names that are only in the contigs database from the
                # splits of interest:
                [splits_of_interest.remove(s) for s in splits_only_in_contigs_db]

            if not len(splits_of_interest):
                raise ConfigError("No split names were left to work with :( Genes you are interested in are not occurring on any contig that "
                                  "survived the `anvi-profile` step :/ Sorry!")

            args.split_names_of_interest = splits_of_interest

        profile_db = dbops.ProfileSuperclass(args)

        gene_coverages_path = args.output_file_prefix +  '-GENE-COVERAGES.txt'
        gene_detection_path = args.output_file_prefix +  '-GENE-DETECTION.txt'

        filesnpaths.is_output_file_writable(gene_coverages_path)
        filesnpaths.is_output_file_writable(gene_detection_path)

        gene_coverages_file = open(gene_coverages_path, 'w')
        gene_detection_file = open(gene_detection_path, 'w')

        header_coverages = ['gene_callers_id'] + profile_db.p_meta['samples']
        header_detection = ['gene_callers_id'] + profile_db.p_meta['samples']

        gene_coverages_file.write('\t'.join(header_coverages) + '\n')
        gene_detection_file.write('\t'.join(header_detection) + '\n')

        def write_to_file():
            for gene_callers_id in profile_db.gene_level_coverage_stats_dict:
                line_coverages = [gene_callers_id]
                line_detection = [gene_callers_id]

                for sample_name in profile_db.p_meta['samples']:
                    line_coverages.append(profile_db.gene_level_coverage_stats_dict[gene_callers_id][sample_name]['mean_coverage'])
                    line_detection.append(profile_db.gene_level_coverage_stats_dict[gene_callers_id][sample_name]['detection'])

                gene_coverages_file.write('\t'.join(map(str, line_coverages)) + '\n')
                gene_detection_file.write('\t'.join(map(str, line_detection)) + '\n')

            profile_db.gene_level_coverage_stats_dict = {}
            profile_db.split_coverage_values_per_nt_dict = {}

        if genes_of_interest:
            profile_db.init_gene_level_coverage_stats_dicts(callback=write_to_file, callback_interval=1000, gene_caller_ids_of_interest={int(g) for g in genes_of_interest})
        else:
            profile_db.init_gene_level_coverage_stats_dicts(callback=write_to_file, callback_interval=1000)

        gene_coverages_file.close()
        gene_detection_file.close()

        run.info('Gene coverages', gene_coverages_path)
        run.info('Gene detection', gene_detection_path)
    except ConfigError as e:
        print(e)
        sys.exit(-1)
    except FilesNPathsError as e:
        print(e)
        sys.exit(-1)


def get_args():
    parser = ArgumentParser(description=__description__)

    groupA = parser.add_argument_group('DATABASES', "Anvi'o databases to read from")
    groupA.add_argument(*anvio.A('profile-db'), **anvio.K('profile-db'))
    groupA.add_argument(*anvio.A('contigs-db'), **anvio.K('contigs-db'))

    groupB = parser.add_argument_group('OUTPUT', "Define a prefix for your output files")
    groupB.add_argument(*anvio.A('output-file-prefix'), **anvio.K('output-file-prefix', {'required': True}))

    groupC = parser.add_argument_group('GENES', "Gene calls you want to work with. Without these parameters "
                                       "anvi'o will report everything it finds in the profile database (please) "
                                       "note that the reported genes will only include those that occur in contigs "
                                       "that were taken into consideration during `anvi-profile`, which means "
                                       "if there was a length cutoff for profiling, genes that occur in contigs "
                                       "shorter than that cutoff will not appear in your output.")
    groupC.add_argument(*anvio.A('gene-caller-id'), **anvio.K('gene-caller-id'))
    groupC.add_argument(*anvio.A('genes-of-interest'), **anvio.K('genes-of-interest'))

    return parser.get_args(parser)

if __name__ == '__main__':
    main()

#!/usr/bin/env python
# -*- coding: utf-8

import os
import sys

import anvio
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths
import anvio.taxonomyops.trna as trnataxonomyops

from anvio.errors import ConfigError, FilesNPathsError


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['meren']
__requires__ = ['profile-db', 'contigs-db', 'trna-taxonomy', 'collection', 'bin', 'metagenomes', 'dna-sequence']
__provides__ = ['genome-taxonomy', 'genome-taxonomy-txt']
__resources__ = []
__description__ = ("Estimates taxonomy at genome and metagenome level using tRNA sequences.")


@terminal.time_program
def main():
    args = get_args()

    try:
        if args.metagenomes:
            raise ConfigError("This program can't do this yet.")

        if args.dna_sequence:
            if args.anticodon_for_metagenome_mode:
                anticodon = args.anticodon_for_metagenome_mode
                genomic_sequence = args.dna_sequence
            else:
                # learn about the anticodon
                import anvio.trnaidentifier as trnaidentifier
                profiler = trnaidentifier.Profiler()

                sequence_profile = profiler.profile_gene(args.dna_sequence)
                if not sequence_profile.predicted_profile:
                    raise ConfigError("As far as anvi'o knows, this sequence does not look like a tRNA sequence :/")

                # Allow for args.dna_sequence to be a tRNA transcript: remove the 3'-CCA acceptor
                # regardless of whether it is encoded in the gene or added post-transcriptionally.
                genomic_sequence = args.dna_sequence[: -3] if sequence_profile.has_encoded_acceptor else args.dna_sequence
            anticodon = sequence_profile.predicted_profile.anticodon_seq

            if not anticodon:
                raise ConfigError("Anvi'o couldn't identify the anticodon in this sequence. Perhaps it is not a proper tRNA sequence?")

            # setup search
            fasta_formatted_sequence = f'>{anticodon}\n{genomic_sequence}'

            # go
            log_file_path = filesnpaths.get_temp_file_path(just_the_path=True)
            trnataxonomyops.PopulateContigsDatabaseWithTRNATaxonomy(args).get_gene_estimation_output(anticodon,
                                                                                                     fasta_formatted_sequence,
                                                                                                     log_file_path=log_file_path,
                                                                                                     show_all_hits=True)

            if not anvio.DEBUG:
                os.remove(log_file_path)
        else:
            t = trnataxonomyops.TRNATaxonomyEstimatorSingle(args)
            t.estimate()
    except ConfigError as e:
        print(e)
        sys.exit(-1)
    except FilesNPathsError as e:
        print(e)
        sys.exit(-2)


def get_args():
    from anvio.argparse import ArgumentParser

    parser = ArgumentParser(description=__description__)

    groupA = parser.add_argument_group('INPUT #1', "The minimum you must provide this program is a contigs database. In which case\
                                                    anvi'o will attempt to estimate taxonomy for all contigs in it, assuming that\
                                                    the contigs database represents a single genome. If the contigs database is actually\
                                                    a metagenome, you should use the `--metagenome` flag to explicitly declare that.")
    groupA.add_argument(*anvio.A('contigs-db'), **anvio.K('contigs-db', {'required': False}))
    groupA.add_argument(*anvio.A('metagenome-mode'), **anvio.K('metagenome-mode'))

    groupB = parser.add_argument_group('INPUT #2', "In addition, you can also point out a profile database. In which case you also must\
                                                    provide a collection name. When you do that anvi'o will offer taxonomy estimates for\
                                                    each bin in your collection.")
    groupB.add_argument(*anvio.A('profile-db'), **anvio.K('profile-db', {'required': False}))
    groupB.add_argument(*anvio.A('collection-name'), **anvio.K('collection-name'))

    groupC = parser.add_argument_group('INPUT #3', "You can also work with a metagenomes file, assuming that you have multiple metagenomes\
                                                    with or without associated mapping results, and anvi'o would generate a singe output\
                                                    file for all.")
    groupC.add_argument(*anvio.A('metagenomes'), **anvio.K('metagenomes'))

    groupD = parser.add_argument_group('INPUT #4', "Ad hoc sequence search. No contigs databases, no profiles. The lazy stayla. Please note\
                                                    that if you use parameters defined under this optin, none of the other standard\
                                                    parameters for this program will be taken into consideration.")
    groupD.add_argument(*anvio.A('dna-sequence'), **anvio.K('dna-sequence'))
    groupD.add_argument(*anvio.A('max-num-target-sequences'), **anvio.K('max-num-target-sequences'))

    groupE = parser.add_argument_group('OUTPUT AND FORMATTING', "Anvi'o will do its best to offer you some fancy output tables for your viewing \
                                                                 pleasure by default. But in addition to that, you can ask the resulting information \
                                                                 to be stored in a TAB-delimited file (which is a much better way to include the \
                                                                 results in your study as supplementary information, or work with these results \
                                                                 using other analysis tools such as R). Depending on the mode you are running this \
                                                                 program, anvi'o may ask you to use an 'output file prefix' rather than an 'output \
                                                                 file path'.")
    groupE.add_argument(*anvio.A('output-file'), **anvio.K('output-file'))
    groupE.add_argument(*anvio.A('per-anticodon-output-file'), **anvio.K('per-anticodon-output-file'))
    groupE.add_argument(*anvio.A('output-file-prefix'), **anvio.K('output-file-prefix'))
    groupE.add_argument(*anvio.A('taxonomic-level'), **anvio.K('taxonomic-level', {'default': None}))
    groupE.add_argument(*anvio.A('matrix-format'), **anvio.K('matrix-format', {'help': "If you want the reports to look like sparse matrices whenever "
                                                                                   "possible, declare this flag. Matrices are especially good to use "
                                                                                   "when you are working with internal/external genomes since they can "
                                                                                   "show you quickly the distribution of each taxon across all metagenomes "
                                                                                   "in programs like EXCEL. WELL TRY IT AND SEE."}))
    groupE.add_argument(*anvio.A('raw-output'), **anvio.K('raw-output'))

    groupF = parser.add_argument_group('PERFORMANCE', "We are not sure if allocating more threads for this operation will change anything.\
                                                       But hey. One can try.")
    groupF.add_argument(*anvio.A('num-threads'), **anvio.K('num-threads'))

    groupF = parser.add_argument_group('AUTHORITY', "Assert your dominance.")
    groupF.add_argument(*anvio.A('anticodon-for-metagenome-mode'), **anvio.K('anticodon-for-metagenome-mode'))
    groupF.add_argument(*anvio.A('report-anticodon-frequencies'), **anvio.K('report-anticodon-frequencies'))
    groupF.add_argument(*anvio.A('just-do-it'), **anvio.K('just-do-it'))

    groupG = parser.add_argument_group('ADVANCED', "Very pro-like stuff.")
    groupG.add_argument(*anvio.A('simplify-taxonomy-information'), **anvio.K('simplify-taxonomy-information'))
    groupG.add_argument(*anvio.A('compute-anticodon-coverages'), **anvio.K('compute-anticodon-coverages'))
    groupG.add_argument(*anvio.A('update-profile-db-with-taxonomy'), **anvio.K('update-profile-db-with-taxonomy'))

    groupH = parser.add_argument_group('BORING', "Options that you will likely never need.")
    groupH.add_argument(*anvio.A('taxonomy-database'), **anvio.K('taxonomy-database'))

    return parser.get_args(parser)


if __name__ == '__main__':
    main()

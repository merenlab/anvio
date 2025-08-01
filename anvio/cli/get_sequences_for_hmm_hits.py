#!/usr/bin/env python
# -*- coding: utf-8
"""Get sequences for all HMM hits in a given bin.

   This program takes a profile database, a collection ID, and a bin name, and an
   HMM source, and returnes sequences of HMM hits. This program is useful when you
   want to get actual sequencs for each single-copy gene hit in a particular genome
   bin.

  You want to play with it? This is how you could quickly test it:

  Downloaded the anvi'o data pack for the infant gut data, which is here:

    https://ndownloader.figshare.com/files/8252861

  Unpack it and went into it:

    tar -zxvf INFANTGUTTUTORIAL.tar.gz && cd INFANT-GUT-TUTORIAL

  Import the collection `merens`:

    anvi-import-collection additional-files/collections/merens.txt -p PROFILE.db -c CONTIGS.db -C merens

  Then I run the program `anvi-get-sequences-for-hmm-hits` in the anvi'o master this way:

    anvi-get-sequences-for-hmm-hits -p PROFILE.db \
                                    -c CONTIGS.db \
                                    -C merens \
                                    -o OUTPUT.fa \
                                    --hmm-source Campbell_et_al \
                                    --gene-names Ribosomal_L27,Ribosomal_L28,Ribosomal_L3 \
                                    --return-best-hit \
                                    --get-aa-sequences \
                                    --concatenate
"""

import os
import sys

import anvio
import anvio.utils as utils
import anvio.hmmops as hmmops
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths
import anvio.ccollections as ccollections
import anvio.hmmopswrapper as hmmopswrapper
import anvio.genomedescriptions as genomedescriptions

from anvio.dbops import ContigsSuperclass, ContigsDatabase
from anvio.argparse import ArgumentParser
from anvio.errors import ConfigError, FilesNPathsError


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['meren']
__requires__ = ['contigs-db', 'profile-db', 'external-genomes', 'internal-genomes', 'hmm-source', "hmm-hits"]
__provides__ = ['genes-fasta', 'concatenated-gene-alignment-fasta']
__description__ = "Get sequences for HMM hits from many inputs"
__resources__ = [("A tutorial on anvi'o phylogenomics workflow", "http://merenlab.org/2017/06/07/phylogenomics/"),
                 ("A detailed application of phylogenomics to place a new genome on a tree", "http://merenlab.org/data/parcubacterium-in-hbcfdna/")]


@terminal.time_program
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
    progress = terminal.Progress()
    run = terminal.Run()
    P = terminal.pluralize

    if args.external_genomes and (args.contigs_db or args.profile_db):
        raise ConfigError("If you are interested in using a list of external genomes, you shouldn't be using a contigs db or "
                          "profile db prameters.")

    if not (args.internal_genomes or args.external_genomes or args.contigs_db or args.profile_db):
        raise ConfigError("You gotta give this program some input files :/ Come on.")

    if args.list_defline_variables:
        defline_data_dict = hmmops.SequencesForHMMHits(None).defline_data_dict
        defline_format = hmmops.SequencesForHMMHits(None).defline_format

        run.warning(f"Here are the variables you can use to provide a user-defined defline template: ")
        for key in defline_data_dict.keys():
            run.info_single("{%s}" % key)
        run.info_single("Remember, by default, anvi'o will use the following template to format the deflines of "
                        "FASTA files it produces whenever possible.", level=0, nl_before=1, nl_after=1, mc='red')

        run.info_single("%s" % defline_format, level=0, nl_after=1, cut_after=0)
        sys.exit()

    if args.list_hmm_sources or args.list_available_gene_names:
        if args.contigs_db:
            if args.list_hmm_sources:
                ContigsDatabase(args.contigs_db).list_available_hmm_sources()
            else:
                hmmops.SequencesForHMMHits(args.contigs_db).list_available_gene_names(sources=[s.strip() for s in args.hmm_sources.split(',')] if args.hmm_sources else [])
        elif args.external_genomes or args.internal_genomes:
            genomedescriptions.GenomeDescriptions(args).list_HMM_info_and_quit()

        sys.exit()

    if args.concatenate_genes:
        if not args.return_best_hit:
            raise ConfigError("If you want your genes to be concatenated into a multi-alignment file, you must also ask for "
                              "the best hit (using the `--return-best-hit`) flag to avoid issues if there are more than one "
                              "hit for a gene in a given genome. Anvi'o could have set this flag on your behalf, but it just "
                              "is not that kind of a platform :/")

    if not args.concatenate_genes and args.partition_file:
        raise ConfigError("Partition files are only relevant when you use the flag `--concatenate-genes`.")

    filesnpaths.is_output_file_writable(args.partition_file) if args.partition_file else None

    if args.concatenate_genes:
        # test whether we know about the aligner early on
        from anvio.drivers import Aligners
        Aligners().select(args.align_with)

    if args.max_num_genes_missing_from_bin and not args.gene_names and not args.list_available_gene_names:
        raise ConfigError("You can only use --max-num-genes-missing-from-bin flag if you already know what gene names you are "
                          "interested in (just to make sure you know what you are doing).")

    for param, value in [('--max-num-genes-missing-from-bin', args.max_num_genes_missing_from_bin),
                         ('--min-num-bins-gene-occurs', args.min_num_bins_gene_occurs)]:
        if value is not None:
            try:
                value = int(value)
                assert(value >= 0)
            except:
                raise ConfigError("For obvious reasons, anvi'o expects the parameter %s to be a non-negative integer :/" % param)

    hmm_sources = set([s.strip() for s in args.hmm_sources.split(',')]) if args.hmm_sources else set([])

    # the following if/else block either uses SequencesForHMMHits or the wrapper class
    # SequencesForHMMHitsWrapperForMultipleContigs (if there are multiple contigs files)
    # to get sequences, and construct `splits_dict` depending on the input files. it is
    # shitty code, and can be improved.
    if args.external_genomes or args.internal_genomes:
        s = hmmopswrapper.SequencesForHMMHitsWrapperForMultipleContigs(args, hmm_sources)
        splits_dict = s.splits_dict
    else:
        info_table = hmmops.SequencesForHMMHits(args.contigs_db).hmm_hits_info

        # let's quickly check whether we have all the hmm_sources the user may have
        # requested has anything to do with the ones we have in the database
        if hmm_sources:
            missing_hmm_sources = [s for s in hmm_sources if s not in info_table]
            if(missing_hmm_sources):
                raise ConfigError("At least one of the HMM sources you requested are missing form the HMMs the contigs database "
                                  "knows about :/ Here they are: '%s'" % (', '.join(missing_hmm_sources)))

        if args.list_available_gene_names:
            hmmops.SequencesForHMMHits(args.contigs_db).list_available_gene_names(sources=list(hmm_sources))

        if (args.profile_db and not args.collection_name) or (args.collection_name and not args.profile_db):
            raise ConfigError("You can't use this program with a profile database but without a collection name,\
                               and vice versa, you also can't use a collection if you didn't provide a profile database. Yes. Because.")

        if args.profile_db:
            utils.is_profile_db_and_contigs_db_compatible(args.profile_db, args.contigs_db)
            splits_dict = ccollections.GetSplitNamesInBins(args).get_dict()
            run.info('Init', '%d splits in %d bin(s)' % (sum([len(v) for v in list(splits_dict.values())]), len(splits_dict)))
        else:
            contigs_db = ContigsSuperclass(args, r = run, p = progress)
            contigs_db_name = os.path.basename(args.contigs_db[:-3])
            splits_dict = {contigs_db_name: list(contigs_db.splits_basic_info.keys())}

        defline_format = args.defline_format if args.defline_format else None
        s = hmmops.SequencesForHMMHits(args.contigs_db, sources=hmm_sources, defline_format=defline_format)

    CHK = lambda: exec('raise ConfigError("Your selections returned an empty list of genes to work with :/")') if not len(hmm_sequences_dict) else None

    hmm_sequences_dict = s.get_sequences_dict_for_hmm_hits_in_splits(splits_dict, return_amino_acid_sequences=args.get_aa_sequences)
    CHK()

    run.info('Sources', f"{', '.join(hmm_sources)}")
    run.info('Hits', '%d HMM hits for %d source(s)' % (len(hmm_sequences_dict), len(s.sources)))

    # if user requested AA sequences, let's check if all or some of them are empty
    if args.get_aa_sequences:
        hits_with_empty_aa_seqs = [h for h in hmm_sequences_dict if not hmm_sequences_dict[h]['sequence']]
        if hits_with_empty_aa_seqs:
            if len(hits_with_empty_aa_seqs) == len(hmm_sequences_dict):
                raise ConfigError("You requested amino acid sequences with the `--get-aa-sequences`, but none of the "
                                  "genes for your requested HMM source(s) have AA sequences associated with them. This often "
                                  "happens with ribosomal RNA genes, for example. Basically, the only way to get sequences for "
                                  "these HMM hits is to get rid of the `--get-aa-sequences` flag.")
            else:
                gene_names = [hmm_sequences_dict[h]['gene_name'] for h in hits_with_empty_aa_seqs]
                gene_names_str = ", ".join(gene_names)
                run.warning(f"Some of the HMM hits you requested do not have amino acid sequences associated with them. "
                            f"Their entries in the output FASTA file will be empty. Here are the gene names of each hit "
                            f"that is missing an AA sequence: {gene_names_str}")

    # keep track of bins removed from the analysis results due to various filters:
    bins_removed_for_any_reason = set([])

    # figure out gene names.. if the user provided a file, use that, otherwhise parse gene names out of the comma-separated text
    if args.gene_names and filesnpaths.is_file_exists(args.gene_names, dont_raise=True):
        gene_names = [g.strip() for g in open(args.gene_names, 'r').readlines()] if args.gene_names else []
    else:
        gene_names = [g.strip() for g in args.gene_names.split(',')] if args.gene_names else []

    run.info('Genes of interest', f"{', '.join(gene_names) if gene_names else None}")

    if len(gene_names):
        hmm_sequences_dict, bins_removed = s.filter_hmm_sequences_dict_for_to_only_include_specific_genes(hmm_sequences_dict, gene_names)
        CHK()

        if len(bins_removed):
            bins_removed_for_any_reason.update(bins_removed)

        run.info('Filtered hits', '%d hits remain after filtering for %d gene(s)' % (len(hmm_sequences_dict), len(gene_names)))

    if args.ignore_genes_longer_than > 0:
        hmm_sequences_dict, gene_calls_removed, bins_removed = s.filter_hmm_sequences_dict_for_genes_that_are_too_long(hmm_sequences_dict, int(args.ignore_genes_longer_than))
        CHK()

        if len(bins_removed):
            bins_removed_for_any_reason.update(bins_removed)

        if len(gene_calls_removed):
            run.info('Filtered hits', '%d hits remain after filtering for %d genes longer than %d' % (len(hmm_sequences_dict), len(gene_calls_removed), args.ignore_genes_longer_than), nl_before=1)

    if args.max_num_genes_missing_from_bin is not None:
        hmm_sequences_dict, bins_removed = s.filter_hmm_sequences_dict_for_bins_that_lack_more_than_N_genes(hmm_sequences_dict, gene_names, int(args.max_num_genes_missing_from_bin))
        CHK()

        if len(bins_removed):
            bins_removed_for_any_reason.update(bins_removed)
            run.info('Filtered hits', '%d hits remain after filtering for `--max-num-genes-missing-from-bin` flag' % (len(hmm_sequences_dict)))

    if args.min_num_bins_gene_occurs is not None:
        hmm_sequences_dict, genes_removed = s.filter_hmm_sequences_dict_from_genes_that_occur_in_less_than_N_bins(hmm_sequences_dict, int(args.min_num_bins_gene_occurs))
        CHK()

        if len(genes_removed):
            run.info('Filtered hits', '%d hits remain after filtering for `--min-num-bins-gene-occurs` flag' % (len(hmm_sequences_dict)))

            run.warning("The `--min-num-bins-gene-occurs` parameter caused the removal of %d genes from your analysis because "
                        "they occurred in less than %d bins/genomes in your analysis. This is the list of genes that gon' "
                        "buhbye: %s." % (len(genes_removed), int(args.min_num_bins_gene_occurs), ', '.join(genes_removed)))

            # update the gene names variable .. this is such a mess :( "WHO WROTE THIS SHIT CODE", yelled Meren, looking
            # at this cursed main function in his office on a Saturday night, knowing very well who did it. thanks
            # to his lack of shame he said to himself "well, I guess it is OK if it stays like this for now".
            if(gene_names):
                gene_names = [g for g in gene_names if g not in genes_removed]

    if args.return_best_hit:
        run.warning("You requested only the best hits to be reported, which means, if, say, there are more than one RecA "
                    "hits in a bin for a given HMM source, only the one with the lowest e-value will be kept, and others "
                    "will be removed from your final results.")

        if not args.profile_db:
            run.warning("You requested to get only the best hits, but you did not provide a profile database. At this point "
                        "anvi'o just hopes you know what you are doing. Since this is like the zone of 'potentially a terrible "
                        "idea but it may be quite relevant when done right'.")

        hmm_sequences_dict = s.filter_hmm_sequences_dict_for_splits_to_keep_only_best_hits(hmm_sequences_dict)
        CHK()

        run.info('Filtered hits', '%d hits remain after removing weak hits for multiple genes' % (len(hmm_sequences_dict)))

    if args.unique_genes:
        run.warning("You asked anvi'o for each gene that has multiple hits to report only the best hit, which means, if, "
                    "say, there are multiple models in your HMM source that match to a single gene, only the one with the "
                    "lowest e-value will be kept, and others will be removed from your final results.")

        hmm_sequences_dict = s.filter_hmm_sequences_dict_to_keep_only_unique_gene_hits(hmm_sequences_dict)
        CHK()

        run.info('Filtered hits', '%d hits remain after removing weak hits for multiple models' % (len(hmm_sequences_dict)))

    if len(bins_removed_for_any_reason):
            run.warning(f"Please note that a total of {P('bin', len(bins_removed_for_any_reason))} from your final results "
                        f"were removed while anvi'o was applying all the filters to your HMM hits. Please carefully "
                        f"review the logs above to make sure you have a good grasp on what happened, and you're happy "
                        f"with the reults.", header='CONDOLENCES FROM DAAN TO YOU FOR YOUR BINS THAT WENT 💀')

    if args.separator:
        separator = args.separator
    else:
        separator = 'XXX' if args.get_aa_sequences else 'NNN'

    # make sure an output file name is provided.
    if not args.output_file:
        args.output_file = 'sequences-for-hmm-hits.fa'

    # the magic is happening here:
    s.store_hmm_sequences_into_FASTA(hmm_sequences_dict,
                                     args.output_file,
                                     concatenate_genes=args.concatenate_genes,
                                     partition_file_path=args.partition_file,
                                     separator=separator,
                                     wrap=None if args.no_wrap else 120,
                                     genes_order=list(gene_names) if len(gene_names) else None,
                                     align_with=args.align_with,
                                     just_do_it=args.just_do_it)

    run.info('Mode', 'AA sequences' if args.get_aa_sequences else 'DNA sequences', mc='green')
    run.info('Genes are concatenated', args.concatenate_genes)
    run.info('Output', args.output_file)


def get_args():
    parser = ArgumentParser(description=__description__)

    groupA = parser.add_argument_group('INPUT OPTION #1: CONTIGS DB', "There are multiple ways to access to sequences. Your first option is to\
                                        provide a contigs database, and call it a day. In this case the program will return you\
                                        everything from it.")
    groupA.add_argument(*anvio.A('contigs-db'), **anvio.K('contigs-db', {'required': False}))

    groupB = parser.add_argument_group('INPUT OPTION #2: CONTIGS DB + PROFLIE DB', "You can also work with anvi'o profile databases and collections\
                                        stored in them. If you go this way, you still will need to provide a contigs database. If you\
                                        just specify a collection name, you will get hits from every bin in it. You can also use\
                                        the bin name or bin ids file parameters to specify your interest more precisely.")
    groupB.add_argument(*anvio.A('profile-db'), **anvio.K('profile-db', {'required': False}))
    groupB.add_argument(*anvio.A('collection-name'), **anvio.K('collection-name'))
    groupB.add_argument(*anvio.A('bin-id'), **anvio.K('bin-id'))
    groupB.add_argument(*anvio.A('bin-ids-file'), **anvio.K('bin-ids-file'))

    groupC = parser.add_argument_group('INPUT OPTION #3: INT/EXTERNAL GENOMES FILE', "Yes. You can alternatively use as input an internal or external\
                                        genomes file, or both of them together. If you have multiple contigs databases without any profile\
                                        database, you can use the external genomes file. So if you just have a bunch of FASTA files and nothing else,\
                                        this is what you need. In contrast, if you want to access to genes in bins described in collections\
                                        stored in anvi'o profile databases, then you can use internal genomes file route. Or you can mix the two,\
                                        because why not. There is not much room for excuses here.")
    groupC.add_argument(*anvio.A('external-genomes'), **anvio.K('external-genomes'))
    groupC.add_argument(*anvio.A('internal-genomes'), **anvio.K('internal-genomes'))

    groupD = parser.add_argument_group('HMM STUFF', "This is where you can specify an HMM source, and/or a list of genes to filter\
                                        your results.")
    groupD.add_argument(*anvio.A('hmm-sources'), **anvio.K('hmm-sources'))
    groupD.add_argument(*anvio.A('gene-names'), **anvio.K('gene-names'))
    groupD.add_argument(*anvio.A('list-hmm-sources'), **anvio.K('list-hmm-sources'))
    groupD.add_argument(*anvio.A('list-available-gene-names'), **anvio.K('list-available-gene-names'))

    groupE = parser.add_argument_group('THE OUTPUT', "Where should the output go. It will be a FASTA file, and you better give it\
                                        a nice name..")
    groupE.add_argument(*anvio.A('output-file'), **anvio.K('output-file'))
    groupE.add_argument(*anvio.A('no-wrap'), **anvio.K('no-wrap'))
    groupE.add_argument(*anvio.A('list-defline-variables'), **anvio.K('list-defline-variables'))
    groupE.add_argument(*anvio.A('defline-format'), **anvio.K('defline-format', {'default': None}))

    groupF = parser.add_argument_group('THE ALPHABET', "The sequences are reported in DNA alphabet, but you can also get them\
                                        translated just like all the other cool kids.")
    groupF.add_argument(*anvio.A('get-aa-sequences'), **anvio.K('get-aa-sequences'))

    groupG = parser.add_argument_group('PHYLOGENOMICS? K!', "If you want, you can get your sequences concatanated. In this case\
                                       anwi'o will use muscle to align every homolog, and concatenate them the order you specified\
                                       using the `gene-names` argument. Each concatenated sequence will be separated from the other\
                                       ones by the `separator`.")
    groupG.add_argument(*anvio.A('concatenate-genes'), **anvio.K('concatenate-genes'))
    groupG.add_argument(*anvio.A('partition-file'), **anvio.K('partition-file'))
    groupG.add_argument(*anvio.A('max-num-genes-missing-from-bin'), **anvio.K('max-num-genes-missing-from-bin'))
    groupG.add_argument(*anvio.A('min-num-bins-gene-occurs'), **anvio.K('min-num-bins-gene-occurs'))
    groupG.add_argument(*anvio.A('align-with'), **anvio.K('align-with'))
    groupG.add_argument(*anvio.A('separator'), **anvio.K('separator', {'help': 'A word that will be used to\
                                  sepaate concatenated gene sequences from each other (IF you are using this\
                                  program with `--concatenate-genes` flag). The default is "XXX" for amino\
                                  acid sequences, and "NNN" for DNA sequences'}))

    groupH = parser.add_argument_group('OPTIONAL', "Everything is optional, but some options are more optional than others.")
    groupH.add_argument(*anvio.A('return-best-hit'), **anvio.K('return-best-hit'))
    groupH.add_argument(*anvio.A('ignore-genes-longer-than'), **anvio.K('ignore-genes-longer-than'))
    groupH.add_argument(*anvio.A('unique-genes'), **anvio.K('unique-genes'))
    groupH.add_argument(*anvio.A('just-do-it'), **anvio.K('just-do-it'))

    return parser.get_args(parser)


if __name__ == '__main__':
    main()

#!/usr/bin/env python
# -*- coding: utf-8

import os
import sys

import anvio
import anvio.utils as utils
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths
import anvio.ccollections as ccollections
import anvio.genomedescriptions as genomedescriptions

from anvio.completeness import Completeness
from anvio.errors import ConfigError, FilesNPathsError
from anvio.dbops import ContigsSuperclass


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['meren']
__requires__ = ["contigs-db", "profile-db", "external-genomes", "collection"]
__provides__ = ["completion"]
__description__ = "Estimate completion and redundancy using domain-specific single-copy core genes"


pp = terminal.pretty_print


@terminal.time_program
def main():
    args = get_args()

    A = lambda x: args.__dict__[x] if x in args.__dict__ else None
    contigs_db_path = A('contigs_db')
    profile_db_path = A('profile_db')
    collection_name = A('collection_name')
    list_collections = A('list_collections')
    output_file_path = A('output_file')
    external_genomes = A('external_genomes')
    internal_genomes = A('internal_genomes')

    try:
        if contigs_db_path:
            if external_genomes or internal_genomes:
                raise ConfigError("You can't work with a single contigs database and a set of external genomes "
                                  "at the same time.")

            if not profile_db_path and list_collections:
                raise ConfigError("When there is no profile database involved in your analysis, there is no "
                                  "collections to list really :/")

            if profile_db_path:
                utils.is_profile_db_and_contigs_db_compatible(profile_db_path, contigs_db_path)

        if external_genomes or internal_genomes:
            if profile_db_path:
                raise ConfigError("If you want to work with a set of external genomes, there is no use of "
                                  "declaring a profile database :/")
            if collection_name:
                raise ConfigError("When working with external genomes, collection names are automatically "
                                  "determined from the name column of the external genomes file.")

            if list_collections:
                raise ConfigError("You are working with external genomes, there are no collections to list :/")


        if output_file_path:
            filesnpaths.is_output_file_writable(output_file_path)

        if contigs_db_path:
            print_for_a_single_contigs_db(args, contigs_db_path, profile_db_path, collection_name, list_collections)
        elif internal_genomes or external_genomes:
            print_for_internal_or_external_genomes(args)
        else:
            raise ConfigError("We're lost :(")

    except ConfigError as e:
        print(e)
        sys.exit(-1)
    except FilesNPathsError as e:
        print(e)
        sys.exit(-1)


def print_proper(args, collection, completeness, run, progress):
    A = lambda x: args.__dict__[x] if x in args.__dict__ else None
    output_file_path = A('output_file')

    contigs_db = ContigsSuperclass(args, r=terminal.Run(verbose=False), p=terminal.Progress())

    table = []
    header=['bin name', 'domain', 'confidence', '% completion', '% redundancy', 'num_splits', 'total length']
    for bin_id in collection:
        p_completion, p_redundancy, domain, domain_probabilities, info_text, _ = completeness.get_info_for_splits(set(collection[bin_id]), bin_name=bin_id)
        total_length = sum([contigs_db.splits_basic_info[split_name]['length'] for split_name in set(collection[bin_id])])
        num_splits = len(collection[bin_id])
        domain_confidence = domain_probabilities[domain] if domain else 0.0

        table.append([bin_id, '%s' % domain.upper(), '%.1f' % domain_confidence, '%.2f' % p_completion, '%.2f' % p_redundancy, num_splits, total_length])

    anvio.TABULATE(table, header)

    run.info_single("The 'domain' shown in this table is the domain anvi'o predicted for your contigs in a given bin with the "
                    "amount of confidence for that call in the 'domain' column. If the domain is 'mixed', it means it is very "
                    "likely you have contigs in your bin that spans accross multiple domains of life (maybe it is worth a Nobel, "
                    "but more likely it is just garbage). If the domain is 'blank', it means anvi'o did not find enough signal "
                    "from single-copy core genes to determine what domain it should use to estimate the completion and redundancy "
                    "of this bin accurately. You can get much more information about these estimations by running the same command "
                    "line with the additinal flag `--debug`.", nl_before=1, nl_after=1)

    if output_file_path:
        utils.store_array_as_TAB_delimited_file(table, output_file_path, header=header)
        run.info('Results also stored in an output file', output_file_path, nl_after=1)


def print_improper(args, collection, completeness, run, progress):
    A = lambda x: args.__dict__[x] if x in args.__dict__ else None
    output_file_path = A('output_file')

    run.info_single("Hi. I am the function that is called when the completion class is not initialized properly, and the best I can "
                    "do is to dump all the raw the raw C/R estimates from all domains altogether.")

    contigs_db = ContigsSuperclass(args, r=terminal.Run(verbose=False), p=terminal.Progress())

    table = []
    header = ['bin name', 'domain', 'source', '% completion', '% redundancy', 'num_splits', 'total length']

    # we will keep the talbe and header above for the output file, but use a smaller table for
    # printing out information:
    bin_header = ['domain', 'source', '% completion', '% redundancy']
    for bin_id in collection:
        bin_table = []
        _, _, _, _, _, scg_hmm_hits = completeness.get_info_for_splits(set(collection[bin_id]), bin_name=bin_id)
        total_length = sum([contigs_db.splits_basic_info[split_name]['length'] for split_name in set(collection[bin_id])])
        num_splits = len(collection[bin_id])

        run.width = 60
        run.warning(None, header='%s (num_splits: %s; total_length: %s)' % (bin_id, pp(num_splits), pp(total_length)), lc = 'green')

        for domain in scg_hmm_hits:
            for scg_collection_for_domain in scg_hmm_hits[domain]:
                D = scg_hmm_hits[domain][scg_collection_for_domain]
                table.append([bin_id, '%s' % domain.upper(), '%s' % scg_collection_for_domain, '%.2f' % D['percent_completion'], '%.2f' % D['percent_redundancy'], num_splits, total_length])
                bin_table.append(['%s' % domain.upper(), '%s' % scg_collection_for_domain, '%.2f' % D['percent_completion'], '%.2f' % D['percent_redundancy']])

        anvio.TABULATE(bin_table, bin_header)

    if output_file_path:
        utils.store_array_as_TAB_delimited_file(table, output_file_path, header=header)
        run.info('Results also stored in an output file', output_file_path, nl_after=1, nl_before=2)


def print_for_internal_or_external_genomes(args):
    run = terminal.Run(verbose=False if args.concise else True)

    g = genomedescriptions.GenomeDescriptions(args, run=terminal.Run(verbose=False))
    g.load_genomes_descriptions(skip_sanity_check=True)
    genomes = g.genomes

    missing_hmms = [genome_name for genome_name in genomes if not genomes[genome_name]['hmms_for_scgs_were_run']]
    if len(missing_hmms):
        if args.just_do_it:
            for genome_name in missing_hmms:
                genomes.pop(genome_name)

            if not len(genomes):
                raise ConfigError("Anvi'o was ignoring genomes that were missing HMMs for single-cop core genes "
                                  "and left with 0 genomes in total. Great job.")

            run.warning("%d of your genomes in the external genomes file were missing HMMs for single-copy core "
                        "genes, and since you asked anvi'o to `--just-do-it`, it removed these genomes from the "
                        "list of genomes to be considered. You are left with %d." % (len(missing_hmms), len(genomes)), overwrite_verbose=True)
        else:
            if anvio.DEBUG:
                raise ConfigError("%d of your genomes in the external genomes file are missing HMMs for single-copy core "
                                  "genes. Here is the full list of them: %s." % (len(missing_hmms), ','.join(missing_hmms)))
            else:
                raise ConfigError("%d of %d of your genomes in the external genomes file are missing HMMs for single-copy core "
                                  "genes. To see the full list, run the same command with `--debug` flag. If you want avni'o "
                                  "to simply ignore the ones that do not have HMMs, run the same command wiht `--just-do-it` flag." \
                                                    % (len(missing_hmms), len(genomes)))

    A = lambda x: args.__dict__[x] if x in args.__dict__ else None
    output_file_path = A('output_file')

    table = []
    header=['genome name', 'domain', 'confidence', '% completion', '% redundancy', 'num_splits', 'total length']
    for genome in genomes.values():
        p_completion = genome['percent_completion']
        p_redundancy = genome['percent_redundancy']
        domain = genome['scg_domain']
        domain_confidence = genome['scg_domain_confidence']
        total_length = genome['total_length']
        num_splits = genome['num_splits']

        table.append([genome['name'], '%s' % domain.upper(), '%.1f' % domain_confidence, '%.2f' % p_completion, '%.2f' % p_redundancy, num_splits, total_length])

    anvio.TABULATE(table, header)

    run.info_single("The 'domain' shown in this table is the domain anvi'o predicted for your contigs in a given bin with the "
                    "amount of confidence for that call in the 'domain' column. If the domain is 'mixed', it means it is very "
                    "likely you have contigs in your bin that spans accross multiple domains of life (maybe it is worth a Nobel, "
                    "but more likely it is just garbage). If the domain is 'blank', it means anvi'o did not find enough signal "
                    "from single-copy core genes to determine what domain it should use to estimate the completion and redundancy "
                    "of this bin accurately. You can get much more information about these estimations by running the same command "
                    "line with the additinal flag `--debug`.", nl_before=1, nl_after=1)

    if output_file_path:
        utils.store_array_as_TAB_delimited_file(table, output_file_path, header=header)
        run.info('Results also stored in an output file', output_file_path, nl_after=1)


def print_for_a_single_contigs_db(args, contigs_db_path, profile_db_path=None, collection_name=None, list_collections=False):
    run = terminal.Run(verbose=False if args.concise else True)
    progress = terminal.Progress(verbose=False if args.concise else True)

    collections = ccollections.Collections(r=run, p=progress)
    collections.populate_collections_dict(contigs_db_path)
    if profile_db_path:
        collections.populate_collections_dict(profile_db_path)

    if list_collections:
        collections.list_collections()
        sys.exit()

    if profile_db_path and not collection_name:
        raise ConfigError("You can't use this program with a profile database but without a collection name. Yes. Because.")

    if contigs_db_path and not collection_name:
        run.warning("You did not provide a collection name. Anvi'o will assume that what is in your contigs database "
                    "is a single genome (or genome bin).")

    if collection_name and collection_name not in collections.collections_dict:
        raise ConfigError("%s is not a valid collection ID. See a list of available ones with '--list-collections' flag" % collection_name)

    completeness = Completeness(contigs_db_path)
    if not len(completeness.sources):
        raise ConfigError("HMM's were not run for this contigs database :/")

    if collection_name:
        collection = collections.get_collection_dict(collection_name)
    else:
        contigs_db = ContigsSuperclass(args, r=terminal.Run(verbose=False), p=terminal.Progress())
        collection = {os.path.basename(contigs_db_path[:-3]): list(contigs_db.splits_basic_info.keys())}

    if collection_name:
        run.warning(None, header = 'Bins in collection "%s"' % collection_name, lc = 'green')
    else:
        run.warning(None, header = 'Genome in "%s"' % os.path.basename(contigs_db_path), lc = 'green')

    if not completeness.initialized_properly or args.report_all_estimates:
        print_improper(args, collection, completeness, run, progress)
    else:
        print_proper(args, collection, completeness, run, progress)


def get_args():
    from anvio.argparse import ArgumentParser

    parser = ArgumentParser(description=__description__)

    groupA = parser.add_argument_group('MANDATORY INPUT OPTION #1', "Minimum input is an anvi'o contigs database. If you provide nothing else,\
                                                                     anvi'o will assume that it is a single genome (even if it is not), and\
                                                                     give you back what you need.")
    groupA.add_argument(*anvio.A('contigs-db'), **anvio.K('contigs-db', {'required': False}))

    groupB = parser.add_argument_group('MANDATORY INPUT OPTION #2', "Or you can initiate this with an external genomes file.")
    groupB.add_argument(*anvio.A('external-genomes'), **anvio.K('external-genomes'))
    groupB.add_argument(*anvio.A('internal-genomes'), **anvio.K('internal-genomes'))

    groupC = parser.add_argument_group('ADDITIONAL INPUT (OPTIONAL)', "You can also give this program an anvi'o profile database along\
                                                         with a collection name. In which case anvi'o will estimate the completion and\
                                                         redundancy of every bin in this collection. Fun.")
    groupC.add_argument(*anvio.A('profile-db'), **anvio.K('profile-db', {'required': False}))
    groupC.add_argument(*anvio.A('collection-name'), **anvio.K('collection-name'))

    groupD = parser.add_argument_group('PARAMETERS OF CONVENIENCE', "Because life is already very hard as it is.")
    groupD.add_argument(*anvio.A('list-collections'), **anvio.K('list-collections'))
    groupD.add_argument(*anvio.A('report-all-estimates'), **anvio.K('report-all-estimates'))
    groupD.add_argument(*anvio.A('just-do-it'), **anvio.K('just-do-it'))
    groupD.add_argument(*anvio.A('concise'), **anvio.K('concise'))
    groupD.add_argument(*anvio.A('output-file'), **anvio.K('output-file'))

    return parser.get_args(parser)


if __name__ == '__main__':
    main()

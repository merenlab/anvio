#!/usr/bin/env python
# -*- coding: utf-8

import os
import sys
import copy

import anvio
import anvio.utils as utils
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths
import anvio.ccollections as ccollections

from anvio.learning import RF
from anvio.completeness import Completeness
from anvio.errors import ConfigError, FilesNPathsError
from anvio.hmmops import SequencesForHMMHits
from anvio.dbops import ContigsSuperclass


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = ["Tom O. Delmont"]
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['meren']
__description__ = "Screen for genomes to find likely members of CPR"



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
    progress = terminal.Progress()
    pp = terminal.pretty_print

    utils.is_contigs_db(args.contigs_db)

    collections = ccollections.Collections()
    collections.populate_collections_dict(args.contigs_db)
    if args.profile_db:
        collections.populate_collections_dict(args.profile_db)

    if args.list_collections:
        collections.list_collections()
        sys.exit()

    if args.output_file:
        filesnpaths.is_output_file_writable(args.output_file)

    if args.output_file and ( not args.just_do_it and filesnpaths.is_file_exists(args.output_file, dont_raise=True) ):
        raise ConfigError("The output file '%s' already exists. Please either remove it, or provide a different name. "
                          "Anvi'o does not like overwriting your stuff (well, only whenever it is convenient to do so)."\
                                    % args.output_file)

    if args.profile_db and not args.collection_name:
        raise ConfigError("You can't use this program with a profile database but without a collection name. Yes. Because.")

    if not args.collection_name:
        run.warning("You did not provide a collection name. Anvi'o will assume that what is in your contigs database "
                    "is a a single genome (or genome bin).")

    if args.collection_name and args.collection_name not in collections.collections_dict:
        raise ConfigError("%s is not a valid collection ID. See a list of available ones with '--list-collections' flag" % args.collection_name)

    completeness = Completeness(args.contigs_db)
    if not len(completeness.sources):
        raise ConfigError("HMM's were not run for this contigs database :/")

    if not 'Campbell_et_al' in completeness.sources:
        raise ConfigError("This classifier uses Campbell et al. single-copy gene collections, and it is not among the available HMM sources in your "
                           "contigs database :/ Bad news.")

    contigs_db = ContigsSuperclass(args, r = run, p = progress)

    if args.collection_name:
        collection = collections.get_collection_dict(args.collection_name)
    else:
        collection = {os.path.basename(args.contigs_db[:-3]): list(contigs_db.splits_basic_info.keys())}

    s = SequencesForHMMHits(args.contigs_db)
    hmm_hits_per_bin = s.get_hmm_hits_per_bin(collection, source = 'Campbell_et_al')

    rf = RF(args.classifier_object)
    rf.initialize_classifier()

    template = dict(list(zip(rf.features, [0] * len(rf.features))))
    data_dict = {}
    for bin_id in hmm_hits_per_bin:
        data_dict[bin_id] = copy.deepcopy(template)

        for feature in hmm_hits_per_bin[bin_id]:
            data_dict[bin_id][feature] = 1

    predictions = rf.predict(data_dict)

    if args.collection_name:
        run.warning(None, header = 'Bins in collection "%s"' % args.collection_name, lc = 'green')
    else:
        run.warning(None, header = 'Genome in "%s"' % os.path.basename(args.contigs_db), lc = 'green')
    for bin_id in predictions:
        class_prediction = 'CPR' if predictions[bin_id]['CPR'] > predictions[bin_id]['NON-CPR'] else 'NOT CPR'
        class_probability = max(predictions[bin_id].values()) * 100

        p_completion, p_redundancy, domain, domain_probabilities, info_text, results_dict = completeness.get_info_for_splits(set(collection[bin_id]))

        c = results_dict['bacteria']['Campbell_et_al']
        percent_completion = c['percent_completion']
        percent_redundancy = c['percent_redundancy']
        total_length = sum([contigs_db.splits_basic_info[split_name]['length'] for split_name in set(collection[bin_id])])

        if class_probability < args.min_class_probability or \
           percent_completion < args.min_percent_completion or \
           percent_redundancy > args.max_percent_redundancy or \
           total_length < (args.min_genome_size * 10**6):
            class_prediction = 'INCONCLUSIVE'

        if args.report_only_cpr and class_prediction != 'CPR':
            continue
        else:
            result = "%s (Conf: %.0f%%, Size: %s, C/R: %.0f%%/%.0f%%)" % \
                        (class_prediction, class_probability, pp(total_length), percent_completion, percent_redundancy)
            run.info(bin_id, result)
            if args.output_file:
                with open(args.output_file, "a") as f: f.write(result + '\n')


def get_args():
    from anvio.argparse import ArgumentParser

    parser = ArgumentParser(description=__description__)

    parser.add_argument('classifier_object', help = 'Model output generated by anvi-script-gen-CPR-classifier',
                        metavar='CLASSIFIER_FILE')
    parser.add_argument(*anvio.A('contigs-db'), **anvio.K('contigs-db'))
    parser.add_argument(*anvio.A('profile-db'), **anvio.K('profile-db', {'required': False}))
    parser.add_argument(*anvio.A('collection-name'), **anvio.K('collection-name'))
    parser.add_argument(*anvio.A('list-collections'), **anvio.K('list-collections'))
    parser.add_argument('--report-only-cpr', default = False, action = 'store_true', help = 'Include only bins that\
                        look like CPR genomes.')
    parser.add_argument('--min-genome-size', type = float, default = 0.5, help = 'Minimum genome size to consider for CPR in Mbp.\
                        Default is %(default)f')
    parser.add_argument('--min-percent-completion', type = int, default = 50, help = "Minimum percent completion estimate based on\
                        anvi'o default single-copy gene collections. Default is %(default)d")
    parser.add_argument('--max-percent-redundancy', type = int, default = 30, help = "Maxumum percent redundancy or single-copy genes\
                        in an anvi'o bin, or a genome to consider for classification. The default is %(default)d")
    parser.add_argument('--min-class-probability', type = int, default = 75, help = "If the classification confidence is below this\
                        don't bother. Default is %(default)d.")
    parser.add_argument(*anvio.A('output-file'), **anvio.K('output-file'))
    parser.add_argument(*anvio.A('just-do-it'), **anvio.K('just-do-it'))

    return parser.get_args(parser)


if __name__ == '__main__':
    main()

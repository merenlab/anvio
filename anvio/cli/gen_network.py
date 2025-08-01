#!/usr/bin/env python
# -*- coding: utf-8

import os
import sys

import anvio
import anvio.utils as utils
import anvio.dbops as dbops
import anvio.terminal as terminal

from anvio.errors import ConfigError, FilesNPathsError
from anvio.utils import gen_gexf_network_file


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['meren']
__description__ = "Generate a Gephi network for functions based on non-normalized gene coverage values"


class NetworkDescriptonSamples:
    def __init__(self, args, skip_auto_init=False):
        self.args = args
        self.run = terminal.Run()
        self.progress = terminal.Progress()

        self.run.warning("This is not a mature anvi'o program. Please don't use it for anything serious without "
                         "making sure you are on top of its limitations. If your current research question requires "
                         "you to generate network descriptions for your data, please let us know so maybe we can "
                         "work together to develop this tool towards a meaningful direction.")

        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        self.profile_db_path = A('profile_db')
        self.contigs_db_path = A('contigs_db')
        self.annotation_source = A('annotation_source')

        if (not self.profile_db_path) or (not self.contigs_db_path):
            raise ConfigError("This only works if you provide contigs and profile database paths")

        if A('list_annotation_sources'):
            self.contigs_db = dbops.ContigsSuperclass(self.args)
            self.contigs_db.list_function_sources()
            sys.exit()

        if not self.annotation_source:
            raise ConfigError("You must use a single functional annotation source for this to work. If you "
                              "are not sure what is available in your contigs database, run this program with "
                              "the flag '--list-annotation-sources'.")

        self.functions = {'labels': {}}

        self.P = lambda x: os.path.join(os.path.dirname(self.profile_db_path), x)

        if not skip_auto_init:
            self.init()

    def init(self):
        utils.is_profile_db_and_contigs_db_compatible(self.profile_db_path, self.contigs_db_path)

        self.profile_db = dbops.ProfileSuperclass(self.args)
        self.contigs_db = dbops.ContigsSuperclass(self.args)

        if not self.profile_db.p_meta['merged'] or self.profile_db.p_meta['blank']:
            raise ConfigError("The profile database describes a single run. Current implementation of this "
                              "program restricts its use to merged runs. Sorry :/")

        self.profile_db.init_gene_level_coverage_stats_dicts()
        self.contigs_db.init_functions([self.annotation_source])

        self.samples = self.profile_db.p_meta['samples']
        self.genes = list(self.profile_db.gene_level_coverage_stats_dict.keys())

        for gene in self.genes:
            if gene in self.contigs_db.gene_function_calls_dict:
                self.functions['labels'][gene] = self.contigs_db.gene_function_calls_dict[gene][self.annotation_source][1]


        self.samples_dict = self.get_samples_dict()


    def generate_functions_network(self):
        genes_with_functions = sorted(self.functions['labels'].keys())
        self.progress.new('Processing')
        self.progress.update('Generating network description for %d genes w/functions across %d samples ... ' % (len(genes_with_functions), len(self.samples_dict)))

        network_desc_output_path = self.P('SAMPLE-FUNCTION-NETWORK.gexf')
        gen_gexf_network_file(genes_with_functions, self.samples_dict,
                              network_desc_output_path, unit_mapping_dict = self.functions)
        self.progress.end()
        self.run.info('Network for gene functions', network_desc_output_path)


    def get_samples_dict(self, attribute='mean_coverage'):
        samples_dict = {}
        for sample in self.samples:
            samples_dict[sample] = {}

        d = self.profile_db.gene_level_coverage_stats_dict

        for gene_callers_id in d:
            for sample_name in d[gene_callers_id]:
                samples_dict[sample_name][gene_callers_id] = d[gene_callers_id][sample_name][attribute]
        return samples_dict


def main():
    args = get_args()

    try:
        network = NetworkDescriptonSamples(args)
        network.generate_functions_network()
    except ConfigError as e:
        print(e)
        sys.exit(-1)
    except FilesNPathsError as e:
        print(e)
        sys.exit(-2)


def get_args():
    from anvio.argparse import ArgumentParser

    parser = ArgumentParser(description=__description__)

    parser.add_argument(*anvio.A('profile-db'), **anvio.K('profile-db'))
    parser.add_argument(*anvio.A('contigs-db'), **anvio.K('contigs-db'))
    parser.add_argument(*anvio.A('annotation-source'), **anvio.K('annotation-source'))
    parser.add_argument(*anvio.A('list-annotation-sources'), **anvio.K('list-annotation-sources'))

    return parser.get_args(parser)


if __name__ == '__main__':
    main()

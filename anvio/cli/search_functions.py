#!/usr/bin/env python
# -*- coding: utf-8
"""Searches contigs database for a given function"""

import sys

import anvio
import anvio.utils as utils
import anvio.terminal as terminal

from anvio.errors import ConfigError, FilesNPathsError
from anvio.dbops import ContigsSuperclass, PanSuperclass


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['meren']
__requires__ = ['contigs-db', 'genomes-storage-db']
__provides__ = ['functions-txt']
__description__ = ("Search functions in an anvi'o contigs database or genomes storage. Basically, this program "
                   "searches for one or more search terms you define in functional annotations of "
                   "genes in an anvi'o contigs database, and generates multiple reports. The "
                   "default report simply tells you which contigs contain genes with functions "
                   "matching to serach terms you used, useful for viewing in the interface. You "
                   "can also request a much more comprehensive report, which gives you anything "
                   "you might need to know for each hit and serach term")


class SearchResultReporter(object):
    def __init__(self, args):
        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        self.annotation_sources = A('annotation_sources')
        self.list_annotation_sources = A('list_annotation_sources')
        self.basic_report_path = A('output_file') or 'search_results.txt'
        self.full_report_path = A('full_report')
        self.include_sequences = A('include_sequences')
        self.case_sensitive = A('case_sensitive')
        self.exact_match = A('exact_match')
        self.verbose = A('verbose')
        self.args = args
        self.search_mode = None

        self.run = terminal.Run()

        self.contigs_db = A('contigs_db')
        self.genomes_storage = A('genomes_storage')
        self.pan_db = A('pan_db')

        # get search results will fill these
        self.all_item_names = None
        self.matching_item_names_dict = None
        self.verbose_output = None

        # straighten things out
        self.search_terms = [s.strip() for s in A('search_terms').split(A('delimiter'))]

        self.db = self.get_database()

        if self.list_annotation_sources:
            self.db.list_function_sources()
            sys.exit()


    def process(self):
        self.get_search_results()
        self.write_report()

        if self.full_report_path:
            self.write_full_report()


    def get_database(self):
        if self.contigs_db and (self.genomes_storage or self.pan_db):
            raise ConfigError("You can not provide both contigs database and genomes storage")

        if self.contigs_db:
            self.search_mode = 'contigs'
            self.run.info("Searching in Contigs Database", self.contigs_db)
            return ContigsSuperclass(self.args)

        elif self.genomes_storage and self.pan_db:
            self.search_mode = 'gene_clusters'
            self.run.info("Searching in Genomes Storage", self.genomes_storage)
            pan_database = PanSuperclass(self.args)
            pan_database.init_gene_clusters()
            return pan_database

        else:
            raise ConfigError("You did not provide enough arguments to initialize contigs database or pan database. "
                              "To initialize pan database you need to provide both genome storage and pan database.")


    def get_search_results(self):
        self.matching_item_names_dict, self.verbose_output = self.db.search_for_gene_functions(self.search_terms, requested_sources=self.annotation_sources, verbose=self.verbose, case_sensitive=self.case_sensitive, exact_match=self.exact_match)

        self.all_item_names = set([])

        for item_names in list(self.matching_item_names_dict.values()):
            self.all_item_names.update(item_names)


    def write_report(self):
        results_dict = {}

        for item_name in self.all_item_names:
            results_dict[item_name] = dict([(s + '_hits', '') for s in self.search_terms])

            for search_term in self.search_terms:
                if item_name in self.matching_item_names_dict[search_term]:
                    results_dict[item_name][search_term + '_hits'] = search_term

        utils.store_dict_as_TAB_delimited_file(results_dict, self.basic_report_path, headers = [self.search_mode] + [s + '_hits' for s in self.search_terms])
        self.run.info('Items additional data compatible output', self.basic_report_path, nl_before=1)


    def write_full_report(self):
        if self.search_mode == 'contigs':
            header = ['gene_callers_id']
        elif self.search_mode == 'gene_clusters':
            header = ['gene_callers_id', 'genome_name']
        else:
            raise ConfigError("You ended up in a place you should have never ended up. Go back. Go back.")

        header.extend(['source', 'accession', 'function', 'search_term', self.search_mode])

        if self.include_sequences:
            if self.search_mode == 'contigs':
                gene_caller_ids = list(set([e[0] for e in self.verbose_output]))
                _, gene_sequences_dict = self.db.get_sequences_for_gene_callers_ids(gene_caller_ids, include_aa_sequences=True)
                header.extend(['direction', 'rev_compd', 'dna_sequence', 'aa_sequence'])
            elif self.search_mode == 'gene_clusters':
                header.extend(['dna_sequence', 'aa_sequence'])

        report = open(self.full_report_path, 'w')
        report.write('\t'.join(header) + '\n')
        for entry in self.verbose_output:
            content = [str(item) for item in entry]
            if self.include_sequences:
                if self.search_mode == 'contigs':
                    g = gene_sequences_dict[entry[0]]
                    content.extend([g['direction'],
                                    g['rev_compd'],
                                    g['sequence'],
                                    g['aa_sequence'],])
                elif self.search_mode == 'gene_clusters':
                    # pan results already contains dna and aa sequences
                    pass
            else:
                if self.search_mode == 'gene_clusters':
                    content = content[:-2]

            report.write('\t'.join(content) + '\n')
        report.close()

        self.run.info('Full report', self.full_report_path)


def main():
    args = get_args()

    try:
        SearchResultReporter(args).process()
    except ConfigError as e:
        print(e)
        sys.exit(-1)
    except FilesNPathsError as e:
        print(e)
        sys.exit(-2)


def get_args():
    from anvio.argparse import ArgumentParser

    parser = ArgumentParser(description=__description__)

    groupA = parser.add_argument_group('SEARCH IN', 'Relevant source databases')
    groupA.add_argument(*anvio.A('contigs-db'), **anvio.K('contigs-db', {'required': False}))
    groupA.add_argument(*anvio.A('pan-db'), **anvio.K('pan-db', {'required': False}))
    groupA.add_argument(*anvio.A('genomes-storage'), **anvio.K('genomes-storage', {'required': False}))

    groupB = parser.add_argument_group('SEARCH FOR', 'Relevant terms')
    groupB.add_argument(*anvio.A('search-terms'), **anvio.K('search-terms', {'required': True}))
    groupB.add_argument(*anvio.A('case-sensitive'), **anvio.K('case-sensitive'))
    groupB.add_argument(*anvio.A('exact-match'), **anvio.K('exact-match'))
    groupB.add_argument(*anvio.A('delimiter'), **anvio.K('delimiter'))
    groupB.add_argument(*anvio.A('annotation-sources'), **anvio.K('annotation-sources'))
    groupB.add_argument(*anvio.A('list-annotation-sources'), **anvio.K('list-annotation-sources'))

    groupC = parser.add_argument_group('REPORT', "Anvi'o can report the hits in multiple ways. The output file will be a very simple 2-column\
                                                  TAB-delimited output that is compatible with anvi'o additional data format (so you can give\
                                                  it to the `anvi-interactive` to see which splits contained genes that were matching to your\
                                                  search terms). You can also ask anvi'o to generate a full-report, that contains much more and\
                                                  much helpful information about each hit. Optionally you can even ask the gene sequences to\
                                                  appear in this report.")
    groupC.add_argument(*anvio.A('output-file'), **anvio.K('output-file'))
    groupC.add_argument(*anvio.A('full-report'), **anvio.K('full-report'))
    groupC.add_argument(*anvio.A('include-sequences'), **anvio.K('include-sequences'))
    groupC.add_argument(*anvio.A('verbose'), **anvio.K('verbose'))

    return parser.get_args(parser)


if __name__ == '__main__':
    main()

#!/usr/bin/env python
# -*- coding: utf-8

import re
import sys
import argparse
from anvio.argparse import ArgumentParser

import anvio
import anvio.dbops as dbops
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths

from anvio.dbinfo import DBInfo as dbi
from anvio.errors import ConfigError, FilesNPathsError
from anvio.tables.miscdata import TableForItemAdditionalData
from anvio.dbinfo import is_profile_db_and_contigs_db_compatible
from anvio.utils.database import get_all_item_names_from_the_database
from anvio.utils.files import store_dict_as_TAB_delimited_file
from anvio.utils.sequences import rev_comp


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['meren']
__requires__ = ["profile-db", "contigs-db", "genes-db"]
__provides__ = ['misc-data-items', 'misc-data-layers']
__description__ = "A program to find one or more sequence motifs in contig or gene sequences, and store their frequencies"


class SequenceMotifSearch:
    """A class to work with sequence motifs"""
    def __init__(self, args, run=terminal.Run(), progress=terminal.Progress(), only_get_instance=False):
        self.args = args
        self.run = run
        self.progress = progress
        self.P = terminal.pluralize

        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        self.contigs_db_path = A("contigs_db")
        self.profile_db_path = A("profile_db")
        self.genes_db_path = A("genes_db")
        self.output_file_path = A("output_file")
        self.store_in_db = A("store_in_db")
        self.motifs = set([s.strip().upper() for s in A("motifs").split(',') if len(s.strip())]) if A("motifs") else set([])

        self.sanity_check()

        self.run.info('Contigs database', self.contigs_db_path)
        self.run.info('Profile database', self.profile_db_path)
        self.run.info('Genes database', self.genes_db_path)
        self.run.info(f'Motifs ({len(self.motifs)} found)', ', '.join(self.motifs), mc='green')

        # this will be filled by the class member functions
        self.motif_frequencies = {}

        if not only_get_instance:
            self.search_motifs()
            self.store_motif_frequencies()


    def sanity_check(self):
        if not self.contigs_db_path:
            raise ConfigError("You must provide a contigs database to this program.")

        dbi(self.contigs_db_path, expecting='contigs')
        dbi(self.profile_db_path, expecting='profile') if self.profile_db_path else None
        dbi(self.genes_db_path, expecting='genes') if self.genes_db_path else None

        if self.profile_db_path and self.genes_db_path:
            raise ConfigError("You should use this program with either a profile database, or a "
                              "genes database, and not both :/")

        if not (self.profile_db_path or self.genes_db_path) and self.store_in_db:
            raise ConfigError("The flag `--store-in-db` is only relevant if you have provide a profile or "
                              "a genes database. If all you have is a contigs database, then all you will "
                              "get is a flat TAB-delimited output file for your frequencies.")

        if not (self.profile_db_path or self.genes_db_path) and not self.output_file_path:
            raise ConfigError("If all you have is a contigs database, then you must at least use the "
                              "`--output-file` parameter so anvi'o would have a place to store all the "
                              "motif frequencies.")

        if self.output_file_path:
            filesnpaths.is_output_file_writable(self.output_file_path)

        if self.profile_db_path:
            is_profile_db_and_contigs_db_compatible(self.profile_db_path, self.contigs_db_path)

        # let's make sure motif sequences are actual DNA/RNA sequences
        character_regex = re.compile(r'^[ACGTNUactgnu]+$')
        bad_motifs = [m for m in self.motifs if not bool(character_regex.search(m))]
        if len(bad_motifs):
            bad_motifs_str = ', '.join([f"'{m}'" for m in bad_motifs])
            bad_motifs_str += (". Please make sure your motif sequences are like normal people motif sequences "
                               "where they are a combination of the following characters: A, C, G, T, N and/or U")
            if len(self.motifs) == 1:
                raise ConfigError(f"Ehem. Anvi'o doesn't like your motif sequence {bad_motifs_str}.")
            else:
                raise ConfigError(f"Well well well. Anvi'o found {self.P('motif', len(bad_motifs))} among your "
                                  f"sequences that didn't look quite good: {bad_motifs_str}.")


    def search_motifs(self):
        """A wrapper function to search motifs

        Returns
        =======
        self.motif_frequencies : dict
            A dictionary of all motifs and their frequencies. Items of this dictionary will
            depend on the source of data. I.e., it will be gene caller ids if the user sent
            a genes database, and split names if the user sent a profile database as parameter.
            If the user did not mention a profile or genes database, then the search will be
            performed on contig sequences and the item names will be contig names.
        """

        if self.profile_db_path:
            self.search_motifs_in_split_sequences()
        elif self.genes_db_path:
            self.search_motifs_in_gene_sequences()
        else:
            self.search_motifs_in_contig_sequences()

        self.run.warning(None, header="MOTIF FREQUENCIES", lc='green')
        for motif in self.motifs:
            self.run.info(motif, sum([v[motif] for v in self.motif_frequencies.values()]))

        return self.motif_frequencies


    def get_motif_frequencies(self, sequences_dict):
        """Literally searches motifs in a sequences dictionary, and fills in `self.motif_frequencies`"""

        if not len(self.motifs):
            raise ConfigError("There are no motifs to search here. Weird.")

        self.motif_frequencies = {}
        for motif in self.motifs:
            for item_name in sequences_dict:
                if item_name not in self.motif_frequencies:
                    self.motif_frequencies[item_name] = {}

                self.motif_frequencies[item_name][motif] = sequences_dict[item_name]['sequence'].count(motif) + \
                                                           sequences_dict[item_name]['sequence'].count(rev_comp(motif))

    def store_motif_frequencies(self):
        """Stores motif frequencies into databases or output files

        If the user have profile and genes databases and asked results to be stored in
        databse, this funciton will take care of that.
        """
        if not len(self.motif_frequencies):
            raise ConfigError("Something is wrong, the motif frequencies dictionary is empty :/")

        args = None
        header_item_name = None

        if self.profile_db_path:
            header_item_name = "split_name"
            args = argparse.Namespace(profile_db=self.profile_db_path)
        elif self.genes_db_path:
            header_item_name="gene_callers_id"
            args = argparse.Namespace(genes_db=self.genes_db_path)
        else:
            header_item_name = "contig_name"

        if (self.profile_db_path or self.genes_db_path):
            if self.store_in_db:
                # here we will get a copy of the motif frequencies dict so we can store them in the
                # db with the `motif_` prefix:
                M = {}
                for item_name in self.motif_frequencies:
                    M[item_name] = {}
                    for motif in self.motifs:
                        M[item_name][f"motif_{motif}"] = self.motif_frequencies[item_name][motif]

                args.skip_check_names = True
                args.just_do_it = True
                TableForItemAdditionalData(args, r=terminal.Run(verbose=False)).add(M, [f'motif_{m}' for m in list(self.motifs)])
                self.run.info("Motif frequencies stored in db", "YES", mc='green', nl_before=1)
            else:
                self.run.info("Motif frequencies stored in db", "No (but could've been with `--store-in-db`)", nl_before=1)

        if self.output_file_path:
            store_dict_as_TAB_delimited_file(self.motif_frequencies, self.output_file_path, headers=[header_item_name] + sorted(list(self.motifs)))
            self.run.info("Motif frequencies stored in output file", self.output_file_path, mc='green')
        else:
            self.run.info("Motif frequencies stored in output file", "No (but could've been with `--output-file`)")


    def search_motifs_in_contig_sequences(self):
        contigs_db = dbops.ContigsSuperclass(self.args, r=self.run, p=self.progress)
        contigs_db.init_contig_sequences()

        return self.get_motif_frequencies(contigs_db.contig_sequences)


    def search_motifs_in_split_sequences(self):
        split_names_in_profile_db = get_all_item_names_from_the_database(self.profile_db_path, run=self.run)

        contigs_db = dbops.ContigsSuperclass(self.args, r=self.run, p=self.progress)
        contigs_db.init_split_sequences(split_names_of_interest=split_names_in_profile_db)

        return self.get_motif_frequencies(contigs_db.split_sequences)


    def search_motifs_in_gene_sequences(self):
        gene_caller_ids = get_all_item_names_from_the_database(self.genes_db_path)

        contigs_db = dbops.ContigsSuperclass(self.args, r=self.run, p=self.progress)
        gene_caller_ids_list, sequences_dict = contigs_db.get_sequences_for_gene_callers_ids(gene_caller_ids_list=list(gene_caller_ids))

        return self.get_motif_frequencies(sequences_dict)


def main():
    args, unknown = get_args()

    try:
        SequenceMotifSearch(args)
    except ConfigError as e:
        print(e)
        sys.exit(-1)
    except FilesNPathsError as e:
        print(e)
        sys.exit(-2)


def get_args():
    parser = ArgumentParser(description=__description__)

    groupA = parser.add_argument_group('SEQUENCES', "A contigs database, essentially, to search for motifs.")
    groupA.add_argument(*anvio.A('contigs-db'), **anvio.K('contigs-db'))

    groupB = parser.add_argument_group('OPTIONAL DBs', "This program can store the frequencies of your motifs into profile "
                                        "or gene databases. See the online documentation at the end of the help menu for details.")
    groupB.add_argument(*anvio.A('profile-db'), **anvio.K('profile-db', {'required': False}))
    groupB.add_argument(*anvio.A('genes-db'), **anvio.K('genes-db', {'required': False}))

    groupC = parser.add_argument_group('MOTIFS', "Sequences to search for..")
    groupC.add_argument('--motifs', required=True, help="The motif sequence. You can search for more than one, in which case"
                                        "you should use comma (',') to separate them from each other.")

    groupD = parser.add_argument_group('OUTPUT', "Output options. The output file is the obvious option. But if you provided a profile "
                                        "or genes database, AND if you use the flag `--store-in-db`, then anvi'o will also store the "
                                        "motif frequencies in your databases as items additional data")
    groupD.add_argument(*anvio.A('output-file'), **anvio.K('output-file'))
    groupD.add_argument(*anvio.A('store-in-db'), **anvio.K('store-in-db'))

    return parser.parse_known_args()


if __name__ == "__main__":
    main()

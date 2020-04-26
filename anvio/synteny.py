#!/usr/bin/env python # -*- coding: utf-8
# pylint: disable=line-too-long

"""Classes to work with ngrams of contig functions.

These are classes to deconstruct loci into ngrams. They are used
to analyze conserved genes and synteny structures across loci.
"""

import sys
import pandas as pd

from collections import Counter

import anvio
import anvio.dbops as dbops
import anvio.utils as utils
import anvio.tables as t
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths
import anvio.ccollections as ccollections
## FIXME: Will need to change if accepting genome-storage instead of external-genomes
import anvio.genomestorage as genomestorage
import anvio.genomedescriptions as genomedescriptions

from anvio.dbops import PanDatabase
from anvio.errors import ConfigError



__author__ = "Developers of anvi'o (see AUTHORS.txt)"
__copyright__ = "Copyleft 2015-2018, the Meren Lab (http://merenlab.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "Matthew Schechter"
__email__ = "mschechter@uchicago.edu"


class NGram(object):
    """class for counting NGrams

    anvi-analyze-synteny is designed to work with a group of similar loci, where each locus is a contig (which
    can lie in any number of contigs dbs)

    Parameters
    ==========
    args : argparse.Namespace
        For examples, arguments accepted by anvi-analyze-syntenty

    skip_sanity_check : bool, False
        If True, sanity_check will not be called.

    Notes
    =====
    - Currently the design assumes that each locus is a contig. In the future we have plans to expand this to
      compare genomes, not just loci, to each other. If that behavior is desired in the current design, each
      genome should be a single contig.
    """

    def __init__(self, args, run=terminal.Run(), progress=terminal.Progress(), skip_sanity_check=False):
        """Parses arguments and run sanity_check"""

        self.args = args
        self.run = run
        self.progress = progress

        # Parse arguments
        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        self.annotation_source = A('annotation_source')
        self.window_range = A('ngram_window_range') or "2:3"
        self.is_in_unknowns_mode = A('analyze_unknown_functions')
        self.output_file = A('output_file')
        self.skip_init_functions = A('skip_init_functions')


        self.pan_db_path = A('pan_db')
        if self.pan_db_path:
            self.pan_db = PanDatabase(self.pan_db_path)

            self.p_meta = self.pan_db.meta

            self.p_meta['creation_date'] = utils.get_time_to_date(self.p_meta['creation_date']) if 'creation_date' in self.p_meta else 'unknown'
            self.p_meta['genome_names'] = sorted([s.strip() for s in self.p_meta['external_genome_names'].split(',') + self.p_meta['internal_genome_names'].split(',') if s])
            self.p_meta['num_genomes'] = len(self.p_meta['genome_names'])
            self.genome_names = self.p_meta['genome_names']
            self.gene_clusters_gene_alignments_available = self.p_meta['gene_alignments_computed']
        else:
            self.pan_db = None

        self.genomes_storage_path = A('genomes_storage')

        if self.pan_db:
            self.genomes_storage = genomestorage.GenomeStorage(self.genomes_storage_path,
                                                               self.p_meta['genomes_storage_hash'],
                                                               genome_names_to_focus=self.p_meta['genome_names'],
                                                               skip_init_functions=self.skip_init_functions,
                                                               run=self.run,
                                                               progress=self.progress)
        else:
            self.genomes_storage = genomestorage.GenomeStorage(self.genomes_storage_path,
                                                               skip_init_functions=self.skip_init_functions,
                                                               run=self.run,
                                                               progress=self.progress)

        self.list_annotation_sources = A('list_annotation_sources')

        self.gene_function_source_set = self.genomes_storage.db.get_table_as_dataframe('gene_function_calls').source.unique()

        if self.list_annotation_sources:
            self.run.info('Available functional annotation sources', ', '.join(self.gene_function_source_set))
            sys.exit()
 
        # This houses the ngrams' data
        self.ngram_attributes_list = []

        self.num_contigs_in_external_genomes_with_genes = len(self.genomes_storage.get_all_genome_names())

        if not skip_sanity_check:
            self.sanity_check()

        # unless we are in debug mode, let's keep things quiet.
        if anvio.DEBUG:
            self.run_object = terminal.Run()
        else:
            self.run_object = terminal.Run(verbose=False)


    def sanity_check(self):
        """Sanity_check will confirm input for NGram class"""

        # checking if the annotation source is common accross all contigs databases
        ## FIXME: Will need to change if accepting genome-storage instead of external-genomes
        if self.annotation_source not in self.gene_function_source_set:
            raise ConfigError("The annotation source you requested does not appear to be in all of "
                              "the contigs databases from the external-genomes file. "
                              "Please confirm your annotation-source and that all contigs databases have it :)")

        if not self.args.output_file:
            raise ConfigError("You should provide an output file name.")

        # checking window-range input
        if ":" not in self.window_range:
            raise ConfigError("anvi'o would love to slice and dice your loci, but the "
                              "Format of window_range must be x:y (e.g. Window sizes 2 to 4 would be denoted as: 2:4)")

        try:
            self.window_range = [int(n) for n in self.window_range.split(":")]
        except ValueError:
            raise ConfigError("anvi'o would love to slice and dice your loci, but the "
                              "window-ranges need to be integers :)")

        if self.window_range[0] > self.window_range[1]:
            raise ConfigError("anvi'o would love to slice and dice your loci, but the "
                              "window-range needs to be from small to big :)")

        # Must contain 2 integers for window
        if len(self.window_range) > 2 or not isinstance(self.window_range[0], int) or not isinstance(self.window_range[1], int):
            raise ConfigError("anvi'o would love to slice and dice your loci, but... the "
                              "window_range must only contain 2 integers and be formated as x:y (e.g. Window sizes 2 to 4 would be denoted as: 2:4)")

        # Loop through each contigs db, test that each contig contains at least as many genes as max window size and confirm every contig has annotations
        for contigs_db_name in self.genomes_storage.get_genomes_dict():

            gene_caller_ids = self.genomes_storage.get_gene_caller_ids(contigs_db_name)
            num_genes = len(gene_caller_ids)

            if self.window_range[1] > num_genes:
                print("asdf")
                raise ConfigError("The largest window size you requested (%d) is larger than the number of genes found on this genome: %s" % \
                    (self.window_range[1], contigs_db_name))

    def populate_genes(self):
        """Iterates through all contigs and use self.count_synteny to count all ngrams in that contig.

        This populates the self.ngram_attributes_list, where each element looks like:
        (ngram, count, contigs_db_name, contig_name, n)

        """
        # Get gene cluster info from panDB
        if self.pan_db:
            gene_cluster_frequencies_dataframe = self.pan_db.db.get_table_as_dataframe('gene_clusters')

        self.run.warning("Anvi'o is now looking for Ngrams in your contigs!", lc='green')

        self.run.info_single("What do we say to loci that appear to have no coherent synteny patterns...? Not today! ⚔️", nl_before=1, nl_after=1)

        genes_and_functions_list = []
        for contigs_db_name in self.genomes_storage.get_genomes_dict():

            # Get list of genes-callers-ids
            gene_caller_ids_list = list(self.genomes_storage.get_gene_caller_ids(contigs_db_name))

            # Create dicts for annotate Ngrams
            gene_function_call_df = self.genomes_storage.db.get_table_as_dataframe('gene_function_calls')
            gene_caller_id_to_function_dict = self.get_genes_and_functions_dict(contigs_db_name, gene_function_call_df)

            if self.pan_db:
                gene_caller_id_to_gene_cluster_dict = self.get_gene_cluster_dict(contigs_db_name, gene_cluster_frequencies_dataframe)

            # Iterate over range of window sizes and run synteny algorithm to count occurrences of ngrams in a contig
            for n in range(self.window_range[0], self.window_range[1]):
                ngram_counts_dict = self.count_synteny(gene_caller_ids_list, n)

                # add results of this window size, contig pairing to ngram attributes
                for ngram, count in ngram_counts_dict.items():
                    if self.pan_db and self.annotation_source:
                        ngram_annotation = tuple([gene_caller_id_to_function_dict[g] for g in ngram])
                        ngram_gene_clusters = tuple([gene_caller_id_to_gene_cluster_dict[g] for g in ngram])
                        self.ngram_attributes_list.append([ngram, ngram_annotation, ngram_gene_clusters, count, contigs_db_name, n])
                    elif self.annotation_source and not self.pan_db:
                        ngram_annotation = tuple([gene_caller_id_to_function_dict[g] for g in ngram])
                        self.ngram_attributes_list.append([ngram, ngram_annotation, count, contigs_db_name, n]) 

    def convert_to_df(self):
        """Takes self.ngram_attributes_list and returns a pandas dataframe"""

        ngram_count_df_list = []

        for ngram_attribute in self.ngram_attributes_list:
            ngram = "::".join(map(str, list(ngram_attribute[0])))
            ngram_function = "::".join(map(str, list(ngram_attribute[1])))
            if self.pan_db and self.annotation_source:
                ngram_gene_clusters = "::".join(map(str, list(ngram_attribute[2])))
                df = pd.DataFrame(columns=['ngram', 'ngram_functions', 'ngram_gene_clusters', 'count', 'contig_db_name', 'N', 'number_of_loci'])
                df = df.append({'ngram': ngram,
                                'ngram_functions': ngram_function,
                                'ngram_gene_clusters': ngram_gene_clusters,
                                'count': ngram_attribute[3],
                                'contig_db_name': ngram_attribute[4],
                                'N':ngram_attribute[5],
                                'number_of_loci':self.num_contigs_in_external_genomes_with_genes}, ignore_index=True)
            elif not self.pan_db and self.annotation_source:
                df = pd.DataFrame(columns=['ngram', 'ngram_functions', 'count', 'contig_db_name', 'N', 'number_of_loci'])
                df = df.append({'ngram': ngram,
                                'ngram_functions': ngram_function,
                                'count': ngram_attribute[2],
                                'contig_db_name': ngram_attribute[3],
                                'N':ngram_attribute[4],
                                'number_of_loci':self.num_contigs_in_external_genomes_with_genes}, ignore_index=True)
            ngram_count_df_list.append(df)

        ngram_count_df_final = pd.concat(ngram_count_df_list)

        return ngram_count_df_final


    def report_ngrams_to_user(self):
        """Counts ngrams per contig and reports as tab-delimited file"""

        self.populate_genes()
        df = self.convert_to_df()
        df.to_csv(self.output_file, sep = '\t', index=False)
        self.run.info("Ngram table", self.output_file)


    def get_genes_and_functions_dict(self, contigs_db_name, gene_function_call_df):
        """This method will extract a list of gene attributes from each contig within a contigsDB.

        Returns
        =======
        output : list of lists
            first element is gene_caller_id, second is function accession, third is the contig name
        """

        # get contigsDB
        gene_function_call_df_filtered = gene_function_call_df[(gene_function_call_df['genome_name'] == contigs_db_name) & (gene_function_call_df['source'] == self.annotation_source)]
        gene_callers_id_to_accession_dict = gene_function_call_df_filtered[['gene_callers_id','accession']].set_index('gene_callers_id')['accession'].to_dict()
        
        gene_caller_ids_list = self.genomes_storage.get_gene_caller_ids(contigs_db_name)

        # Make list of lists containing gene attributes. If there is not annotation add one in!
        genes_and_functions_list = [] # List of lists [gene-caller-id, accessions, contig-name]
        counter = 0
        for gene_callers_id in gene_caller_ids_list: 
            list_of_gene_attributes = []

            if gene_callers_id in gene_callers_id_to_accession_dict:
                accession = gene_callers_id_to_accession_dict[gene_callers_id]
                accession = accession.replace(" ","")
                list_of_gene_attributes.extend((gene_callers_id, accession))
                genes_and_functions_list.append(list_of_gene_attributes)
            else:
                # adding in "unknown annotation" if there is none
                accession = "unknown-function"
                list_of_gene_attributes.extend((counter, accession))
                genes_and_functions_list.append(list_of_gene_attributes)
            counter = counter + 1

        gene_caller_id_to_accession_dict = {}
        for entry in genes_and_functions_list:
            gene_caller_id_to_accession_dict[entry[0]] = entry[1] 

        return gene_caller_id_to_accession_dict

    def get_gene_cluster_dict(self, contigs_db_name, gene_cluster_frequencies_dataframe):
        gene_cluster_frequencies_dataframe_filtered = gene_cluster_frequencies_dataframe[gene_cluster_frequencies_dataframe['genome_name'] == contigs_db_name]
        gene_callers_id_to_gene_cluster_id_dict = gene_cluster_frequencies_dataframe_filtered[['gene_caller_id','gene_cluster_id']].set_index('gene_caller_id')['gene_cluster_id'].to_dict()

        gene_caller_ids_list = self.genomes_storage.get_gene_caller_ids(contigs_db_name)

        # Make list of lists containing gene cluster attributes. If there is not annotation add one in!
        genes_cluster_list = [] # List of lists [gene-caller-id, gene-cluster-id, contig-name]
        counter = 0
        for gene_callers_id in gene_caller_ids_list: 
            list_of_gene_attributes = []

            if gene_callers_id in gene_callers_id_to_gene_cluster_id_dict:
                gene_cluster_id = gene_callers_id_to_gene_cluster_id_dict[gene_callers_id]
                gene_cluster_id = gene_cluster_id.replace(" ","")
                list_of_gene_attributes.extend((gene_callers_id, gene_cluster_id))
                genes_cluster_list.append(list_of_gene_attributes)
            else:
                # adding in "unknown annotation" if there is none
                gene_cluster_id = "no-gene-cluster-annotation"
                list_of_gene_attributes.extend((counter, gene_cluster_id))
                genes_cluster_list.append(list_of_gene_attributes)
            counter = counter + 1

        gene_caller_id_to_gene_cluster_dict = {}
        for entry in genes_cluster_list:
            gene_caller_id_to_gene_cluster_dict[entry[0]] = entry[1] 

        return gene_caller_id_to_gene_cluster_dict


    def count_synteny(self, function_list, n):
        """This method will count NGrams in contigs

        This method will interate through a dict of contigs {contig_name: genes_and_functions_list} count NGrams
        in each contig using a sliding window of size N. The final output will be a dictionary {ngram:count}

        Parameters
        ==========
        function_list : list
            A list of gene functions as they appear in the contig

        n : int
            A window size to extract a ngram

        Returns
        =======
        NGramFreq_dict : dict
            A dict where each key is list of genes that represent an ngram and the value is the count in that contig

        Notes
        =====
        Future goal: Need to return counts for 1 contig at a time and
        give back a dictionary with contig {contig_name: {ngram:count}}
        """

        ngram_frequencies_dict = Counter({})
        for i in range(0, len(function_list) - n + 1):
            # window = sorted(function_list[i:i + n])
            # extract window
            window = function_list[i:i + n]

            # order the window based arbitrarily based on the first gene in the synteny being smallest
            original_order = window
            flipped_order = window[::-1]
            if original_order[0] < flipped_order[0]:
                window = original_order
            else:
                window = flipped_order

            ngram = tuple(window)

            if not self.is_in_unknowns_mode and "unknown-function" in ngram: # conditional to record NGrams with unk functions
                continue
            else:
                # if ngram is not in dictionary add it, if it is add + 1
                ngram_frequencies_dict[ngram] +=  1

        return ngram_frequencies_dict


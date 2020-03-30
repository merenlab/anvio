#!/usr/bin/env python
# -*- coding: utf-8
# pylint: disable=line-too-long

"""Classes to work with ngrams of contig functions.

These are classes to deconstruct loci into ngrams. They are used
to analyze conserved genes and synteny structures across loci.
"""

import pandas as pd

from collections import Counter

import anvio
import anvio.dbops as dbops
import anvio.utils as utils
import anvio.tables as t
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths
import anvio.genomedescriptions as genomedescriptions

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
        self.in_in_unknowns_mode = A('analyze_unknown_functions')
        self.output_file = A('output_file')
        self.genomes = genomedescriptions.GenomeDescriptions(self.args)
        self.genomes.load_genomes_descriptions(init=False)

        # This houses the ngrams' data
        self.ngram_attributes_list = []

        self.num_contigs_in_external_genomes_with_genes = 0

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
        if self.annotation_source not in self.genomes.function_annotation_sources:
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
        # Set the self.num_contigs_in_external_genomes_with_genes variable
        for contigs_db_name in self.genomes.external_genomes_dict:
            # extract contigsDB path
            contigs_db_path = self.genomes.external_genomes_dict[contigs_db_name]["contigs_db_path"]

            # Confirm appropriate window size
            genes_and_functions_list = self.get_genes_and_functions_from_contigs_db(contigs_db_path)
            num_genes = len(genes_and_functions_list)

            if self.window_range[1] > num_genes:
                raise ConfigError("The largest window size you requested (%d) is larger than the number of genes found on this contig: %s" % \
                    (self.window_range[1], self.genomes.external_genomes_dict[contigs_db_name]["name"]))

            # Confirm all contigs have annotations
            contigs_db = dbops.ContigsDatabase(contigs_db_path)
            genes_in_contigs = contigs_db.db.get_table_as_dataframe('genes_in_contigs')
            all_contigs_in_db = contigs_db.db.get_table_as_dataframe('contigs_basic_info')['contig']

            self.num_contigs_in_external_genomes_with_genes += genes_in_contigs['contig'].nunique()

            contigs_without_genes = []
            for contig in all_contigs_in_db:
                if contig not in set(genes_in_contigs['contig'].unique()):
                    contigs_without_genes.append(contig)

            if len(contigs_without_genes):
                self.run.warning("Just so you know, %d contigs in %s had no genes. Here are the the first 5: %s" % \
                                  (len(contigs_without_genes), contigs_db_name, ''.join([str(x) for x in contigs_without_genes[:5]])))


    def populate_genes(self):
        """Iterates through all contigs and use self.count_synteny to count all ngrams in that contig.

        This populates the self.ngram_attributes_list, where each element looks like:
        (ngram, count, contigs_db_name, contig_name, n)

        """

        genes_and_functions_list = []
        for contigs_db_name in self.genomes.external_genomes_dict:
            # Extract file path
            contigs_db_path = self.genomes.external_genomes_dict[contigs_db_name]["contigs_db_path"]

            # Get list of genes and functions
            genes_and_functions_list = self.get_genes_and_functions_from_contigs_db(contigs_db_path)

            # Get unique list of the contigs from this contigsDB (there could be more than one)
            contigs_set = set(([entry[2] for entry in genes_and_functions_list]))

            for contig_name in contigs_set:
                contig_function_list = []

                # Extract gene_callers_id and gene_function_accession from genes_and_functions_list using contig_name
                # genes_and_functions_list has multiple contigs gene info in it, so we filter one contig at a time
                # (contig_ID = name of contig in genes_and_functions_list)
                for gene_callers_id, gene_function_accession, contig_ID in genes_and_functions_list:
                    if contig_name == contig_ID:
                        contig_function_list.append([gene_callers_id,gene_function_accession])

                # Iterate over range of window sizes and run synteny algorithm to count occurrences of ngrams in a contig
                for n in range(self.window_range[0],self.window_range[1]):
                    function_list = [entry[1] for entry in contig_function_list]
                    ngram_counts_dict = self.count_synteny(function_list, n)

                    # add results of this window size, contig pairing to ngram attributes
                    for ngram, count in ngram_counts_dict.items():
                        self.ngram_attributes_list.append([ngram, count, contigs_db_name, contig_name, n])


    def convert_to_df(self):
        """Takes self.ngram_attributes_list and returns a pandas dataframe"""

        ngram_count_df_list = []

        for ngram_attribute in self.ngram_attributes_list:
            ngram = "::".join(map(str, list(ngram_attribute[0])))
            df = pd.DataFrame(columns=['ngram', 'count', 'contig_db_name','contig_name','N','number_of_loci'])
            df = df.append({'ngram': ngram,
                            'count': ngram_attribute[1],
                            'contig_db_name': ngram_attribute[2],
                            'contig_name':ngram_attribute[3],
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


    def get_genes_and_functions_from_contigs_db(self, contigs_db_path):
        """This method will extract a list of gene attributes from each contig within a contigsDB.

        Returns
        =======
        output : list of lists
            first element is gene_caller_id, second is function accession, third is the contig name
        """

        # get contigsDB
        contigs_db = dbops.ContigsDatabase(contigs_db_path)
        
        # extract contigs names
        genes_in_contigs = contigs_db.db.get_table_as_dict(t.genes_in_contigs_table_name)
        
        # extract annotations and filter for the sources designated by user using self.annotation_source
        annotations_dict = contigs_db.db.get_table_as_dict(t.gene_function_calls_table_name)
        annotations_dict = utils.get_filtered_dict(annotations_dict, 'source', set([self.annotation_source]))

        # Make dict with gene-caller-id:accession
        gene_callers_id_to_accession_dict = {entry['gene_callers_id']: entry['accession'] for entry in annotations_dict.values()}

        # Make list of lists containing gene attributes. If there is not annotation add one in!
        genes_and_functions_list = [] # List of lists [gene-caller-id, accessions, contig-name]
        counter = 0
        for gene_callers_id in genes_in_contigs: 
            list_of_gene_attributes = []

            if gene_callers_id in gene_callers_id_to_accession_dict:
                accession = gene_callers_id_to_accession_dict[gene_callers_id]
                accession = accession.replace(" ","")
                contig_name = genes_in_contigs[gene_callers_id]['contig']
                list_of_gene_attributes.extend((gene_callers_id, accession, contig_name))
                genes_and_functions_list.append(list_of_gene_attributes)
            else:
                # adding in "unknown annotation" if there is none
                accession = "unknown-function"
                contig_name = genes_in_contigs[counter]['contig']
                list_of_gene_attributes.extend((counter,accession,contig_name))
                genes_and_functions_list.append(list_of_gene_attributes)
            counter = counter + 1

        return genes_and_functions_list


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

            if not self.in_in_unknowns_mode and "unknown-function" in ngram: # conditional to record NGrams with unk functions
                continue
            else:
                # if ngram is not in dictionary add it, if it is add + 1
                ngram_frequencies_dict[ngram] +=  1

        return ngram_frequencies_dict


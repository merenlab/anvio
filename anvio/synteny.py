#!/usr/bin/env python
# -*- coding: utf-8
# pylint: disable=line-too-long
"""
Classes to work with ngrams of contig functions.

These are classes to deconstruct loci into ngrams. They will be used
to analyze conserved genes and synteny structures across loci.
"""

import pandas as pd
from collections import Counter
import anvio
import anvio.tables as t
import anvio.utils as utils
import anvio.dbops as dbops
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
    def __init__(self, args, run=terminal.Run(), progress=terminal.Progress(), skip_sanity_check=False):
        """
        __init__ will parse arguments and run sanity_check

        """

        self.args = args
        self.run = run
        self.progress = progress

        # Parse arguments
        #----------------
        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        self.external_genomes = A('external_genomes')
        self.annotation_source = A('annotation_source')
        self.window_range = A('window_range')
        self.in_in_unknowns_mode = A('analyze_unknown_functions')
        self.external_genomes = utils.get_TAB_delimited_file_as_dictionary(self.external_genomes)
        self.external_genomes_number = len(self.external_genomes)
        self.output_file = A('output_file')

        # Run sanity_check
        #-----------------
        if not skip_sanity_check:
            self.sanity_check()

        # unless we are in debug mode, let's keep things quiet.
        if anvio.DEBUG:
            self.run_object = terminal.Run()
        else:
            self.run_object = terminal.Run(verbose=False)


    def sanity_check(self):
        """
        sanity_check will confirm input for NGram class
        """

        # checking if the annotation source is common accross all contigs databases
        #----------------------------

        g = genomedescriptions.GenomeDescriptions(self.args)
        g.load_genomes_descriptions(init=False)

        if self.annotation_source not in g.function_annotation_sources:
            raise ConfigError("The annotation source you requested does not appear to be in the "
                              "contigs database :/ "
                              "Please confirm you are only providing one annotation source :)")

        if not self.args.output_file:
            raise ConfigError("You should provide an output file name.")

        # checking window-range input
        #----------------------------

        # Must contain ":"
        if ":" not in self.window_range:
            raise ConfigError("anvi'o would love to slice and dice your loci, but... the "
                              "Format of window_range must be x:y (e.g. Window sizes 2 to 4 would be denoted as: 2:4)")

        # Must contain 2 integers for window
        self.window_range = [int(n) for n in self.window_range.split(":")] # parse self.window_range
        if self.window_range[0] > self.window_range[1]:
            raise ConfigError("anvi'o would love to slice and dice your loci, but... the "
                              "window-range needs to be from small to")
        if len(self.window_range) > 2 or not isinstance(self.window_range[0], int) or not isinstance(self.window_range[1], int):
            raise ConfigError("anvi'o would love to slice and dice your loci, but... the "
                              "window_range must only contain 2 integers and be formated as x:y (e.g. Window sizes 2 to 4 would be denoted as: 2:4)")

    def populate_genes(self):
        """
        populate_genes() will iterate through all contigs and use self.count_synteny to count all ngrams in that contig.

        The output of this method will be a list of ngram_attributes_list containing: ngram (ngram string), 
                                                                                      count (count of ngram in the contig), 
                                                                                      contigs_db_name (name of contigsDB), 
                                                                                      contig_name (nameof contig IN contigsDB), 
                                                                                      n (ngram size)
        """
        genes_and_functions_list = []
        ngram_attributes_list = []
        self.num_contigs_in_external_genomes = 0
        for contigs_db_name in self.external_genomes:
            # Extract file path
            contigs_db_path = self.external_genomes[contigs_db_name]["contigs_db_path"]

            # Get list of genes and functions
            genes_and_functions_list = self.get_genes_and_functions_from_contigs_db(contigs_db_path)

            # Get unique list of the contigs from this contigsDB (there could be more than one)
            contigs_list = set(([entry[2] for entry in genes_and_functions_list]))

            # Calculate TOTAL number of contigs within external-genomes files (there may be more than one per contigsDB)
            self.num_contigs_in_external_genomes = self.num_contigs_in_external_genomes + len(contigs_list)

            # iterate through list of contigs and make dictionary 'contig_name': list_of_functions
            contigs_dict = {}
            for contig_name in contigs_list:
                contig_function_list = []
                for gci_accession_contigname in genes_and_functions_list:
                    if contig_name == gci_accession_contigname[2]:
                        contig_function_list.append([gci_accession_contigname[0],gci_accession_contigname[1]])
                # sanity check to see that window size is smaller than loci length
                if len(contig_function_list) < self.window_range[1]:
                    raise ConfigError("anvi'o noticed that one of your ngram window sizes is larger than the number of genes in a locus! "
                      "please change your --window-range so that the largest n size is smaller than the smallest locus :)")
                contigs_dict[contig_name] = contig_function_list

            # Iterate over range of window sizes and run synteny algorithm to count occurrences of ngrams
                for n in range(self.window_range[0],self.window_range[1]):
                    ngram_count_df_list_dict = self.count_synteny(contigs_dict, n)
                    # make list of ngram attributes
                    for ngram,count in ngram_count_df_list_dict.items():
                        inidvidual_ngram_attributes_list = [ngram, count, contigs_db_name, contig_name, n]
                        ngram_attributes_list.append(inidvidual_ngram_attributes_list)
        return ngram_attributes_list


    def convert_to_df(self):
        """
        convert_to_df() will take the list of ngram_attributes_list from opulate_genes() and return a pandas dataframe
        """
        ngram_attributes_list = self.populate_genes()

        ngram_count_df_list = []

        for ngram_attribute in ngram_attributes_list:
            ngram = "::".join(map(str, list(ngram_attribute[0])))
            df = pd.DataFrame(columns=['ngram', 'count', 'contigDB','contig_name','N','number_of_loci'])
            df = df.append({'ngram': ngram, 
                            'count': ngram_attribute[1], 
                            'contigDB': ngram_attribute[2], 
                            'contig_name':ngram_attribute[3],
                            'N':ngram_attribute[4],
                            'number_of_loci':self.num_contigs_in_external_genomes}, ignore_index=True)
            ngram_count_df_list.append(df)

        ngram_count_df_final = pd.concat(ngram_count_df_list)

        return ngram_count_df_final

    def report_ngrams_to_user(self):
        """
        This method will save the pandas dataframe from convert_to_df for the user :)
        """
        df = self.convert_to_df()
        df.to_csv(self.output_file, sep = '\t',index=False)


    def get_genes_and_functions_from_contigs_db(self, contigs_db_path):
        """
        This method will extract a list of gene attributes from each contig within a contigsDB.

        The list will contain: gci (gene-caller-id), 
                               accession (gene function accession, e.g. COG1234), 
                               contig_name (name of contig)
        """
        # get contigsDB
        contigs_db = dbops.ContigsDatabase(contigs_db_path)
        # extract contigs names
        genes_in_contigs = contigs_db.db.get_table_as_dict(t.genes_in_contigs_table_name)
        # extract annotations and filter for the sources designated by user using self.annotation_source
        annotations_dict = contigs_db.db.get_table_as_dict(t.gene_function_calls_table_name)
        annotations_dict = utils.get_filtered_dict(annotations_dict, 'source', set([self.annotation_source]))

        # Make dict with gene-caller-id:accession
        gci_to_accession_dict = {entry['gene_callers_id']: entry['accession'] for entry in annotations_dict.values()}

        # Make list of lists containing gene attributes. If there is not annotation add one in!
        genes_and_functions_list = [] # List of lists [gene-caller-id, accessions, contig-name]
        counter = 0
        for gci in genes_in_contigs: # gci = gene-caller-id
            list_of_gene_attributes = []
            if gci in gci_to_accession_dict:
                accession = gci_to_accession_dict[gci]
                accession = accession.replace(" ","")
                contig_name = genes_in_contigs[gci]['contig']
                list_of_gene_attributes.extend((gci, accession, contig_name))
                genes_and_functions_list.append(list_of_gene_attributes)
            else: # adding in "unknown annotation" if there is none
                accession = "unknown-function"
                contig_name = genes_in_contigs[counter]['contig']
                list_of_gene_attributes.extend((counter,accession,contig_name))
                genes_and_functions_list.append(list_of_gene_attributes)
            counter = counter + 1

        return genes_and_functions_list


    def count_synteny(self, contigs_dict, n):
        """
        This method will interate through a dict of contigs {contig_name: genes_and_functions_list} count NGrams 
        in each contig using a sliding window of size N. The final output will be a dictionary {ngram:count} 

        CURRENTLY: count_synteny will break if there is more than one contig per contigsDB (i.e. NGramFreq_dict will
        be reset at the loop and only hold counts for the last contig in the contigsDB)

        Future goal:
        Need to return counts for 1 contig at a time and
        give back a dictionary with contig {contig_name: {ngram:count}}
        """
        for contig_name in sorted(contigs_dict.keys()):
            contig_gci_function_list = contigs_dict[contig_name]
            function_list = [entry[1] for entry in contig_gci_function_list]

            NGramFreq_dict = {}
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
                    if ngram in NGramFreq_dict:
                        NGramFreq_dict[ngram] +=  1
                    else:
                        NGramFreq_dict[ngram] = 1
            return NGramFreq_dict


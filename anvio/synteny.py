
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
        __init__ will parse arguments, run sanity_check, and run the driver method
        of this class, populate_genes.

        """

        self.args = args
        self.run = run
        self.progress = progress

        self.genes = {}
        self.contigs_list = {}

        # Parse arguments
        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        self.external_genomes = A('external_genomes')
        self.annotation_source = A('annotation_source')
        self.window_range = A('window_range')
        self.is_in_unkowns_mode = A('analyze_unknown_functions')
        self.external_genomes = utils.get_TAB_delimited_file_as_dictionary(self.external_genomes)
        self.output_file = A('output_file')

        # Run main methods
        if not skip_sanity_check:
            self.sanity_check()

        # unless we are in debug mode, let's keep things quiet.
        if anvio.DEBUG:
            self.run_object = terminal.Run()
        else:
            self.run_object = terminal.Run(verbose=False)


    def sanity_check(self):

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
        self.window_range = [int(n) for n in self.window_range.split(":")]
        if len(self.window_range) > 2 or not isinstance(self.window_range[0], int) or not isinstance(self.window_range[1], int):
            raise ConfigError("anvi'o would love to slice and dice your loci, but... the "
                              "window_range must only contain 2 integers and be formated as x:y (e.g. Window sizes 2 to 4 would be denoted as: 2:4)")

        # FIXME: add sanity check where config error is raised if the window size is larger than the loci length

    def populate_genes(self):
        genes_and_functions_list = []
        ngram_count_df_list = []
        ngram_count_df = pd.DataFrame(columns=['ngram', 'count', 'contigDB', 'contig_name', 'N'])
        final_list = []
        counter = 0
        # Iterate through contigsDBs
        for contigs_db_name in self.external_genomes:
            # Extract file path
            contigs_db_path = self.external_genomes[contigs_db_name]["contigs_db_path"]

            # Get list of genes and functions
            genes_and_functions_list = self.get_genes_and_functions_from_contigs_db(contigs_db_path)

            # Get unique list of the contigs from this contigsDB (there could be more than one)
            contigs_list = set(([entry[2] for entry in genes_and_functions_list]))

            # iterate through list of contigs and make dictionary 'contig_name': list_of_functions
            contigs_dict = {}
            for contig_name in contigs_list:
                contig_function_list = []
                for i in genes_and_functions_list:
                    if contig_name == i[2]:
                        contig_function_list.append([i[0],i[1]])
                contigs_dict[contig_name] = contig_function_list

            # Iterate over range of window sizes and run synteny algorithm to count occurrences of ngrams
                for n in range(self.window_range[0],self.window_range[1]):
                    ngram_count_df_list_dict = self.count_synteny(contigs_dict, n)
                    df = pd.DataFrame(list(ngram_count_df_list_dict.items()), columns= ['ngram','count'])
                    df['contigDB'] = contigs_db_name
                    df['contig_name'] = contig_name
                    df['N'] = n
                    ngram_count_df_list.append(df)

                # counter = counter + 1
                # if counter == 2:
                #     break
       
        ngram_count_df_final = pd.concat(ngram_count_df_list)
        ngram_count_df_final.to_csv(self.output_file, sep = '\t',index=False)


        # print(ngram_count_df_final)



    def get_genes_and_functions_from_contigs_db(self, contigs_db_path):
        # get contigsDB
        contigs_db = dbops.ContigsDatabase(contigs_db_path)
        # extract contigs names
        genes_in_contigs = contigs_db.db.get_table_as_dict(t.genes_in_contigs_table_name)
        # extract annotations and filter for the sources designated by user
        annotations_dict = contigs_db.db.get_table_as_dict(t.gene_function_calls_table_name)
        annotations_dict = utils.get_filtered_dict(annotations_dict, 'source', set([self.annotation_source]))

        # Make dict with gene-caller-id:accession
        gene_to_function_dict = {entry['gene_callers_id']: entry['accession'] for entry in annotations_dict.values()}

        # Make list of lists containing gene attributes. If there is not annotation add one in!
        genes_and_functions_list = [] # List of lists [gene-caller-id, accessions, contig-name]
        counter = 0
        for gci in genes_in_contigs:
            list_of_gene_attributes = []
            if gci in gene_to_function_dict:
                accession = gene_to_function_dict[gci]
                accession = accession.replace(" ","")
                contig_name = genes_in_contigs[gci]['contig']
                list_of_gene_attributes.extend((gci, accession, contig_name))
                genes_and_functions_list.append(list_of_gene_attributes)
            else: # adding in "unknown annotation" if there is none
                # accession = "unknown-function-" + "{:06d}".format(1) + str(counter) # add leading 0
                accession = "unknown-function"
                contig_name = genes_in_contigs[counter]['contig']
                list_of_gene_attributes.extend((counter,accession,contig_name))
                genes_and_functions_list.append(list_of_gene_attributes)
            counter = counter + 1

        return genes_and_functions_list


    def count_synteny(self, contigs_dict, n):
        """
        Need to return counts for 1 contig at a time and
        give back a dictionary with contig {name: {ngram:count}}
        """
        # if self.is_in_unkowns_mode:
        #     print("as;ldkjfas;lkdfas;dlkfa;slkdfl;ksd")
        # n = int(n)
        for contig_name in sorted(contigs_dict.keys()):
            contig_gci_function_list = contigs_dict[contig_name]
            function_list = [entry[1] for entry in contig_gci_function_list]

            kFreq = {}
            for i in range(0, len(function_list) - n + 1):
                window = sorted(function_list[i:i + n])
                ngram = "::".join(map(str, list(window)))
                if not self.is_in_unkowns_mode and "unknown-function" in ngram: # conditional to record unk functions
                    continue
                else:
                    # if ngram is not in dictionary add it
                    # if it is add + 1
                    if ngram in kFreq:
                        kFreq[ngram] +=  1
                    else:
                        kFreq[ngram] = 1
            return kFreq


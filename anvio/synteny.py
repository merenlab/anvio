# -*- coding: utf-8
# pylint: disable=line-too-long
"""
    Classes to work with ngrams.

    These are classes to deconstruct loci into ngrams. They will be used
    to analyze conserved genes and synteny structures across loci.
"""

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
    def __init__(self, args, run=terminal.Run(), progress=terminal.Progress()):
        """
        To test this code run:

            >>> ./run_pangenome_tests.sh
            >>> cd anvio/anvio/tests/sandbox/test-output/pan_test
            >>> anvi-analyze-synteny -e external-genomes.txt  \
                                     --annotation-source COG_FUNCTION \
                                     --window-size 3 \
                                     -o test
        """

        self.args = args
        self.run = run
        self.progress = progress

        self.genes = {}
        self.contigs_list = {}

        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        self.external_genomes = A('external_genomes')
        self.annotation_source = A('annotation_source')
        self.window_size = A('window_size')

        self.external_genomes = utils.get_TAB_delimited_file_as_dictionary(self.external_genomes)

        self.sanity_check()

        self.populate_genes()


    def sanity_check(self):
        #FIXME: do the remaining sanity checks here.

        # checking if the annotation source is common accross all contigs databases
        g = genomedescriptions.GenomeDescriptions(self.args)
        g.load_genomes_descriptions(init=False)

        if self.annotation_source not in g.function_annotation_sources:
            raise ConfigError("The annotation source you requested does not appear to be in the\
                               contigs database :/\
                               Please confirm you are only providing one annotation source :)")

        # Confirm only one annotation source
        # if len(self.annotation_source) > 1:
        #     raise ConfigError("You can only provide one annotation source  for this\
        #                         anvi-analyze-synteny.")

        if not self.args.output_file:
            raise ConfigError("You should provide an output file name.")


    def populate_genes(self):
        genes_and_functions_list = []
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
            for item in contigs_list:
                contig_function_list = []
                for i in genes_and_functions_list:
                    if item == i[2]:
                        contig_function_list.append([i[0],i[1]])
                contigs_dict[item] = contig_function_list

            # Run synteny algorithm and count occurrences of ngrams
            print(self.count_synteny(contigs_dict))

            counter = counter + 1
            if counter == 2:
                break

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
                accession =accession.split(", ")
                accession = sorted(accession)
                accession = "-".join(map(str, list(accession)))
                contig_name = genes_in_contigs[gci]['contig']
                list_of_gene_attributes.append(gci)
                list_of_gene_attributes.append(accession)
                list_of_gene_attributes.append(contig_name)
                genes_and_functions_list.append(list_of_gene_attributes)
            else: # adding in "unknown annotation" if there is none
                gci = counter
                accession = "unknown-function-" + str(counter)
                contig_name = genes_in_contigs[gci]['contig']
                list_of_gene_attributes.append(gci)
                list_of_gene_attributes.append(accession)
                list_of_gene_attributes.append(contig_name)
                genes_and_functions_list.append(list_of_gene_attributes)
            counter = counter + 1

        return genes_and_functions_list


    def count_synteny(self, contigs_dict):
        """
        Need to make an example where there is more than one contig in a contigsDB
        """
        k = self.window_size
        for key in sorted(contigs_dict.keys()):
            gci_function = contigs_dict[key]
            genes = [entry[1] for entry in gci_function]

            kFreq = {}
            # k = 3
            for i in range(0, len(genes) - k + 1):
                window = sorted(genes[i:i + k])
                ngram = "_".join(map(str, list(window)))
                # print(ngram)
                # if ngram is not in dictionary add it
                # if it is add + 1
                if ngram in kFreq:
                    kFreq[ngram] +=  1
                else:
                    kFreq[ngram] = 1
            return kFreq


        # for i in range(0, len(genes) - k + 1):
        #     print(i)
        #     # extract window
        #     window = sorted(genes[i:i + k])
        #     ngram = "_".join(map(str, list(window)))
        #     # if ngram is not in dictionary add it
        #     # if it is add + 1
        #     if ngram in kFreq:
        #         kFreq[ngram] +=  1
        #     else:
        #         kFreq[ngram] = 1
        #     return kFreq
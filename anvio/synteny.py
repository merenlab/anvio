# -*- coding: utf-8
# pylint: disable=line-too-long
"""
    Classes to work with ngrams.
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
            >>> anvi-analyze-synteny -e external-genomes.txt  --annotation-source COG_CATEGORY
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
                               contigs database :/")


    def populate_genes(self):
        genes_and_functions_list = []
        counter = 0
        for contigs_db_name in self.external_genomes:
            contigs_db_path = self.external_genomes[contigs_db_name]["contigs_db_path"]

            # contigs_list = [] 
            contigs_dict = {}


            genes_and_functions_list = self.get_genes_and_functions_from_contigs_db(contigs_db_path)
            # print(type(genes_and_functions_list))
            # contigs_list.append([entry[2] for entry in genes_and_functions_list])  
            # Get unique list of the contigs from this contigsDB (there could be more than one)
            contigs_list = set(([entry[2] for entry in genes_and_functions_list]))
            for item in contigs_list:
                contig_function_list = []
                for i in genes_and_functions_list:
                    if item == i[2]:
                        contig_function_list.append([i[0],i[1]])
                contigs_dict[item] = contig_function_list

            print(self.count_synteny(contigs_dict))


            counter = counter + 1
            if counter == 1:
                break
            # for item in genes_and_functions_list:
                # contig_name = genes_and_functions_list[item][2]
                # function_and_genecallerid = (genes_and_functions_list[item][0],genes_and_functions_list[item][1])

                # if contig_name 
                # print (contig_name)

        # print(genes_and_functions_list)

            # do something with genes_and_functions_list
            # counts = self.count_synteny(genes_and_functions_list)

    def get_genes_and_functions_from_contigs_db(self, contigs_db_path):
        contigs_db = dbops.ContigsDatabase(contigs_db_path)

        genes_in_contigs = contigs_db.db.get_table_as_dict(t.genes_in_contigs_table_name)
        annotations_dict = contigs_db.db.get_table_as_dict(t.gene_function_calls_table_name)
        annotations_dict = utils.get_filtered_dict(annotations_dict, 'source', set([self.annotation_source]))

        # FIXME: turn genes_and_functions down below into a dictionary where keys are contig names and values are list of gene caller ids and function names.

        genes_and_functions = [(entry['gene_callers_id'], 
                                entry['function'], 
                                genes_in_contigs[entry['gene_callers_id']]['contig']) for entry in annotations_dict.values()]

        return genes_and_functions


    def count_synteny(self, contigs_dict):
        """
        Need to make an example where there is more than one contig in a contigsDB
        """
        for key in sorted(contigs_dict.keys()):
            gci_function = contigs_dict[key]
            genes = [(entry[1]) for entry in gci_function]

            kFreq = {}
            k = 3
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


    def driveSynteny(self):
        for key,value in self.genes:
            self.genes[key] = self.count_synteny(value)

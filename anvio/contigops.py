# -*- coding: utf-8
# pylint: disable=line-too-long
"""Classes and functions for handling, storing, and retrieving atomic data
   from contigs and splits. Also includes classes to deal with external
   contig data such as GenbankToAnvio."""

import os
import re
import io
import gzip
import string
import argparse

from Bio import SeqIO

import anvio
import anvio.tables as t
import anvio.utils as utils
import anvio.terminal as terminal
import anvio.constants as constants
import anvio.filesnpaths as filesnpaths

from anvio.sequence import Coverage
from anvio.errors import ConfigError
from anvio.variability import ColumnProfile
from anvio.variability import VariablityTestFactory


__author__ = "Developers of anvi'o (see AUTHORS.txt)"
__copyright__ = "Copyleft 2015-2018, the Meren Lab (http://merenlab.org/)"
__credits__ = ["Mike Lee", "Faruk Uzun"]
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "A. Murat Eren"
__email__ = "a.murat.eren@gmail.com"
__status__ = "Development"

OK_CHARS_FOR_ORGANISM_NAME = string.ascii_letters + string.digits + '_'
OK_CHARS_FOR_ACCESSION = OK_CHARS_FOR_ORGANISM_NAME + '.'


variability_test_class_default = VariablityTestFactory(params={'b': 2, 'm': 1.45, 'c': 0.05})
variability_test_class_null = VariablityTestFactory(params=None) # get everything for every coverage level


def gen_split_name(parent_name, order):
    return '_'.join([parent_name, 'split', '%05d' % (order + 1)])


def get_atomic_data_dicts(sample_id, contigs):
    """Takes a list of contigops.Contig objects, and returns contigs and splits atomic data
       dictionaries"""
    atomic_data_contigs = {}
    atomic_data_splits = {}

    # this loop will get atomic_data information from Contig instanes and store them into the db
    # at once. this was broken down into about 10 functions, but this structure seems to be the most efficient
    # although it looks crappy:
    for contig in contigs:
        contig_atomic_data = contig.get_atomic_data_dict()

        for split in contig.splits:
            atomic_data_contigs[split.name] = {'contig': contig.name}
            for atomic_data_field in t.atomic_data_table_structure[1:]:
                atomic_data_contigs[split.name][atomic_data_field] = contig_atomic_data[atomic_data_field]

        # contig is done, deal with splits in it:
        for split in contig.splits:
            split_atomic_data = split.get_atomic_data_dict()
            atomic_data_splits[split.name] = {'contig': split.name}
            for atomic_data_field in t.atomic_data_table_structure[1:]:
                atomic_data_splits[split.name][atomic_data_field] = split_atomic_data[atomic_data_field]

    return atomic_data_splits, atomic_data_contigs


class Contig:
    def __init__(self, name):
        self.name = name
        self.sequence = None
        self.parent = None
        self.splits = []
        self.length = 0
        self.abundance = 0.0
        self.coverage = Coverage()

        self.min_coverage_for_variability = 10
        self.skip_SNV_profiling = False
        self.report_variability_full = False
        self.ignore_orphans = True
        self.max_coverage_depth = constants.max_depth_for_coverage
        self.codon_frequencies_dict = {}


    def get_atomic_data_dict(self):
        d = {'std_coverage': self.coverage.std,
             'mean_coverage': self.coverage.mean,
             'mean_coverage_Q2Q3': self.coverage.mean_Q2Q3,
             'max_normalized_ratio': 1.0,
             'relative_abundance': 1.0,
             'detection': self.coverage.detection,
             'abundance': self.abundance,
             'variability': sum(s.auxiliary.variation_density for s in self.splits) if not self.skip_SNV_profiling else None,
             '__parent__': None}

        return d


    def analyze_coverage(self, bam):
        contig_coverage = []

        counter = 1
        for split in self.splits:
            split.coverage = Coverage()
            split.coverage.run(bam, split, 
                            ignore_orphans=self.ignore_orphans, 
                            max_coverage_depth=self.max_coverage_depth)
            contig_coverage.extend(split.coverage.c)

            counter += 1

        self.coverage.process_c(contig_coverage)


    def analyze_auxiliary(self, bam):
        counter = 1
        for split in self.splits:
            split.auxiliary = Auxiliary(split,
                                        bam,
                                        parent_outlier_positions=self.coverage.outlier_positions,
                                        min_coverage=self.min_coverage_for_variability,
                                        report_variability_full=self.report_variability_full,
                                        ignore_orphans=self.ignore_orphans,
                                        max_coverage_depth=self.max_coverage_depth)

            counter += 1


class Split:
    def __init__(self, name, sequence, parent, order, start=0, end=0):
        self.name = name
        self.sequence = sequence
        self.parent = parent
        self.end = end
        self.order = order
        self.start = start
        self.length = end - start
        self.explicit_length = 0
        self.abundance = 0.0
        self.column_profiles = {}
        self.auxiliary = None

    def get_atomic_data_dict(self):
        d = {'std_coverage': self.coverage.std,
             'mean_coverage': self.coverage.mean,
             'mean_coverage_Q2Q3': self.coverage.mean_Q2Q3,
             'max_normalized_ratio': 1.0,
             'relative_abundance': 1.0,
             'detection': self.coverage.detection,
             'abundance': self.abundance,
             'variability': self.auxiliary.variation_density if self.auxiliary else None,
             '__parent__': self.parent}

        return d


class Auxiliary:
    def __init__(self, split, bam, parent_outlier_positions, 
                 min_coverage=10, 
                 report_variability_full=False, 
                 ignore_orphans=True,
                 max_coverage_depth=constants.max_depth_for_coverage):
        self.v = []
        self.rep_seq = ''
        self.split = split
        self.variation_density = 0.0
        self.parent_outlier_positions = parent_outlier_positions
        self.competing_nucleotides = {}
        self.min_coverage = min_coverage
        self.column_profile = self.split.column_profiles
        self.report_variability_full = report_variability_full
        self.ignore_orphans = ignore_orphans
        self.max_coverage_depth = max_coverage_depth

        self.run(bam)


    def run(self, bam):
        ratios = []

        for pileupcolumn in bam.pileup(self.split.parent, self.split.start, self.split.end, 
                                    ignore_orphans=self.ignore_orphans, max_depth=self.max_coverage_depth):

            pos_in_contig = pileupcolumn.pos
            if pos_in_contig < self.split.start or pos_in_contig >= self.split.end:
                continue

            valid_nts = [pileupread.alignment.seq[pileupread.query_position] for pileupread in pileupcolumn.pileups if not pileupread.is_del and not pileupread.is_refskip]

            coverage = len(valid_nts)
            if coverage < self.min_coverage:
                continue

            column = ''.join(valid_nts)

            pos_in_split = pos_in_contig - self.split.start
            base_in_contig = self.split.sequence[pos_in_split]

            cp = ColumnProfile(column,
                               reference=base_in_contig,
                               coverage=coverage,
                               split_name=self.split.name,
                               pos=pos_in_split,
                               test_class=variability_test_class_null if self.report_variability_full else variability_test_class_default).profile

            if cp['worth_reporting']:
                ratios.append((cp['departure_from_reference'], cp['coverage']), )
                cp['pos_in_contig'] = pos_in_contig
                cp['cov_outlier_in_split'] = pos_in_split in self.split.coverage.outlier_positions
                cp['cov_outlier_in_contig'] = pos_in_contig in self.parent_outlier_positions
                self.column_profile[pos_in_contig] = cp

        # variation density = gene_callers_idber of SNVs per kb
        self.variation_density = len(ratios) * 1000.0 / self.split.length

        for i in range(self.split.start, self.split.end):
            if i in self.column_profile:
                self.rep_seq += self.column_profile[i]['reference']
                self.v.append(self.column_profile[i]['departure_from_reference'])
                self.competing_nucleotides[self.column_profile[i]['pos']] = self.column_profile[i]['competing_nts']
            else:
                self.rep_seq += 'N'
                self.v.append(0)


class GenbankToAnvioWrapper:
    """An ad hoc class that takes a metadata file generated by `ncbi-genome-download`"""

    def __init__(self, args, run=terminal.Run(), progress=terminal.Progress()):
        self.args = args

        self.run = run
        self.progress = progress

        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        self.metadata_file_path = A('metadata')
        self.output_directory_path = os.path.abspath(A('output_dir') or os.path.curdir)
        self.output_fasta_descriptor = A('output_fasta_txt') or os.path.join(self.output_directory_path, 'fasta-input.txt')
        self.exclude_gene_calls_from_fasta_txt = A('exclude_gene_calls_from_fasta_txt')


    def sanity_check(self):
        filesnpaths.is_file_tab_delimited(self.metadata_file_path)

        if os.path.exists(self.output_directory_path):
            filesnpaths.is_output_dir_writable(self.output_directory_path)
        else:
            filesnpaths.gen_output_directory(self.output_directory_path)

        filesnpaths.is_output_file_writable(self.output_fasta_descriptor)


    def process(self):
        self.sanity_check()

        self.run.info('Input metadata file', self.metadata_file_path)
        self.run.info('Output directory', self.output_directory_path)

        columns = utils.get_columns_of_TAB_delim_file(self.metadata_file_path)
        if 'organism_name' not in columns or 'local_filename' not in columns:
            raise ConfigError("The metadata file you provided does not look like a metadata\
                               file output from the program `ncbi-genome-download` :/ Why?\
                               Because anvi'o expects that file to have at least the following\
                               two columns in it: 'organism_name' and 'local_filename'.")

        metadata = utils.get_TAB_delimited_file_as_dictionary(self.metadata_file_path)

        for entry in metadata:
            if not os.path.exists(metadata[entry]['local_filename']):
                raise ConfigError("At least one of the files in your metadata input does not seem to be\
                                   where they think they are :/ Please make sure the entry %s and others\
                                   point to proper local file paths..." % entry)

        self.run.info('Num entries in metadata', len(metadata))

        output_fasta_dict = {}
        self.progress.new("GenBank to anvi'o", progress_total_items=len(metadata))
        for entry in metadata:
            self.progress.increment()
            self.progress.update('Processing %s ...' % entry)

            # set the organism name and accession id and clean them from weird
            # characters.
            organism_name = metadata[entry]['organism_name']
            for char in [c for c in organism_name if c not in OK_CHARS_FOR_ORGANISM_NAME]:
                organism_name = organism_name.replace(char, '_')

            accession_id = entry
            for char in [c for c in accession_id if c not in OK_CHARS_FOR_ACCESSION]:
                accession_id = accession_id.replace(char, '_')

            final_name = '_'.join([organism_name, accession_id])

            args = argparse.Namespace(input_genbank=metadata[entry]['local_filename'],
                                      output_file_prefix=os.path.join(self.output_directory_path, final_name))
            g = GenbankToAnvio(args, run=terminal.Run(verbose=False), progress=terminal.Progress(verbose=False))

            if final_name in output_fasta_dict:
                raise ConfigError("The final name '%s' for your genome has alrady been used by\
                                   another one :/ This should never happen unless your metadata\
                                   contains entries with identical accession numbers...")
            output_fasta_dict[final_name] = g.process()

        self.progress.end()

        headers = ['name', 'path', 'gene_functional_annotation']
        if not self.exclude_gene_calls_from_fasta_txt:
            headers.append('external_gene_calls')

        utils.store_dict_as_TAB_delimited_file(output_fasta_dict, self.output_fasta_descriptor, headers=headers)

        self.run.info('Output FASTA descriptor', self.output_fasta_descriptor)


class GenbankToAnvio:
    """A class to deal with GenBank files. The initial code was implemented by Mike Lee"""

    def __init__(self, args, run=terminal.Run(), progress=terminal.Progress()):
        self.args = args

        self.run = run
        self.progress = progress

        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        self.input_genbank_path = A('input_genbank')
        self.output_file_prefix = A('output_file_prefix')
        self.output_fasta_path = A('output_fasta')
        self.output_functions_path = A('output_functions')
        self.output_gene_calls_path = A('output_gene_calls')
        self.source = A('annotations_source') or 'NCBI_PGAP'
        self.version = A('annotations_source') or 'v4.6'

        # gene callers id start from 0. you can change your instance
        # prior to processing the genbank file to start from another
        # value
        self.gene_callers_id = 0

        # dumping gene if noted as these in the "note" section of the call
        self.note_terms_to_exclude = ["frameshifted", "internal stop", "incomplete"]

        # dumping gene if "location" section contains any of these: "join" means the
        # gene call spans multiple contigs; "<" or ">" means the gene call runs off a contig
        self.location_terms_to_exclude = ["join", "<", ">"]


    def sanity_check(self):
        if self.output_file_prefix and (self.output_fasta_path or self.output_functions_path or self.output_gene_calls_path):
            raise ConfigError("Your arguments contain an output file prefix, and other output file paths. You can either\
                               define a prefix, and the output files would be named accordingly (such as 'PREFIX-extenral-gene-calls',\
                               'PREFIX-external-functions.txt', and 'PREFIX-contigs.fa'), ORRR you can set output file names\
                               or paths for each of these files independently. You can also leave it as is for default file names to\
                               be used. But you can't mix everything together and confuse us here.")

        self.output_fasta_path = self.output_fasta_path or 'contigs.fa'
        self.output_functions_path = self.output_functions_path or 'external-functions.txt'
        self.output_gene_calls_path = self.output_gene_calls_path or 'external-gene-calls.txt'

        if self.output_file_prefix:
            J = lambda x: '-'.join([self.output_file_prefix, x])
            self.output_fasta_path = J(self.output_fasta_path)
            self.output_functions_path = J(self.output_functions_path)
            self.output_gene_calls_path = J(self.output_gene_calls_path)

        filesnpaths.is_output_file_writable(self.output_fasta_path)
        filesnpaths.is_output_file_writable(self.output_functions_path)
        filesnpaths.is_output_file_writable(self.output_gene_calls_path)
        filesnpaths.is_file_exists(self.input_genbank_path)

        files_already_exist = [f for f in [self.output_fasta_path, self.output_functions_path, self.output_gene_calls_path] if os.path.exists(f)]
        if len(files_already_exist):
            raise ConfigError("Some of the output files already exist :/ Anvi'o feels uneasy about simply overwriting\
                               them and would like to outsource that risk to you. Please either use different output\
                               file names, or delete these files and come back: '%s'" % (', '.join(files_already_exist)))


    def process(self):
        self.sanity_check()

        output_fasta = {}
        output_gene_calls = {}
        output_functions = {}
        num_genbank_records_processed = 0
        num_genes_found = 0
        num_genes_reported = 0
        num_genes_with_functions = 0

        try:
            if self.input_genbank_path.endswith('.gz'):
                genbank_file_object = SeqIO.parse(io.TextIOWrapper(gzip.open(self.input_genbank_path, 'r')), "genbank")
            else:
                genbank_file_object = SeqIO.parse(open(self.input_genbank_path, "r"), "genbank")
        except Exception as e:
            raise ConfigError("Someone didn't like your unput 'genbank' file :/ Here's what they said\
                               about it: '%s'." % e)

        for genbank_record in genbank_file_object:
            num_genbank_records_processed += 1

            output_fasta[genbank_record.name] = str(genbank_record.seq)

            genes = [gene for gene in genbank_record.features if gene.type =="CDS"] # focusing on features annotated as "CDS" by NCBI's PGAP

            for gene in genes:
                num_genes_found += 1
                location = str(gene.location)

                # dumping gene if "location" section contains any of these terms set above: "join" means the gene call spans multiple contigs; "<" or ">" means the gene call runs off a contig
                if any(exclusion_term in location for exclusion_term in self.location_terms_to_exclude):
                    continue

                if "note" in gene.qualifiers:
                    note = str(gene.qualifiers["note"][0])

                    # dumping gene if noted as any of these in the "note" section set above
                    if any(exclusion_term in note for exclusion_term in self.note_terms_to_exclude):
                        continue

                # dumping if overlapping translation frame
                if "transl_except" in gene.qualifiers:
                    continue

                # dumping if gene declared a pseudogene
                if "pseudo" in gene.qualifiers:
                    continue

                # cleaning up gene coordinates to more easily parse:
                location = location.replace("[", "")
                location = re.sub('](.*)', '', location)
                location = location.split(":")

                start = location[0] # start coordinate
                end = location[1] # end coordinate


                # setting direction to "f" or "r":
                if gene.strand == 1:
                    direction="f"
                else:
                    direction="r"

                # for accession, storing protein id if it has one, else the the locus tag, else "None"
                if "protein_id" in gene.qualifiers:
                    accession = gene.qualifiers["protein_id"][0]
                elif "locus_tag" in gene.qualifiers:
                    accession = gene.qualifiers["locus_tag"][0]
                else:
                    accession = "None"

                # storing gene product annotation
                function = gene.qualifiers["product"][0]

                # if present, adding gene name to product annotation:
                if "gene" in gene.qualifiers:
                    gene_name=str(gene.qualifiers["gene"][0])
                    function = function + " (" + gene_name + ")"

                output_gene_calls[self.gene_callers_id] = {'contig': genbank_record.name,
                                                           'start': start,
                                                           'stop': end,
                                                           'direction': direction,
                                                           'partial': 0,
                                                           'source': self.source,
                                                           'version': self.version}
                num_genes_reported += 1

                # not writing gene out to functions table if no annotation
                if "hypothetical protein" not in function:
                    output_functions[self.gene_callers_id] = {'source': self.source,
                                                              'accession': accession,
                                                              'function': function,
                                                              'e_value': 0}
                    num_genes_with_functions += 1

                # increment the gene callers id fo rthe next
                self.gene_callers_id += 1

        if num_genbank_records_processed == 0:
            raise ConfigError("It seems there was no records in your input genbank file :/ Are you sure you\
                               gave the right file path that actually resolves to a genbank formatted\
                               text file?")

        self.run.info('Num GenBank entries processed', num_genbank_records_processed)
        self.run.info('Num gene records found', num_genes_found)
        self.run.info('Num genes reported', num_genes_reported, mc='green')
        self.run.info('Num genes with functins', num_genes_with_functions, mc='green', nl_after=1)

        # time to write these down:
        utils.store_dict_as_FASTA_file(output_fasta,
                                       self.output_fasta_path,
                                       wrap_from=None)
        self.run.info('FASTA file path', self.output_fasta_path)

        if len(output_gene_calls):
            utils.store_dict_as_TAB_delimited_file(output_gene_calls,
                                                   self.output_gene_calls_path,
                                                   headers=["gene_callers_id", "contig", "start", "stop", "direction", "partial", "source", "version"])
            self.run.info('External gene calls file', self.output_gene_calls_path)

            utils.store_dict_as_TAB_delimited_file(output_functions,
                                                   self.output_functions_path,
                                                   headers=['gene_callers_id', 'source', 'accession', 'function', 'e_value'])
            self.run.info('TAB-delimited functions', self.output_functions_path)
        else:
            self.output_gene_calls_path = None
            self.output_functions_path = None
            self.run.warning("Anvi'o couldn't find any gene calles in the GenBank file, hence you will get\
                              no output files for external gene calls or functions :/ We hope you can\
                              survive this terrible terrible news :(")

        self.run.info_single('Mmmmm ☘ ', nl_before=1, nl_after=1)

        return {'external_gene_calls': self.output_gene_calls_path,
                'gene_functional_annotation': self.output_functions_path,
                'path': self.output_fasta_path}

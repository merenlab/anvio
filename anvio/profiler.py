# coding: utf-8
# pylint: disable=line-too-long
"""Provides the necessary class to profile BAM files."""

import gc
import os
import sys
import copy
import shutil
import argparse
import numpy as np

# multiprocess is a fork of multiprocessing that uses the dill serializer instead of pickle
# using the multiprocessing module directly results in a pickling error in Python 3.10 which
# goes like this:
#
#   >>> AttributeError: Can't pickle local object 'SOMEFUNCTION.<locals>.<lambda>' multiprocessing
#
import multiprocess as multiprocessing

from collections import OrderedDict

import anvio
import anvio.tables as t
import anvio.dbops as dbops
import anvio.utils as utils
import anvio.bamops as bamops
import anvio.terminal as terminal
import anvio.contigops as contigops
import anvio.constants as constants
import anvio.clustering as clustering
import anvio.filesnpaths as filesnpaths
import anvio.auxiliarydataops as auxiliarydataops
import anvio.streamingops as streamingops

from anvio.errors import ConfigError
from anvio.tables.views import TablesForViews
from anvio.tables.indels import TableForIndels
from anvio.tables.miscdata import TableForLayerAdditionalData
from anvio.tables.variability import TableForVariability
from anvio.tables.codonfrequencies import TableForCodonFrequencies


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "A. Murat Eren"
__email__ = "a.murat.eren@gmail.com"

null_progress = terminal.Progress(verbose=False)
null_run = terminal.Run(verbose=False)
pp = terminal.pretty_print


class BAMProfilerQuick:
    """A class for rapid profiling of BAM files that produce text files rather than profile dbs"""

    def __init__(self, args, skip_sanity_check=False, run=terminal.Run(), progress=terminal.Progress()):
        self.args = args
        self.run = run
        self.progress = progress

        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        self.bam_file_paths = A('bam_files')
        self.contigs_db_path = A('contigs_db')
        self.output_file_path = A('output_file')
        self.gene_level_stats = A('gene_mode')
        self.gene_caller = A('gene_caller')
        self.gene_caller_ids = A('gene_caller_ids')
        self.genes_of_interest_file = A('genes_of_interest')
        self.report_minimal = A('report_minimal')
        self.collection_txt_path = A('collection_txt')

        if not self.gene_caller:
            self.gene_caller = utils.get_default_gene_caller(self.contigs_db_path)

        if not skip_sanity_check:
            self.sanity_check()

        self.run.info('Contigs DB', self.contigs_db_path)
        self.run.info('Num BAM files', len(self.bam_file_paths))
        self.run.info('Reporting', 'MINIMAL' if self.report_minimal else 'EVERYTHING', mc="red" if self.report_minimal else "green")

        # if requested, load genes of interest
        self.gene_ids_of_interest = set([])
        id_list = []
        if self.genes_of_interest_file:
            id_list = [g.strip() for g in open(self.genes_of_interest_file, 'r').readlines()]
        elif self.gene_caller_ids:
            id_list = self.gene_caller_ids.split(',')
        nonnumeric = [i for i in id_list if not i.isnumeric()]
        self.gene_ids_of_interest = set(int(i) for i in id_list if i not in nonnumeric)
        if nonnumeric:
            raise ConfigError("Some of the gene caller IDs you requested do not look like gene caller IDs. Here they are "
                             f"so you can remove them from your request: {', '.join(nonnumeric)}")
        
        # to be filled later if necessary
        self.contigs_basic_info = {}
        self.gene_calls_per_contig = {}
        self.collection_bins = OrderedDict()
        self.contig_names_to_process = []


    def sanity_check(self):
        if not self.contigs_db_path:
            raise ConfigError("You need to provide an anvi'o contigs database for this to work :/")

        utils.is_contigs_db(self.contigs_db_path)

        if not len(self.bam_file_paths):
            raise ConfigError("You need to provide at least one BAM file for this to work.")

        if not self.output_file_path:
            raise ConfigError("Please provide an output file path.")

        filesnpaths.is_output_file_writable(self.output_file_path, ok_if_exists=False)
        if self.collection_txt_path:
            if self.gene_level_stats:
                raise ConfigError("The flag `--collection-txt` is only available when reporting contig-level "
                                  "stats. Please drop --gene-mode or skip the collection.")
            filesnpaths.is_file_plain_text(self.collection_txt_path)

        if not self.gene_level_stats and (self.gene_caller_ids or self.genes_of_interest_file):
            raise ConfigError("You requested genes of interest but not --gene-mode, which doesn't make sense. So "
                              "we are stopping you right here, just to make sure you know what you be doin'")
        if self.gene_caller_ids and self.genes_of_interest_file:
            raise ConfigError("The parameter --gene-caller-ids is not compatible with --genes-of-interest. "
                              "Please provide only one.")
        elif self.genes_of_interest_file:
            filesnpaths.is_file_exists(self.genes_of_interest_file)

        # find all the bad BAM files
        self.progress.new("Sanity checking BAM files")
        self.progress.update('...')
        bad_bam_files = [f for f in self.bam_file_paths if not filesnpaths.is_file_bam_file(f, dont_raise=True)]
        if len(bad_bam_files):
            self.progress.reset()
            raise ConfigError(f"Not all of your BAM files look like BAM files. Here is the list that "
                              f"samtools didn't like: {', '.join(bad_bam_files)}")
        self.progress.end()

        if self.gene_level_stats:
            contigs_db = dbops.ContigsDatabase(self.contigs_db_path)
            genes_are_called = contigs_db.meta['genes_are_called']
            gene_callers_list = contigs_db.meta['gene_callers']
            contigs_db.disconnect()

            if not genes_are_called:
                raise ConfigError("There are no gene calls in this contigs database :/ You can't use the flag "
                                  "`--gene-mode`. Yes.")

            gene_callers = [tpl[0] for tpl in gene_callers_list]
            if self.gene_caller not in gene_callers:
                dbops.ContigsDatabase(self.contigs_db_path).list_gene_caller_sources()
                raise ConfigError(f"The gene caller '{self.gene_caller}' is not among those that are found in "
                                  f"the contigs database (and shown above for your convenience).")



    def process(self):
        contigs_db = dbops.ContigsDatabase(self.contigs_db_path)

        self.run.warning(None, header="CONTIGS DB", lc="green")
        self.run.info("Number of contigs", pp(contigs_db.meta['num_contigs']))
        self.run.info("Total num nucleotides", pp(contigs_db.meta['total_length']))
        self.run.info("Genes are called", "Yes" if contigs_db.meta['genes_are_called'] else "No :/")

        self.run.info("Reporting mode", "GENES" if self.gene_level_stats else "CONTIGS", mc='green')
        if self.gene_level_stats:
            self.run.info("Gene caller", self.gene_caller)
            self.run.info("Number of genes", pp([tpl[1] for tpl in contigs_db.meta['gene_callers'] if tpl[0] == self.gene_caller][0]))
            self.run.info("Gene caller IDs of interest", "ALL GENES" if not self.gene_ids_of_interest else ",".join([str(i) for i in self.gene_ids_of_interest]))


        self.progress.new('Reading data into memory')
        self.progress.update('Contigs basic info table ...')
        self.contigs_basic_info = contigs_db.db.get_table_as_dict(t.contigs_info_table_name)
        self.progress.end()

        contigs_db.disconnect()

        self.contig_names_to_process = list(self.contigs_basic_info.keys())

        if self.collection_txt_path:
            self.init_collection_bins()

        if self.gene_level_stats:
            self.recover_gene_data()

        self.report_stats()


    def recover_gene_data(self):
        # if we are working with genes, we need to first read the gene calls table into
        # memory
        self.progress.new('Reading data into memory')
        self.progress.update('Gene calls table ...')
        contigs_db = dbops.ContigsDatabase(self.contigs_db_path)
        genes_where_clause = f"source == '{self.gene_caller}'"
        if self.gene_ids_of_interest:
            genes_where_clause += f" and gene_callers_id in ({', '.join([str(g) for g in self.gene_ids_of_interest])})" 
        self.genes_in_contigs = contigs_db.db.get_some_rows_from_table_as_dict(t.genes_in_contigs_table_name, where_clause=genes_where_clause, error_if_no_data=True)

        # little sanity check to make sure we got all the gene calls the user asked for
        if self.gene_ids_of_interest:
            missing_ids = [str(g) for g in self.gene_ids_of_interest if g not in self.genes_in_contigs]
            if missing_ids:
                self.progress.reset()
                raise ConfigError(f"Some of the gene caller IDs you requested were not found in the contigs database. "
                                  f"As this could be Very Bad News depending on what you are doing, we'll stop the show "
                                  f"here and give you a chance to double check those genes (and remove them from your list "
                                  f"of genes of interest, if you find that they actually aren't the IDs you were looking for). "
                                  f"Here are the affected gene caller IDs: {', '.join(missing_ids)}")

        # and then update contigs basic info to easily track gene calls in a given contig
        # for reporting purposes
        contigs_with_genes = set()
        self.progress.update('Updating contigs basic info ...')
        for gene_callers_id in self.genes_in_contigs:
            contig_name = self.genes_in_contigs[gene_callers_id]['contig']
            contigs_with_genes.add(contig_name)

            if 'gene_caller_ids' not in self.contigs_basic_info[contig_name]:
                self.contigs_basic_info[contig_name]['gene_caller_ids'] = set([gene_callers_id])
            else:
                self.contigs_basic_info[contig_name]['gene_caller_ids'].add(gene_callers_id)

        self.progress.end()

        contigs_db.disconnect()

        self.contig_names_to_process = list(contigs_with_genes)
        self.run.info('Contigs with gene calls', len(self.contig_names_to_process))
        if len(self.contig_names_to_process) < len(self.contigs_basic_info):
            self.run.info('Contigs ignored (no genes)', len(self.contigs_basic_info) - len(self.contig_names_to_process))


    def init_collection_bins(self):
        """Load contig-to-bin mappings from a collection-txt and precompute bin metadata.

        The standard collection-txt has two columns (contig name, bin name). Here we validate
        that every contig listed exists in the contigs database and compute per-bin length and
        GC-content so downstream coverage summaries can operate at the bin/genome level.
        """

        collection_entries = utils.get_TAB_delimited_file_as_dictionary(self.collection_txt_path,
                                                                        no_header=True,
                                                                        column_names=['contig', 'bin'])

        self.collection_bins = OrderedDict()
        contig_to_bin = {}
        bins = OrderedDict()
        missing_contigs = set([])

        for contig_name, entry in collection_entries.items():
            bin_name = entry['bin']

            if contig_name in contig_to_bin:
                raise ConfigError(f"The contig '{contig_name}' appears more than once in the collection-txt. "
                                  f"Each contig can only belong to a single bin for `anvi-profile-blitz`.")
            contig_to_bin[contig_name] = bin_name

            if contig_name not in self.contigs_basic_info:
                missing_contigs.add(contig_name)
                continue

            if bin_name not in bins:
                bins[bin_name] = []
            bins[bin_name].append(contig_name)

        if missing_contigs:
            missing_contigs = sorted(list(missing_contigs))
            missing_contig_examples = "', '".join(missing_contigs[:5])
            if len(missing_contigs) > 5:
                missing_contig_examples += f"', and {len(missing_contigs) - 5} more"
            raise ConfigError(f"The collection-txt lists contigs that are not present in your contigs database :/ "
                              f"Here is an example: '{missing_contig_examples}'). Someone's gotta do something "
                              f"about it.")

        if not len(bins):
            raise ConfigError("None of the contigs listed in your collection-txt are present in the contigs database o_O "
                              "There is nothing to report (apart from this).")

        for bin_name, contigs in bins.items():
            bin_length = 0
            weighted_gc = 0

            for contig_name in contigs:
                contig_length = int(self.contigs_basic_info[contig_name]['length'])
                gc_content = float(self.contigs_basic_info[contig_name]['gc_content'])
                bin_length += contig_length
                weighted_gc += gc_content * contig_length

            gc_content = (weighted_gc / bin_length) if bin_length else 0

            self.collection_bins[bin_name] = {'contigs': contigs,
                                              'length': bin_length,
                                              'gc_content': gc_content}

        self.contig_names_to_process = [contig for data in self.collection_bins.values() for contig in data['contigs']]

        self.run.info('Collection', self.collection_txt_path)
        self.run.info('Bins in collection', len(self.collection_bins))
        self.run.info('Contigs in collection', len(self.contig_names_to_process))
        if len(self.contig_names_to_process) < len(self.contigs_basic_info):
            self.run.info('Contigs ignored (not in collection)', len(self.contigs_basic_info) - len(self.contig_names_to_process))


    def summarize_bin_coverages(self, coverage_arrays):
        """Concatenate contig coverages so bin metrics rely on nucleotide-level values.

        Bin-level detection and coverage can not be inferred from per-contig averages, so we
        stitch all contig coverage arrays together and let downstream statistics treat the bin
        as a single genome.
        """
        if not len(coverage_arrays):
            return None

        return np.concatenate(coverage_arrays)


    def _update_progress(self, mem_tracker, mem_usage, mem_diff, bam_index, contig_index, num_contigs, total_num_bam_files):
        if contig_index == 1 or contig_index % 100 == 0:
            if mem_tracker.measure():
                mem_usage = mem_tracker.get_last()
                mem_diff = mem_tracker.get_last_diff()

            self.progress.increment(increment_to=((bam_index * num_contigs) + contig_index))
            self.progress.update(f"BAM {bam_index + 1}/{pp(total_num_bam_files)} :: CONTIG {pp(contig_index)}/"
                                 f"{pp(num_contigs)} :: MEMORY 🧠 {mem_usage} ({mem_diff})")

        return mem_usage, mem_diff


    def report_stats(self):
        """Iterates through bam files, reports contigs stats"""

        contig_names = self.contig_names_to_process if len(self.contig_names_to_process) else list(self.contigs_basic_info.keys())
        num_contigs = len(contig_names)

        if not num_contigs:
            raise ConfigError("There are no contigs to profile. Please check your inputs.")

        # the total number of items we will process is equal to the number of contigs that will
        # have to be processed for each BAM file
        total_num_items = len(self.bam_file_paths) * num_contigs
        total_num_bam_files = len(self.bam_file_paths)
        reporting_bins = len(self.collection_bins) > 0
        reporting_genes = self.gene_level_stats and not reporting_bins

        self.progress.new("Bleep bloop", progress_total_items=total_num_items)
        self.progress.update('...')

        mem_tracker = terminal.TrackMemory(at_most_every=5)
        mem_usage, mem_diff = mem_tracker.start()

        with open(self.output_file_path, 'w') as output:
            header = self._get_header(reporting_bins, reporting_genes)

            output.write('\t'.join(header) + '\n')

            for i in range(0, total_num_bam_files):
                bam_file_path = self.bam_file_paths[i]

                bam = bamops.BAMFileObject(bam_file_path, 'rb')
                bam_file_name = os.path.splitext(os.path.basename(bam_file_path))[0]

                if reporting_bins:
                    mem_usage, mem_diff = self._process_bins_for_bam(bam, bam_file_name, output, mem_tracker, mem_usage, mem_diff, i, num_contigs, total_num_bam_files)
                elif reporting_genes:
                    mem_usage, mem_diff = self._process_genes_for_bam(bam, bam_file_name, output, contig_names, mem_tracker, mem_usage, mem_diff, i, num_contigs, total_num_bam_files)
                else:
                    mem_usage, mem_diff = self._process_contigs_for_bam(bam, bam_file_name, output, contig_names, mem_tracker, mem_usage, mem_diff, i, num_contigs, total_num_bam_files)

                bam.close()

            self.progress.end()


    def _get_header(self, reporting_bins, reporting_genes):
        if reporting_bins:
            if self.report_minimal:
                return ['bin', 'sample', 'length', 'gc_content', 'num_mapped_reads', 'detection', 'mean_cov']
            else:
                return ['bin', 'sample', 'length', 'gc_content', 'num_mapped_reads', 'detection', 'mean_cov', 'q2q3_cov', 'median_cov', 'min_cov', 'max_cov', 'std_cov']

        elif reporting_genes:
            if self.report_minimal:
                return ['gene_callers_id', 'contig', 'sample', 'length', 'detection', 'mean_cov']
            else:
                return ['gene_callers_id', 'contig', 'sample', 'length', 'num_mapped_reads', 'detection', 'mean_cov', 'q2q3_cov', 'median_cov', 'min_cov', 'max_cov', 'std_cov']

        else:
            if self.report_minimal:
                return ['contig', 'sample', 'length', 'gc_content', 'num_mapped_reads', 'detection', 'mean_cov']
            else:
                return ['contig', 'sample', 'length', 'gc_content', 'num_mapped_reads', 'detection', 'mean_cov', 'q2q3_cov', 'median_cov', 'min_cov', 'max_cov', 'std_cov']

        raise ConfigError("This function reached a point it should have never :(")


    def _process_bins_for_bam(self, bam, bam_file_name, output, mem_tracker, mem_usage, mem_diff, bam_index, num_contigs, total_num_bam_files):
        """Process bin-level reporting for a single BAM file."""

        contigs_processed = 0

        for bin_name, bin_data in self.collection_bins.items():
            bin_num_reads = 0
            total_positions = 0
            total_cov = 0
            total_detected_positions = 0

            combined_coverage = None
            cursor = 0

            for contig_name in bin_data['contigs']:
                contigs_processed += 1
                mem_usage, mem_diff = self._update_progress(mem_tracker, mem_usage, mem_diff, bam_index, contigs_processed, num_contigs, total_num_bam_files)

                coverage_obj = self._get_contig_coverage(bam, contig_name)
                if coverage_obj is None:
                    continue

                bin_num_reads += coverage_obj.num_reads

                # this part of the code could be written in a much more easy-to-read fashion. but it is the
                # way it is just to achieve maximum memory use efficiency .. especially in minimal mode.
                # when --report-minimal is used in collection mode, the coverages will  accumulate on the fly,
                # so per-bin nucleotide arrays will not be kept in memory. and for full reports, bin coverages
                # are written into a single preallocated numpy array per bin (without any concatenation as
                # data accummulates) which reduces peak RAM. while this code makes sure minimal mode is as
                # fast and as memory-efficient as possible, the full bin mode will still needs one combined
                # coverage array per bin to compute detection and quantiles, but at least it avoids
                # double-buffering overhead during copy :/
                if self.report_minimal:
                    total_positions += len(coverage_obj.c)
                    total_cov += np.sum(coverage_obj.c)
                    total_detected_positions += np.count_nonzero(coverage_obj.c)
                else:
                    if combined_coverage is None:
                        combined_coverage = np.empty(bin_data['length'], dtype=coverage_obj.c.dtype)

                    combined_coverage[cursor:cursor + len(coverage_obj.c)] = coverage_obj.c
                    cursor += len(coverage_obj.c)

            if self.report_minimal:
                if not total_positions:
                    continue

                mean = total_cov / total_positions
                detection = total_detected_positions / total_positions
                self._write_bin_stats_minimal(output, bin_name, bam_file_name, bin_data, mean, detection, bin_num_reads)
                continue

            if combined_coverage is None or cursor == 0:
                continue

            combined_coverage = combined_coverage[:cursor]

            self._write_bin_stats(output, bin_name, bam_file_name, bin_data, combined_coverage, bin_num_reads)

        return mem_usage, mem_diff


    def _process_contigs_for_bam(self, bam, bam_file_name, output, contig_names, mem_tracker, mem_usage, mem_diff, bam_index, num_contigs, total_num_bam_files):
        """Process contig-level reporting for a single BAM file."""

        contigs_processed = 0

        for contig_name in contig_names:
            contigs_processed += 1
            mem_usage, mem_diff = self._update_progress(mem_tracker, mem_usage, mem_diff, bam_index, contigs_processed, num_contigs, total_num_bam_files)

            coverage_obj = self._get_contig_coverage(bam, contig_name)
            if coverage_obj is None:
                continue

            self._write_contig_stats(output, contig_name, bam_file_name, coverage_obj)

        return mem_usage, mem_diff


    def _process_genes_for_bam(self, bam, bam_file_name, output, contig_names, mem_tracker, mem_usage, mem_diff, bam_index, num_contigs, total_num_bam_files):
        """Process gene-level reporting for a single BAM file."""

        contigs_processed = 0

        for contig_name in contig_names:
            contigs_processed += 1
            mem_usage, mem_diff = self._update_progress(mem_tracker, mem_usage, mem_diff, bam_index, contigs_processed, num_contigs, total_num_bam_files)

            if 'gene_caller_ids' not in self.contigs_basic_info[contig_name]:
                continue

            coverage_obj = self._get_contig_coverage(bam, contig_name)
            if coverage_obj is None:
                continue

            for gene_callers_id in self.contigs_basic_info[contig_name]['gene_caller_ids']:
                g = self.genes_in_contigs[gene_callers_id]
                gc = coverage_obj.c[g['start']:g['stop']]

                if self.report_minimal:
                    num_mapped_reads = None
                else:
                    num_mapped_reads = 0
                    for r in bam.fetch_only(contig_name, start=g['start'], end=g['stop']):
                        num_mapped_reads += 1

                self._write_gene_stats(output, gene_callers_id, contig_name, bam_file_name, g['stop'] - g['start'], num_mapped_reads, gc)

        return mem_usage, mem_diff


    def _get_contig_coverage(self, bam, contig_name):
        coverage_obj = bamops.Coverage()

        try:
            coverage_obj.run(bam, contig_name, read_iterator='fetch', skip_coverage_stats=True)
        except:
            return None

        return coverage_obj


    def _write_bin_stats(self, output, bin_name, bam_file_name, bin_data, coverage_array, num_reads):
        if self.report_minimal:
            mean = np.mean(coverage_array)
            detection = np.sum(coverage_array > 0) / len(coverage_array)
            self._write_bin_stats_minimal(output, bin_name, bam_file_name, bin_data, mean, detection, num_reads)
        else:
            C = utils.CoverageStats(coverage_array, skip_outliers=True)
            output.write(f"{bin_name}\t"
                         f"{bam_file_name}\t"
                         f"{bin_data['length']}\t"
                         f"{float(bin_data['gc_content']):.3}\t"
                         f"{num_reads}\t"
                         f"{C.detection:.4}\t"
                         f"{C.mean:.4}\t"
                         f"{C.mean_Q2Q3:.4}\t"
                         f"{C.median}\t"
                         f"{C.min}\t"
                         f"{C.max}\t"
                         f"{C.std:.4}\n")


    def _write_bin_stats_minimal(self, output, bin_name, bam_file_name, bin_data, mean, detection, num_reads):
        output.write(f"{bin_name}\t"
                     f"{bam_file_name}\t"
                     f"{bin_data['length']}\t"
                     f"{float(bin_data['gc_content']):.3}\t"
                     f"{num_reads}\t"
                     f"{detection:.4}\t"
                     f"{mean:.4}\n")


    def _write_contig_stats(self, output, contig_name, bam_file_name, coverage_obj):
        if self.report_minimal:
            mean = np.mean(coverage_obj.c)
            detection = np.sum(coverage_obj.c > 0) / len(coverage_obj.c)
            output.write(f"{contig_name}\t"
                         f"{bam_file_name}\t"
                         f"{self.contigs_basic_info[contig_name]['length']}\t"
                         f"{float(self.contigs_basic_info[contig_name]['gc_content']):.3}\t"
                         f"{coverage_obj.num_reads}\t"
                         f"{detection:.4}\t"
                         f"{mean:.4}\n")
        else:
            C = utils.CoverageStats(coverage_obj.c, skip_outliers=True)
            output.write(f"{contig_name}\t"
                         f"{bam_file_name}\t"
                         f"{self.contigs_basic_info[contig_name]['length']}\t"
                         f"{float(self.contigs_basic_info[contig_name]['gc_content']):.3}\t"
                         f"{coverage_obj.num_reads}\t"
                         f"{C.detection:.4}\t"
                         f"{C.mean:.4}\t"
                         f"{C.mean_Q2Q3:.4}\t"
                         f"{C.median}\t"
                         f"{C.min}\t"
                         f"{C.max}\t"
                         f"{C.std:.4}\n")


    def _write_gene_stats(self, output, gene_callers_id, contig_name, bam_file_name, gene_length, num_mapped_reads, coverage_array):
        if self.report_minimal:
            mean = np.mean(coverage_array)
            detection = np.sum(coverage_array > 0) / len(coverage_array)
            output.write(f"{gene_callers_id}\t"
                         f"{contig_name}\t"
                         f"{bam_file_name}\t"
                         f"{gene_length}\t"
                         f"{detection:.4}\t"
                         f"{mean:.4}\n")
        else:
            GC = utils.CoverageStats(coverage_array, skip_outliers=True)
            output.write(f"{gene_callers_id}\t"
                         f"{contig_name}\t"
                         f"{bam_file_name}\t"
                         f"{gene_length}\t"
                         f"{num_mapped_reads}\t"
                         f"{GC.detection:.4}\t"
                         f"{GC.mean:.4}\t"
                         f"{GC.mean_Q2Q3:.4}\t"
                         f"{GC.median}\t"
                         f"{GC.min}\t"
                         f"{GC.max}\t"
                         f"{GC.std:.4}\n")


class BAMProfiler(dbops.ContigsSuperclass):
    """Creates an über class for BAM file operations"""

    def __init__(self, args, r=terminal.Run(width=50), p=terminal.Progress()):
        self.args = args
        self.progress = p
        self.run = r

        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        self.input_file_path = A('input_file')
        self.contigs_db_path = A('contigs_db')
        self.serialized_profile_path = A('serialized_profile')
        self.output_directory = A('output_dir')
        self.list_contigs_and_exit = A('list_contigs')
        self.min_contig_length = A('min_contig_length') or 0
        self.max_contig_length = A('max_contig_length') or sys.maxsize
        self.min_coverage_for_variability = A('min_coverage_for_variability')
        self.contigs_shall_be_clustered = A('cluster_contigs')
        self.skip_hierarchical_clustering = A('skip_hierarchical_clustering')
        self.sample_id = A('sample_name')
        self.report_variability_full = A('report_variability_full')
        self.skip_edges = A('skip_edges')
        self.overwrite_output_destinations = A('overwrite_output_destinations')
        self.skip_SNV_profiling = A('skip_SNV_profiling')
        self.skip_INDEL_profiling = A('skip_INDEL_profiling')
        self.profile_SCVs = A('profile_SCVs')
        self.processing_chunk_size = A('processing_chunk_size')  # None means use split length
        self.min_percent_identity = A('min_percent_identity')
        self.fetch_filter = A('fetch_filter')
        self.gen_serialized_profile = A('gen_serialized_profile')
        self.distance = A('distance') or constants.distance_metric_default
        self.linkage = A('linkage') or constants.linkage_method_default
        self.num_threads = int(A('num_threads') or 1)
        self.queue_size = int(A('queue_size') if A('queue_size') is not None else 0)
        self.write_buffer_size_per_thread = int(A('write_buffer_size_per_thread') if A('write_buffer_size_per_thread') is not None else 500)
        self.write_buffer_size = self.write_buffer_size_per_thread * self.num_threads
        self.total_length_of_all_contigs = 0
        self.total_coverage_values_for_all_contigs = 0
        self.total_reads_kept = 0
        self.description_file_path = A('description')

        if self.fetch_filter:
            valid_fetch_filters = [k for k in constants.fetch_filters.keys() if k]

            if self.fetch_filter not in valid_fetch_filters:
                raise ConfigError(f"Your fetch filter '{self.fetch_filter}' is not among those anvi'o knows about :/ If you "
                                  f"would like to try again with a different one, here is a list for your consideration: "
                                  f"{', '.join(valid_fetch_filters)}.")

        # these are the views this class will be filling in:
        self.essential_data_fields_for_anvio_profiles = copy.deepcopy(constants.essential_data_fields_for_anvio_profiles)

        # but the views should not include the view `variability` if SNVs are not going to be profiled:
        if 'variability' in self.essential_data_fields_for_anvio_profiles and self.skip_SNV_profiling:
            self.essential_data_fields_for_anvio_profiles.pop(self.essential_data_fields_for_anvio_profiles.index('variability'))

        # make sure early on that both the distance and linkage is OK.
        clustering.is_distance_and_linkage_compatible(self.distance, self.linkage)

        # whether the profile database is a blank (without any BAM files or reads):
        self.blank = A('blank_profile')

        if self.skip_SNV_profiling and self.skip_edges > 0:
            raise ConfigError(f"You can't ask anvi'o to skip profiling of SNVs and also ask to ignore {self.skip_edges} "
                              f"nucleotides from the beginning and the end of short reads while profiling of SNVs. You "
                              f"either need to drop the `--skip-SNV-profiling` flag, or the `--skip-edges` flag :/")

        if not self.blank and self.contigs_shall_be_clustered and self.skip_hierarchical_clustering:
            raise ConfigError("You are confused, and confusing anvi'o, too. You can't as hierarchical clustering "
                              "to be performed with one flag, and try to skip it with another one :(")

        if self.blank and self.contigs_shall_be_clustered and self.skip_hierarchical_clustering:
            raise ConfigError("So you want to generate a blank profile, and you both want hierarchical clustering "
                              "of your contigs to be performed, and skipped. No.")

        if self.blank and self.contigs_shall_be_clustered:
            raise ConfigError("When the blank profile is asked to be generated, there is no need to ask for the "
                              "hierarchical clustering of contigs. It is going to be done by default. If it is "
                              "not changing anything, why is anvi'o upset with you? Because. Let's don't use flags "
                              "we don't need.")

        if self.blank and not self.skip_hierarchical_clustering:
            self.contigs_shall_be_clustered = True

        if A('contigs_of_interest'):
            filesnpaths.is_file_exists(args.contigs_of_interest)
            self.contig_names_of_interest = set([c.strip() for c in open(args.contigs_of_interest).readlines()\
                                                                           if c.strip() and not c.startswith('#')])
        else:
            self.contig_names_of_interest = set([])

        if self.list_contigs_and_exit:
            self.list_contigs()
            sys.exit()

        if not self.contigs_db_path:
            raise ConfigError("No contigs database, no profilin'. Bye.")

        # store our run object. 'why?', you may ask. well, keep reading.
        my_run = self.run

        # Initialize contigs db
        dbops.ContigsSuperclass.__init__(self, self.args, r=null_run, p=self.progress)
        self.init_contig_sequences(contig_names_of_interest=self.contig_names_of_interest)
        self.contig_names_in_contigs_db = set(self.contigs_basic_info.keys())

        # restore our run object. OK. take deep breath. because you are a good programmer, you have
        # this voice in your head tellin gyou that the tinkering with self.run here feels kind of out
        # of place. that voice is right, but the voice doesn't know the struggles of poor souls that
        # had to resort to a solution like this. You see, we don't want to see any run messages from
        # ContigsSuper in profiler output. But when we pass a `run=null_run` to the class, due to
        # inheritance, it also modifies our own `self.run` with the null one, making the profilesuper
        # go all quiet. so here we basically need to re-engage our `self.run`. But then if the user
        # actually ASKED for ProfileSuper to be quiet, then we can't simply just inherit another Run
        # object and move on with our lives in a single line right here. We in fact need to store our
        # previous run object and then re-engage that one instead of a brand new one. There you have
        # it (yes. funny detail, no one cares, but I do because I have no friends -- IF YOU ARE
        # READING THESE LINES YOU ARE AUTOMATICALLY MY FRIEND, so SORRY):
        self.run = my_run

        self.bam = None
        self.contigs = []

        self.database_paths = {'CONTIGS.db': os.path.abspath(self.contigs_db_path)}

        self.profile_db_path = None

        self.clustering_configs = constants.clustering_configs['blank' if self.blank else 'single']

        # if genes are not called, yet the user is asking for codon frequencies to be profiled, we give
        # a warning and force-turn that flag off.
        if (not self.a_meta['genes_are_called']) and self.profile_SCVs:
            self.run.warning("You asked the codon frequencies to be profiled, but genes were not called "
                             "for your contigs database. Anvi'o is assigning `False` to the --profile-SCVs "
                             "flag, overruling your request like a boss.")
            self.profile_SCVs = False

        # we don't know what we are about
        self.description = None

        # additional layer data will be filled later
        self.layer_additional_keys = []
        self.layer_additional_data = {}


    def init_dirs_and_dbs(self):
        if not self.contigs_db_path:
            raise ConfigError("You can not run profiling without a contigs database. You can create "
                               "one using 'anvi-gen-contigs-database'. Not sure how? Please see the "
                               "tutorial: http://merenlab.org/2015/05/02/anvio-tutorial/")

        if self.description_file_path:
            filesnpaths.is_file_plain_text(self.description_file_path)
            self.description = open(os.path.abspath(self.description_file_path), 'r').read()

        if self.output_directory:
            self.output_directory = filesnpaths.check_output_directory(self.output_directory, ok_if_exists=self.overwrite_output_destinations)
        else:
            output_dir_path = os.path.dirname(os.path.abspath(self.input_file_path))

            if self.sample_id:
                self.output_directory = filesnpaths.check_output_directory(os.path.join(output_dir_path, self.sample_id), ok_if_exists=self.overwrite_output_destinations)
            else:
                raise ConfigError("There is no `self.sample_id`, there is no `self.output_directory` :/ Anvi'o needs an adult :(")

        self.progress.new('Initializing')

        self.progress.update('Creating the output directory ...')
        filesnpaths.gen_output_directory(self.output_directory, self.progress, delete_if_exists=self.overwrite_output_destinations)

        self.progress.update('Creating a new single profile database with contigs hash "%s" ...' % self.a_meta['contigs_db_hash'])
        self.profile_db_path = self.generate_output_destination('PROFILE.db')
        profile_db = dbops.ProfileDatabase(self.profile_db_path)

        if self.blank:
            # if we are about to generate a blank profile, there is no
            # SNV profiling
            self.skip_SNV_profiling = True

        if self.skip_SNV_profiling:
            self.profile_SCVs = False
            self.skip_INDEL_profiling = True
            self.skip_edges = 0

        meta_values = {'db_type': 'profile',
                       'anvio': __version__,
                       'sample_id': self.sample_id,
                       'samples': self.sample_id,
                       'merged': False,
                       'blank': self.blank,
                       'items_ordered': False,
                       'default_view': 'mean_coverage',
                       'min_contig_length': self.min_contig_length,
                       'max_contig_length': self.max_contig_length,
                       'SNVs_profiled': not self.skip_SNV_profiling,
                       'SCVs_profiled': self.profile_SCVs,
                       'INDELs_profiled': not self.skip_INDEL_profiling,
                       'min_percent_identity': self.min_percent_identity or 0,
                       'fetch_filter': self.fetch_filter,
                       'min_coverage_for_variability': self.min_coverage_for_variability,
                       'report_variability_full': self.report_variability_full,
                       'skip_edges_for_variant_profiling': self.skip_edges,
                       'contigs_db_hash': self.a_meta['contigs_db_hash'],
                       'description': self.description if self.description else '_No description is provided_'}
        profile_db.create(meta_values)

        self.progress.update('Creating a new auxiliary database with contigs hash "%s" ...' % self.a_meta['contigs_db_hash'])
        self.auxiliary_db_path = self.generate_output_destination('AUXILIARY-DATA.db')
        self.auxiliary_db = auxiliarydataops.AuxiliaryDataForSplitCoverages(self.auxiliary_db_path,
                                                                            self.a_meta['contigs_db_hash'],
                                                                            create_new=True,
                                                                            run=null_run,
                                                                            progress=null_progress)

        self.progress.end()

        if not self.skip_SNV_profiling:
            self.variable_nts_table = TableForVariability(self.profile_db_path, progress=null_progress)

        if self.profile_SCVs:
            self.variable_codons_table = TableForCodonFrequencies(self.profile_db_path, progress=null_progress)

        if not self.skip_INDEL_profiling:
            self.indels_table = TableForIndels(self.profile_db_path, progress=null_progress)


    def _run(self):
        self.check_args()

        self.set_sample_id()

        self.init_dirs_and_dbs()

        self.run.log_file_path = self.generate_output_destination('RUNLOG.txt')
        self.run.info('Sample name set', self.sample_id)
        self.run.info('Description', 'Found (%d characters)' % len(self.description) if self.description else None)
        self.run.info('Profile DB path', self.profile_db_path, display_only=True)
        self.run.info('Contigs DB path', self.contigs_db_path)
        self.run.info('Contigs DB hash', self.a_meta['contigs_db_hash'])
        self.run.info('Command line', utils.get_cmd_line(), align_long_values=False, nl_after=1)

        self.run.info('Minimum percent identity of reads to be profiled', self.min_percent_identity, mc='green')
        self.run.info('Fetch filter engaged', self.fetch_filter, mc='green', nl_after=1)

        self.run.info('Is merged profile?', False)
        self.run.info('Is blank profile?', self.blank)
        self.run.info('Skip contigs shorter than', self.min_contig_length)
        self.run.info('Skip contigs longer than', self.max_contig_length)
        self.run.info('Perform hierarchical clustering of contigs?', self.contigs_shall_be_clustered, nl_after=1)

        self.run.info('Profile single-nucleotide variants (SNVs)?', not self.skip_SNV_profiling)
        self.run.info('Ancient DNA friendly profiling?', 'True' if self.skip_edges else 'False')
        if self.skip_edges:
            self.run.info(' - How many edge nts ignore for SNV profiling?', self.skip_edges, mc='red')
        self.run.info('Profile single-codon variants (SCVs/+SAAVs)?', self.profile_SCVs)
        self.run.info('Profile insertion/deletions (INDELs)?', not self.skip_INDEL_profiling)
        self.run.info('Minimum coverage to calculate SNVs', self.min_coverage_for_variability)
        self.run.info('Report FULL variability data?', self.report_variability_full)

        self.run.warning("Your minimum contig length is set to %s base pairs. So anvi'o will not take into "
                         "consideration anything below that. If you need to kill this an restart your "
                         "analysis with another minimum contig length value, feel free to press CTRL+C." \
                                                % (pp(self.min_contig_length)))

        if self.max_contig_length < sys.maxsize:
            self.run.warning("Your maximum contig length is set to %s base pairs. Which means anvi'o will remove "
                             "any contigs that are longer than this value." % pp(self.max_contig_length))

        # this is kinda important. we do not run full-blown profile function if we are dealing with a summarized
        # profile...
        if self.blank:
            self.init_mock_profile()

            # creating a null view_data_splits dict:
            for view in self.essential_data_fields_for_anvio_profiles:
                for table_name in [f"{view}_{target}" for target in ['splits', 'contigs']]:
                    TablesForViews(self.profile_db_path).remove(view, table_names_to_blank=[table_name])
                    TablesForViews(self.profile_db_path,
                                   progress=self.progress) \
                                        .create_new_view(view_data=[],
                                                         table_name=table_name,
                                                         skip_sanity_check=True,
                                                         view_name=None if table_name.endswith('contigs') else view)
        elif self.input_file_path:
            self.init_profile_from_BAM()
            if self.num_threads > 1 or self.args.force_multi:
                self.profile_multi_thread()
            else:
                self.profile_single_thread()
        else:
            raise ConfigError("What are you doing? :( Whatever it is, anvi'o will have none of it.")

        # update layer additional data table content
        if self.layer_additional_data:
            self.progress.new("Additional layer data")
            self.progress.update("Updating the profile db ...")
            layer_additional_data_table = TableForLayerAdditionalData(argparse.Namespace(profile_db=self.profile_db_path), r=null_run, p=null_progress)
            layer_additional_data_table.add({self.sample_id: self.layer_additional_data}, self.layer_additional_keys)
            self.progress.end()

            self.run.info("Additional data added to the new profile DB", f"{', '.join(self.layer_additional_keys)}", nl_before=1)

        if self.contigs_shall_be_clustered:
            self.cluster_contigs()

        if self.bam:
            self.bam.close()

        self.run.info_single('Happy 😇', nl_before=1, nl_after=1)

        self.run.quit()


    def generate_variable_codons_table(self):
        if self.skip_SNV_profiling or not self.profile_SCVs:
            return

        for contig in self.contigs:
            for split in contig.splits:
                for gene_callers_id in split.SCV_profiles:
                    # We reorder to the profiles in the order they will appear in the output table
                    split.SCV_profiles[gene_callers_id] = OrderedDict(
                        [(col, split.SCV_profiles[gene_callers_id][col]) for col in t.variable_codons_table_structure]
                    )

                    entries = zip(*split.SCV_profiles[gene_callers_id].values())
                    for entry in entries:
                        self.variable_codons_table.append(entry)

        self.variable_codons_table.store()


    def generate_variable_nts_table(self):
        if self.skip_SNV_profiling:
            return

        for contig in self.contigs:
            for split in contig.splits:
                if split.num_SNV_entries == 0:
                    continue

                # We reorder to the profiles in the order they will appear in the output table
                split.SNV_profiles = OrderedDict(
                    [(col, split.SNV_profiles[col]) for col in t.variable_nts_table_structure]
                )

                entries = zip(*split.SNV_profiles.values())
                for entry in entries:
                    self.variable_nts_table.append(entry)

        self.variable_nts_table.store()


    def generate_indels_table(self):
        if self.skip_INDEL_profiling:
            return

        for contig in self.contigs:
            for split in contig.splits:
                for entry in split.INDEL_profiles.values():
                    self.indels_table.append([self.sample_id] + list(entry.values()))

        self.indels_table.store()


    def store_split_coverages(self):
        for contig in self.contigs:
            for split in contig.splits:
                self.auxiliary_db.append(split.name, self.sample_id, split.coverage.c)

        self.auxiliary_db.store()


    def set_sample_id(self):
        if self.sample_id:
            utils.check_sample_id(self.sample_id)
        else:
            if self.input_file_path:
                self.input_file_path = os.path.abspath(self.input_file_path)
                self.sample_id = os.path.splitext(os.path.basename(self.input_file_path))[0]
                self.sample_id = self.sample_id.replace('-', '_').replace('.', '_').replace(' ', '_')

                if self.sample_id[0] in constants.digits:
                    self.sample_id = 's' + self.sample_id

                if self.fetch_filter:
                    self.sample_id = f"{self.sample_id}_{self.fetch_filter.upper().replace('-', '_').replace('.', '_').replace(' ', '_')}"

                utils.check_sample_id(self.sample_id)
            if self.serialized_profile_path:
                self.serialized_profile_path = os.path.abspath(self.serialized_profile_path)
                self.sample_id = os.path.basename(os.path.dirname(self.serialized_profile_path))


    def check_contigs_without_any_gene_calls(self, contig_names):
        if not self.a_meta['genes_are_called']:
            self.run.warning("The contigs database '%s' does not contain any gene calls. Which means the profiling step "
                             "will not be able to characterize 'gene coverages'. If you are OK with this, anvi'o will be "
                             "OK with it as well." % (self.contigs_db_path))
            return

        contig_names = set(contig_names)
        contigs_without_any_gene_calls = [c for c in contig_names if c not in self.contig_name_to_genes]

        if len(contigs_without_any_gene_calls):
            import random
            P = lambda x: 'are %d contigs' % (x) if x > 1 else 'there is one contig'
            self.run.warning('According to the data generated in the contigs database, there %s in your BAM file '
                             'with 0 gene calls. Which may not be unusual if (a) some of your contigs are very short, '
                             'or (b) your the gene caller was not capable of dealing with the type of data you had. '
                             'If you would like to take a look yourself, here is one contig that is missing any genes: %s"' %\
                                      (P(len(contigs_without_any_gene_calls)), random.choice(contigs_without_any_gene_calls)))


    def list_contigs(self):
        import signal
        signal.signal(signal.SIGPIPE, signal.SIG_DFL)

        self.progress.new('Init')
        self.progress.update('Reading BAM File')
        self.bam = bamops.BAMFileObject(self.input_file_path, 'rb')
        self.bam.fetch_filter = self.fetch_filter
        self.progress.end()

        self.contig_names = self.bam.references
        self.contig_lengths = self.bam.lengths

        utils.check_contig_names(self.contig_names)

        for tpl in sorted(zip(self.contig_lengths, self.contig_names), reverse=True):
            print('%-40s %s' % (tpl[1], pp(int(tpl[0]))))


    def remove_contigs_based_on_min_max_contig_length(self):
        """Removes contigs that are shorter or longer than user defined values"""

        contigs_with_good_lengths = set()

        for i in range(0, len(self.contig_names)):
            if self.contig_lengths[i] >= self.min_contig_length and self.contig_lengths[i] <= self.max_contig_length:
                contigs_with_good_lengths.add(i)

        if not len(contigs_with_good_lengths):
            try:
                filesnpaths.shutil.rmtree(self.output_directory)
            except:
                pass

            raise ConfigError("Anvi'o applied your min/max length criteria for contigs to filter out the bad ones "
                              "and has bad news: not a single contig in your contigs database was greater than %s "
                              "and smaller than %s nts :( So this profiling attempt did not really go anywhere. "
                              "Please remove your half-baked output directory if it is still there: '%s'." \
                                        % (pp(self.min_contig_length), pp(self.max_contig_length), self.output_directory))
        else:
            self.contig_names = [self.contig_names[i] for i in contigs_with_good_lengths]
            self.contig_lengths = [self.contig_lengths[i] for i in contigs_with_good_lengths]
            self.num_contigs = len(self.contig_names)    # we will store these two
            self.total_length = sum(self.contig_lengths) # into the db in a second.

        contigs_with_good_lengths = set(self.contig_names) # for fast access
        self.split_names = set([])
        self.contig_name_to_splits = {}
        for split_name in sorted(self.splits_basic_info.keys()):
            parent = self.splits_basic_info[split_name]['parent']

            if parent not in contigs_with_good_lengths:
                continue

            self.split_names.add(split_name)

            if parent in self.contig_name_to_splits:
                self.contig_name_to_splits[parent].append(split_name)
            else:
                self.contig_name_to_splits[parent] = [split_name]

        self.num_splits = len(self.split_names)


    def init_profile_from_BAM(self):
        filesnpaths.is_file_bam_file(self.input_file_path)

        self.progress.new('Init')
        self.progress.update('Reading BAM File')

        try:
            self.bam = bamops.BAMFileObject(self.input_file_path)
        except Exception as e:
            raise ConfigError(f"Sorry, this BAM file does not look like a BAM file :( Here is "
                              f"the complaint coming from the depths of the codebase: '{e}'.")

        self.bam.fetch_filter = self.fetch_filter
        self.num_reads_mapped = self.bam.mapped
        self.progress.end()

        self.contig_names = self.bam.references
        self.contig_lengths = self.bam.lengths

        if self.min_percent_identity:
            # We will permit as many as 1% of sampled reads to be missing an MD tag, for whatever
            # crazy mapping reason, although it is most probable that if 1 read is missing an MD
            # tag, they all are.
            requirement = lambda x: sum(x) / len(x) >= 0.99

            if not self.bam.reads_have_MD_tags(require=requirement):
                raise ConfigError("Your BAM file has reads that do not have MD "
                                  "tags, which give essential information about the "
                                  "reference bases in alignments. Unfortunately, "
                                  "this is required to use --min-percent-identity. "
                                  "Please remove this flag or redo your mapping.")

        utils.check_contig_names(self.contig_names)

        self.run.info('Input BAM', self.input_file_path)
        self.run.info('Output directory path', self.output_directory, display_only=True, nl_after=1)

        self.run.info('Number of reads in the BAM file', pp(int(self.num_reads_mapped)))
        self.run.info('Number of sequences in the contigs DB', pp(len(self.contig_names)))

        if self.contig_names_of_interest:
            indexes = [self.contig_names.index(r) for r in self.contig_names_of_interest if r in self.contig_names]
            self.contig_names = [self.contig_names[i] for i in indexes]
            self.contig_lengths = [self.contig_lengths[i] for i in indexes]
            self.run.info('Number of contigs selected for analysis', pp(len(self.contig_names)), mc='green')

        # it brings good karma to let the user know what the hell is wrong with their data:
        self.check_contigs_without_any_gene_calls(self.contig_names)

        # check for the -M parameter.
        self.remove_contigs_based_on_min_max_contig_length()

        # let's see whether the user screwed up to follow the simple instructions
        # mentioned here: http://merenlab.org/2015/05/01/anvio-tutorial/#preparation
        for contig_name in self.contig_names:
            if contig_name not in self.contig_names_in_contigs_db:
                raise ConfigError("At least one contig name in your BAM file does not match contig names stored in the "
                                   "contigs database. For instance, this is one contig name found in your BAM file: '%s', "
                                   "and this is another one found in your contigs database: '%s'. You may be using an "
                                   "contigs database for profiling that has nothing to do with the BAM file you are "
                                   "trying to profile, or you may have failed to fix your contig names in your FASTA file "
                                   "prior to mapping, which is described here: %s"\
                                        % (contig_name, self.contig_names_in_contigs_db.pop(), 'http://goo.gl/Q9ChpS'))

        self.run.info('Number of contigs to be conisdered (after -M)', self.num_contigs, display_only=True)
        self.run.info('Number of splits', self.num_splits)
        self.run.info('Number of nucleotides', self.total_length)

        profile_db = dbops.ProfileDatabase(self.profile_db_path, quiet=True)
        profile_db.db.set_meta_value('num_splits', self.num_splits)
        profile_db.db.set_meta_value('num_contigs', self.num_contigs)
        profile_db.db.set_meta_value('total_length', self.total_length)
        profile_db.disconnect()

        self.layer_additional_data['total_reads_mapped'] = self.num_reads_mapped
        self.layer_additional_keys.append('total_reads_mapped')


    def init_mock_profile(self):
        self.progress.new('Init')
        self.progress.update('...')
        self.num_reads_mapped = 0
        self.progress.end()

        self.contig_names = list(self.contigs_basic_info.keys())
        self.contig_lengths = [self.contigs_basic_info[contig_name]['length'] for contig_name in self.contigs_basic_info]
        self.total_length = sum(self.contig_lengths)
        self.num_contigs = len(self.contig_names)

        utils.check_contig_names(self.contig_names)

        self.run.info('input_bam', None)
        self.run.info('output_dir', self.output_directory, display_only=True)

        # check for the -M parameter.
        self.remove_contigs_based_on_min_max_contig_length()

        self.run.info('num_contigs_after_M', self.num_contigs, display_only=True)
        self.run.info('num_contigs', self.num_contigs, quiet=True)
        self.run.info('num_splits', self.num_splits)
        self.run.info('total_length', self.total_length)

        profile_db = dbops.ProfileDatabase(self.profile_db_path, quiet=True)
        profile_db.db.set_meta_value('num_splits', self.num_splits)
        profile_db.db.set_meta_value('num_contigs', self.num_contigs)
        profile_db.db.set_meta_value('total_length', self.total_length)
        profile_db.disconnect()

        self.layer_additional_data['total_reads_mapped'] = self.num_reads_mapped
        self.layer_additional_keys.append('total_reads_mapped')


    def generate_output_destination(self, postfix, directory=False):
        return_path = os.path.join(self.output_directory, postfix)

        if directory == True:
            if os.path.exists(return_path):
                shutil.rmtree(return_path)
            os.makedirs(return_path)

        return return_path


    def populate_gene_info_for_splits(self, contig):
        """Stores info available to the ContigsSuperClass into contigops.Split objects

        Populates Split.per_position_info as a dictionary of arrays, each with length equal to the
        split.
        """

        if self.skip_SNV_profiling:
            return

        info_of_interest = [
            'corresponding_gene_call',
            'codon_order_in_gene',
            'base_pos_in_codon',
            'in_noncoding_gene_call',
            'in_coding_gene_call',
        ]

        if self.profile_SCVs:
            # NOTE The people want SCVs, and that means Split needs to know about gene calls. Rather
            #      than giving each split a dictionary of gene calls, in an ugly but much faster
            #      fashion, we add these 3 pieces of information which provide sufficient gene call
            #      information to calculate SCV. Previously, we parsed the genes_in_contigs dict
            #      which is fine for small databases, but the cost is enormous for very large
            #      (100,000+ contigs) databases. To give perspective, it was about 1s per contig,
            #      which is 30 hours for 100,000 contigs. That code elegantly looked like this:
            #
            #      for split in contig.splits:
            #          gene_ids_in_split = set(gc['gene_callers_id']
            #                                  for key, gc in self.genes_in_splits.items()
            #                                  if key in self.split_name_to_genes_in_splits_entry_ids[split.name])
            #          split.gene_calls = {gene_id: self.genes_in_contigs_dict[gene_id] for gene_id in gene_ids_in_split}
            info_of_interest.extend([
                'forward',
                'gene_start',
                'gene_stop',
            ])

        nt_info = self.get_gene_info_for_each_position(contig.name, info=info_of_interest)

        for split in contig.splits:
            for info in info_of_interest:
                split.per_position_info[info] = nt_info[info][split.start:split.end]


    @staticmethod
    def profile_contig_worker(self, available_index_queue, output_queue):
        import traceback

        bam_file = bamops.BAMFileObject(self.input_file_path)
        bam_file.fetch_filter = self.fetch_filter

        # Attach to shared memory for contig sequences (created by main process)
        # This avoids duplicating sequence data across all worker processes
        if hasattr(self, '_shared_sequence_shm_name') and self._shared_sequence_shm_name:
            try:
                self._worker_sequence_store = streamingops.SharedSequenceStore.from_existing(
                    self._shared_sequence_shm_name,
                    self._shared_sequence_index
                )
                # Replace contig_sequences with a proxy that fetches from shared memory
                # This makes all code that accesses self.contig_sequences work transparently
                self.contig_sequences = self._worker_sequence_store.as_dict_proxy()
            except Exception as e:
                output_queue.put(ConfigError(f"Worker failed to attach to shared memory '{self._shared_sequence_shm_name}': {e}"))
                return
        else:
            # This should not happen in multi-threaded mode - shared memory should always be set up
            output_queue.put(ConfigError(f"Worker missing shared memory attributes. "
                                         f"Has _shared_sequence_shm_name: {hasattr(self, '_shared_sequence_shm_name')}, "
                                         f"Value: {getattr(self, '_shared_sequence_shm_name', 'NOT SET')}"))
            return

        # Attach to shared memory for nt_positions_info (created by main process)
        # This avoids duplicating the large per-nucleotide arrays across workers
        if hasattr(self, '_shared_nt_positions_shm_name'):
            try:
                self._worker_nt_positions_store = streamingops.SharedNtPositionsStore.from_existing(
                    self._shared_nt_positions_shm_name,
                    self._shared_nt_positions_index
                )
                # Replace nt_positions_info with a proxy that fetches from shared memory
                # LazyProperty stores cached data in _lazy_loaded_data[attr_name],
                # so we need to set it there to override the lazy loading
                if not hasattr(self, '_lazy_loaded_data'):
                    self._lazy_loaded_data = {}
                self._lazy_loaded_data['nt_positions_info'] = self._worker_nt_positions_store.as_dict_proxy()
            except Exception as e:
                output_queue.put(ConfigError(f"Worker failed to attach to nt_positions shared memory: {e}"))
                return

        while True:
            try:
                index = available_index_queue.get(True)
                contig_name = self.contig_names[index]
                contig_length = self.contig_lengths[index]

                contig = self.process_contig(bam_file, contig_name, contig_length)
                output_queue.put(contig)

                if contig is not None:
                    # We mark these for deletion the next time garbage is collected
                    for split in contig.splits:
                        del split.coverage
                        del split.auxiliary
                        del split
                    del contig.splits[:]
                    del contig.coverage
                    del contig
            except Exception as e:
                # This thread encountered an error. We send the error back to the main thread which
                # will terminate the job. Include traceback for debugging.
                error_msg = f"{type(e).__name__}: {e}\n\nTraceback:\n{traceback.format_exc()}"
                output_queue.put(ConfigError(error_msg))

        # we are closing this object here for clarity, although we are not really closing it since
        # the code never reaches here and the worker is killed by its parent:
        bam_file.close()
        return


    def process_contig(self, bam_file, contig_name, contig_length):
        """Process a contig using single-pass streaming approach.

        This method uses a single-pass design that extracts coverage, SNVs, and
        INDELs from each read in one iteration through the BAM file. This dramatically
        reduces memory usage compared to the previous multi-pass approach.

        If processing_chunk_size is set, splits larger than the chunk size will be
        processed in smaller chunks to further reduce peak memory usage.
        """
        from anvio.variability import (VariablityTestFactory, ProcessNucleotideCounts,
                                        ProcessCodonCounts, ProcessIndelCounts)

        timer = terminal.Timer(initial_checkpoint_key='Start')

        # Initialize contig (same as process_contig)
        contig = contigops.Contig(contig_name)
        contig.length = contig_length
        contig.split_length = self.a_meta['split_length']
        contig.skip_SNV_profiling = self.skip_SNV_profiling
        timer.make_checkpoint('Initialization done')

        # Create split objects (same as process_contig)
        for split_name in self.contig_name_to_splits[contig_name]:
            s = self.splits_basic_info[split_name]
            # Get split sequence - works with both regular dict and shared memory proxy
            split_sequence = self.contig_sequences[contig_name]['sequence'][s['start']:s['end']]
            split = contigops.Split(split_name, split_sequence, contig_name, s['order_in_parent'], s['start'], s['end'])
            contig.splits.append(split)

        timer.make_checkpoint('Split objects initialized')

        # Populate gene info for splits (same as process_contig)
        self.populate_gene_info_for_splits(contig)
        timer.make_checkpoint('Gene info for split added')

        # Set up variability test class
        variability_test_class_null = VariablityTestFactory(params=None)
        variability_test_class_default = VariablityTestFactory(params={'b': 2, 'm': 1.45, 'c': 0.05})

        # Initialize contig-level coverage array for later slicing
        contig_coverage = np.zeros(contig_length, dtype=np.int32)

        # Process each split
        for split in contig.splits:
            # Determine chunk size: use processing_chunk_size if set and smaller than split
            chunk_size = split.length
            if self.processing_chunk_size and self.processing_chunk_size < split.length:
                chunk_size = self.processing_chunk_size

            # Prepare per-position data for SNV annotation
            additional_per_position_data = {}
            if not self.skip_SNV_profiling and hasattr(split, 'per_position_info'):
                additional_per_position_data = split.per_position_info.copy()

            # Decide read iterator based on percent identity filtering
            if self.min_percent_identity:
                read_iterator = bam_file.fetch_filter_and_trim
                kwargs = {'percent_id_cutoff': self.min_percent_identity}
            else:
                read_iterator = bam_file.fetch_and_trim
                kwargs = {}

            if chunk_size >= split.length:
                # Process entire split at once (standard case)
                coverage_acc = streamingops.CoverageAccumulator(split.length)
                snv_acc = streamingops.SNVAccumulator(split.length) if not self.skip_SNV_profiling else None
                indel_acc = streamingops.INDELAccumulator(split.length) if not self.skip_INDEL_profiling else None

                for read in read_iterator(split.parent, split.start, split.end, **kwargs):
                    evidence = read.extract_variant_evidence(
                        split_start=split.start,
                        skip_SNV=self.skip_SNV_profiling,
                        skip_INDEL=self.skip_INDEL_profiling
                    )
                    coverage_acc.update(evidence['coverage_blocks'])
                    if snv_acc is not None:
                        snv_acc.update_from_vectorized(evidence['vectorized'], split.start)
                    if indel_acc is not None:
                        indel_acc.update(
                            evidence['insertions'],
                            evidence['deletions'],
                            split.name,
                            split.sequence,
                            additional_per_position_data if additional_per_position_data else None
                        )
            else:
                # Process split in chunks to reduce memory usage
                # Initialize arrays/dicts to accumulate chunk results
                full_coverage = np.zeros(split.length, dtype=np.int32)
                full_snv_counts = np.zeros((5, split.length), dtype=np.int32) if not self.skip_SNV_profiling else None
                all_indels = {} if not self.skip_INDEL_profiling else None
                total_reads = 0

                for chunk_start in range(0, split.length, chunk_size):
                    chunk_end = min(chunk_start + chunk_size, split.length)
                    chunk_len = chunk_end - chunk_start

                    # Absolute positions in contig coordinates
                    abs_chunk_start = split.start + chunk_start
                    abs_chunk_end = split.start + chunk_end

                    # Create chunk-sized accumulators
                    chunk_cov_acc = streamingops.CoverageAccumulator(chunk_len)
                    chunk_snv_acc = streamingops.SNVAccumulator(chunk_len) if not self.skip_SNV_profiling else None
                    chunk_indel_acc = streamingops.INDELAccumulator(chunk_len) if not self.skip_INDEL_profiling else None

                    # Fetch reads overlapping this chunk
                    for read in read_iterator(split.parent, abs_chunk_start, abs_chunk_end, **kwargs):
                        evidence = read.extract_variant_evidence(
                            split_start=abs_chunk_start,  # Positions relative to chunk start
                            skip_SNV=self.skip_SNV_profiling,
                            skip_INDEL=self.skip_INDEL_profiling
                        )

                        # Clip coverage blocks to chunk boundaries
                        clipped_blocks = []
                        for block_start, block_end in evidence['coverage_blocks']:
                            # Clip to [0, chunk_len)
                            clipped_start = max(0, block_start)
                            clipped_end = min(chunk_len, block_end)
                            if clipped_start < clipped_end:
                                clipped_blocks.append((clipped_start, clipped_end))
                        chunk_cov_acc.update(clipped_blocks)

                        # Update SNV accumulator (with chunk-relative positions)
                        if chunk_snv_acc is not None:
                            chunk_snv_acc.update_from_vectorized(
                                evidence['vectorized'],
                                abs_chunk_start,
                                position_limit=chunk_len
                            )

                        # Update INDEL accumulator (filter to chunk boundaries)
                        if chunk_indel_acc is not None:
                            # Filter insertions/deletions to chunk boundaries
                            chunk_insertions = [(p, s) for p, s in evidence['insertions'] if 0 <= p < chunk_len]
                            chunk_deletions = [(p, l) for p, l in evidence['deletions'] if 0 <= p < chunk_len]
                            if chunk_insertions or chunk_deletions:
                                chunk_indel_acc.update(
                                    chunk_insertions,
                                    chunk_deletions,
                                    split.name,
                                    split.sequence[chunk_start:chunk_end],
                                    {k: v[chunk_start:chunk_end] for k, v in additional_per_position_data.items()} if additional_per_position_data else None
                                )

                    # Accumulate chunk results into full arrays
                    chunk_result = chunk_cov_acc.finalize()
                    full_coverage[chunk_start:chunk_end] = chunk_result['coverage']
                    total_reads += chunk_result['num_reads']

                    if chunk_snv_acc is not None:
                        chunk_snv_counts = chunk_snv_acc.finalize()
                        full_snv_counts[:, chunk_start:chunk_end] = chunk_snv_counts

                    if chunk_indel_acc is not None:
                        chunk_indels = chunk_indel_acc.finalize()
                        # Adjust positions and merge into all_indels
                        for key, indel in chunk_indels.items():
                            # Create new key with adjusted position
                            adjusted_pos = indel['pos'] + chunk_start
                            new_key = (indel['type'], adjusted_pos, indel.get('sequence', ''), indel['length'])
                            if new_key in all_indels:
                                all_indels[new_key]['count'] += indel['count']
                            else:
                                adjusted_indel = indel.copy()
                                adjusted_indel['pos'] = adjusted_pos
                                all_indels[new_key] = adjusted_indel

                # Create mock accumulators with combined results for compatibility
                # Note: lambdas need 'self' parameter since they become instance methods
                coverage_acc = type('MockCoverageAcc', (), {
                    'finalize': lambda self: {
                        'coverage': full_coverage,
                        'num_reads': total_reads,
                        'min': int(full_coverage.min()) if len(full_coverage) > 0 else 0,
                        'max': int(full_coverage.max()) if len(full_coverage) > 0 else 0,
                        'mean': float(full_coverage.mean()) if len(full_coverage) > 0 else 0.0,
                        'std': float(full_coverage.std()) if len(full_coverage) > 0 else 0.0,
                        'median': float(np.median(full_coverage)) if len(full_coverage) > 0 else 0.0,
                        'detection': float(np.count_nonzero(full_coverage) / len(full_coverage)) if len(full_coverage) > 0 else 0.0,
                    }
                })()
                snv_acc = type('MockSNVAcc', (), {'finalize': lambda self: full_snv_counts})() if full_snv_counts is not None else None
                indel_acc = type('MockINDELAcc', (), {'finalize': lambda self: all_indels})() if all_indels is not None else None

            # Finalize coverage
            cov_result = coverage_acc.finalize()
            split.coverage = bamops.Coverage()
            split.coverage.c = cov_result['coverage']
            split.coverage.min = cov_result['min']
            split.coverage.max = cov_result['max']
            split.coverage.mean = cov_result['mean']
            split.coverage.std = cov_result['std']
            split.coverage.median = cov_result['median']
            split.coverage.detection = cov_result['detection']
            split.coverage.num_reads = cov_result['num_reads']

            # Compute outliers (matching original behavior)
            if len(split.coverage.c) > 0:
                mean_cov = split.coverage.mean
                std_cov = split.coverage.std
                if std_cov > 0:
                    split.coverage.is_outlier = np.abs(split.coverage.c - mean_cov) > (2.5 * std_cov)
                else:
                    split.coverage.is_outlier = np.zeros(len(split.coverage.c), dtype=bool)
                split.coverage.is_outlier_in_parent = split.coverage.is_outlier.copy()
            else:
                split.coverage.is_outlier = np.array([], dtype=bool)
                split.coverage.is_outlier_in_parent = np.array([], dtype=bool)

            # Compute mean_Q2Q3 (matching original behavior)
            if len(split.coverage.c) >= 4:
                sorted_cov = np.sort(split.coverage.c)
                q1_idx = len(sorted_cov) // 4
                q3_idx = (3 * len(sorted_cov)) // 4
                split.coverage.mean_Q2Q3 = np.mean(sorted_cov[q1_idx:q3_idx])
            else:
                split.coverage.mean_Q2Q3 = split.coverage.mean

            # Copy coverage to contig array
            contig_coverage[split.start:split.end] = split.coverage.c

            # Process SNVs if not skipped
            if snv_acc is not None:
                allele_counts_array = snv_acc.finalize()

                # Add coverage outlier info to additional data
                additional_per_position_data['cov_outlier_in_split'] = split.coverage.is_outlier.astype(int)
                additional_per_position_data['cov_outlier_in_contig'] = split.coverage.is_outlier_in_parent.astype(int)

                # Use existing ProcessNucleotideCounts for compatibility
                split_as_index = utils.nt_seq_to_nt_num_array(split.sequence)
                nt_to_array_index = {nt: i for i, nt in enumerate(constants.nucleotides)}

                nt_profile = ProcessNucleotideCounts(
                    allele_counts=allele_counts_array,
                    allele_to_array_index=nt_to_array_index,
                    sequence=split.sequence,
                    sequence_as_index=split_as_index,
                    min_coverage_for_variability=self.min_coverage_for_variability,
                    test_class=variability_test_class_null if self.report_variability_full else variability_test_class_default,
                    additional_per_position_data=additional_per_position_data,
                )
                nt_profile.process()
                split.SNV_profiles = nt_profile.d
                split.num_SNV_entries = len(nt_profile.d.get('coverage', []))

                # Process INDELs if not skipped
                if indel_acc is not None:
                    indels = indel_acc.finalize()
                    # Add pos_in_contig to each indel entry
                    for indel_hash, entry in indels.items():
                        entry['pos_in_contig'] = entry['pos'] + split.start

                    indel_profile = ProcessIndelCounts(
                        indels=indels,
                        coverage=allele_counts_array.sum(axis=0),
                        test_class=variability_test_class_null if self.report_variability_full else variability_test_class_default,
                        min_coverage_for_variability=self.min_coverage_for_variability if not self.report_variability_full else 1,
                    )
                    indel_profile.process()
                    split.INDEL_profiles = indel_profile.indels
                    split.num_INDEL_entries = len(split.INDEL_profiles)

                # Calculate variation density
                if split.num_SNV_entries > 0:
                    split.auxiliary = type('Auxiliary', (), {'variation_density': split.num_SNV_entries * 1000.0 / split.length})()
                else:
                    split.auxiliary = type('Auxiliary', (), {'variation_density': 0.0})()

                # Add redundant data for compatibility (same as process_contig)
                if split.num_SNV_entries > 0:
                    split.SNV_profiles['split_name'] = [split.name] * split.num_SNV_entries
                    split.SNV_profiles['sample_id'] = [self.sample_id] * split.num_SNV_entries
                    split.SNV_profiles['pos_in_contig'] = split.SNV_profiles['pos'] + split.start

                # SCV profiling: targeted second pass for codons containing SNVs
                # This reuses the existing run_SCVs logic which already filters to only
                # process genes that have SNVs, making the second pass efficient.
                if self.profile_SCVs and split.num_SNV_entries > 0:
                    # Create an Auxiliary object configured for SCV-only processing
                    aux = contigops.Auxiliary(
                        split,
                        profile_SCVs=True,
                        skip_INDEL_profiling=True,   # Already done in first pass
                        skip_SNV_profiling=True,     # Already done in first pass
                        min_coverage_for_variability=self.min_coverage_for_variability,
                        report_variability_full=self.report_variability_full,
                        min_percent_identity=self.min_percent_identity,
                        skip_edges=self.skip_edges
                    )

                    # Run SCV profiling (second pass, but only processes genes with SNVs)
                    aux.run_SCVs(bam_file)

                    # Add redundant data for SCV compatibility (same as process_contig)
                    for gene_id in split.SCV_profiles:
                        split.SCV_profiles[gene_id]['sample_id'] = [self.sample_id] * split.num_SCV_entries[gene_id]
                        split.SCV_profiles[gene_id]['corresponding_gene_call'] = [gene_id] * split.num_SCV_entries[gene_id]

        timer.make_checkpoint('Streaming analysis done')

        # Set up contig-level coverage (needed by downstream code)
        contig.coverage = bamops.Coverage()
        contig.coverage.c = contig_coverage
        contig.coverage.process_c(contig_coverage)

        # Clean up per_position_info (same as process_contig)
        del contig.coverage.c  # only split coverage array is needed
        for split in contig.splits:
            if hasattr(split, 'per_position_info'):
                del split.per_position_info

        if anvio.DEBUG:
            timer.gen_report('%s Time Report (streaming)' % contig.name)

        return contig


    def profile_single_thread(self):
        """The main method for anvi-profile when num_threads is 1"""

        bam_file = bamops.BAMFileObject(self.input_file_path)
        bam_file.fetch_filter = self.fetch_filter

        received_contigs = 0

        self.progress.new('Profiling w/1 thread', progress_total_items=self.num_contigs)

        mem_tracker = terminal.TrackMemory(at_most_every=5)
        mem_usage, mem_diff = mem_tracker.start()

        self.progress.update('contigs are being processed ...')
        for index in range(self.num_contigs):
            contig_name = self.contig_names[index]
            contig_length = self.contig_lengths[index]

            contig = self.process_contig(bam_file, contig_name, contig_length)

            self.contigs.append(contig)

            received_contigs += 1

            if mem_tracker.measure():
                mem_usage = mem_tracker.get_last()
                mem_diff = mem_tracker.get_last_diff()

            self.progress.increment(received_contigs)
            msg = '%d/%d contigs ⚙  | MEMORY 🧠  %s (%s)' % (received_contigs, self.num_contigs, mem_usage, mem_diff)
            self.progress.update(msg)

            # Here you're about to witness the poor side of Python (or our use of it). Although
            # we couldn't find any refs to these objects, garbage collecter kept them in the
            # memory. So here we are accessing to the atomic data structures in our split
            # objects to try to relieve the memory by encouraging the garbage collector to
            # realize what's up. Afterwards, we explicitly call the garbage collector
            if self.write_buffer_size > 0 and len(self.contigs) % self.write_buffer_size == 0:
                self.progress.update(f"{received_contigs}/{self.num_contigs} contigs ⚙ | WRITING TO DB 💾 ...")
                self.store_contigs_buffer()
                for c in self.contigs:
                    for split in c.splits:
                        del split.coverage
                        del split.auxiliary
                        del split
                    del c.splits[:]
                    del c.coverage
                    del c
                del self.contigs[:]
                gc.collect()

        self.progress.update(f"{received_contigs}/{self.num_contigs} contigs ⚙ | WRITING TO DB 💾 ...")
        self.store_contigs_buffer()
        self.auxiliary_db.close()

        self.progress.end(timing_filepath='anvio.debug.timing.txt' if anvio.DEBUG else None)

        overall_mean_coverage = 1
        if self.total_length_of_all_contigs != 0:
            overall_mean_coverage = self.total_coverage_values_for_all_contigs / self.total_length_of_all_contigs

        # FIXME: We know this is ugly. You can keep your opinion to yourself.
        if overall_mean_coverage > 0.0:
            # avoid dividing by zero
            dbops.ProfileDatabase(self.profile_db_path).db._exec("UPDATE abundance_splits SET value = value / " + str(overall_mean_coverage) + " * 1.0;")
            dbops.ProfileDatabase(self.profile_db_path).db._exec("UPDATE abundance_contigs SET value = value / " + str(overall_mean_coverage) + " * 1.0;")

        if not self.skip_SNV_profiling:
            self.layer_additional_data['num_SNVs_reported'] = self.variable_nts_table.num_entries
            self.layer_additional_keys.append('num_SNVs_reported')
            self.run.info("Num SNVs reported", self.layer_additional_data['num_SNVs_reported'], nl_before=1)

        if not self.skip_INDEL_profiling:
            self.layer_additional_data['num_INDELs_reported'] = self.indels_table.num_entries
            self.layer_additional_keys.append('num_INDELs_reported')
            self.run.info("Num INDELs reported", self.layer_additional_data['num_INDELs_reported'])

        if self.profile_SCVs:
            self.layer_additional_data['num_SCVs_reported'] = self.variable_codons_table.num_entries
            self.layer_additional_keys.append('num_SCVs_reported')
            self.run.info("Num SCVs reported", self.layer_additional_data['num_SCVs_reported'])

        if self.total_reads_kept != self.num_reads_mapped:
            # Num reads in profile do not equal num reads in bam
            diff = self.num_reads_mapped - self.total_reads_kept
            perc = (1 - self.total_reads_kept / self.num_reads_mapped) * 100

            if self.contig_names_of_interest:
                self.run.warning(f"There were {pp(diff)} reads present in the BAM file that did not end up being used "
                                 f"by the profiler, which corresponds to about {perc}% of all reads. This is "
                                 f"most likely due to the parameter `--contigs-of-interest`. Regardless, anvi'o "
                                 f"thought you should know about this.")
            elif self.fetch_filter:
                self.run.warning(f"There were {pp(diff)} reads present in the BAM file that did not end up being used "
                                 f"by the profiler, which corresponds to about {perc}% of all reads. This is "
                                 f"most likely due to your `--fetch-filter`, which, depending on the filter, can "
                                 f"make use of a very tiny fraction of all reads. If you think this is much more or "
                                 f"much less than what you would have expected, you may want to investigate it.")
            else:
                self.run.warning(f"There were {pp(diff)} reads present in the BAM file that did not end up being used "
                                 f"by the profiler, which corresponds to about {perc}% of all reads. Since you "
                                 f"don't seem to have used parameters such as `--contigs-of-interest` or `--fetch-filter` "
                                 f"anvi'o is Jon Snow here. There are other reasons why some of your reads may end up "
                                 f"not being considered by the profiler. For instance, if pysam encounteres reads that "
                                 f"it can't deal with (i.e., mapped reads in the BAM file with no defined sequences, or "
                                 f"read with defined sequences without mapping). Another reason can be that you have used "
                                 f"a new paramter in the profiler that removes contigs or reads from consideration. "
                                 f"Regardless, if you think this is worth your attention, now you know.")

        self.layer_additional_data['total_reads_kept'] = self.total_reads_kept
        self.layer_additional_keys.append('total_reads_kept')

        self.check_contigs(num_contigs=received_contigs)


    def profile_multi_thread(self):
        """The main method for anvi-profile when num_threads is >1"""

        manager = multiprocessing.Manager()
        available_index_queue = manager.Queue()
        output_queue = manager.Queue(self.queue_size)

        # Create shared memory for contig sequences and nt_positions_info to avoid
        # duplicating them across worker processes. This can save enormous amounts
        # of memory when using many threads.
        self.progress.new('Setting up shared memory')
        self.progress.update('Packing contig sequences into shared memory ...')
        # NOTE: We store the SharedSequenceStore as a local variable, NOT on self.
        # This is because SharedMemory objects don't pickle well, and self gets
        # pickled when passed to worker processes. We only store the name (string)
        # and index (dict) on self, which pickle correctly.
        shared_sequence_store = streamingops.SharedSequenceStore(self.contig_sequences)
        self._shared_sequence_shm_name = shared_sequence_store.shm_name
        self._shared_sequence_index = shared_sequence_store.index

        # Also create shared memory for nt_positions_info (per-nucleotide gene info)
        # This can be 8+ GB and would otherwise be duplicated across all workers
        self.progress.update('Packing nt_positions_info into shared memory ...')
        # Force load the lazy property before we can share it
        _ = self.nt_positions_info
        shared_nt_positions_store = streamingops.SharedNtPositionsStore(self.nt_positions_info)
        self._shared_nt_positions_shm_name = shared_nt_positions_store.shm_name
        self._shared_nt_positions_index = shared_nt_positions_store.index

        # Clear the original data to free memory in the main process
        # before forking workers. Workers will use shared memory instead.
        self.contig_sequences = None
        # Clear the cached nt_positions_info from LazyProperty cache.
        # LazyProperty stores cached values in instance._lazy_loaded_data[attr_name]
        if hasattr(self, '_lazy_loaded_data') and 'nt_positions_info' in self._lazy_loaded_data:
            del self._lazy_loaded_data['nt_positions_info']
        gc.collect()
        self.progress.end()

        # put contig indices into the queue to be read from within the worker
        for i in range(0, self.num_contigs):
            available_index_queue.put(i)

        processes = []
        for i in range(0, self.num_threads):
            processes.append(multiprocessing.Process(target=BAMProfiler.profile_contig_worker, args=(self, available_index_queue, output_queue)))

        for proc in processes:
            proc.start()

        received_contigs = 0

        self.progress.new('Profiling w/%d threads' % self.num_threads, progress_total_items=self.num_contigs)
        self.progress.update('initializing threads ...')

        mem_tracker = terminal.TrackMemory(at_most_every=5)
        mem_usage, mem_diff = mem_tracker.start()

        self.progress.update('contigs are being processed ...')
        while received_contigs < self.num_contigs:
            try:
                contig = output_queue.get()

                if isinstance(contig, Exception):
                    # If thread returns an exception, we raise it and kill the main thread.
                    raise contig

                self.contigs.append(contig)

                received_contigs += 1

                if mem_tracker.measure():
                    mem_usage = mem_tracker.get_last()
                    mem_diff = mem_tracker.get_last_diff()

                self.progress.increment(received_contigs)
                self.progress.update(f"{received_contigs}/{self.num_contigs} contigs ⚙ | MEMORY 🧠  {mem_usage} ({mem_diff}) ...")

                # Here you're about to witness the poor side of Python (or our use of it). Although
                # we couldn't find any refs to these objects, garbage collecter kept them in the
                # memory. So here we are accessing to the atomic data structures in our split
                # objects to try to relieve the memory by encouraging the garbage collector to
                # realize what's up. Afterwards, we explicitly call the garbage collector
                if self.write_buffer_size > 0 and len(self.contigs) % self.write_buffer_size == 0:
                    self.progress.update(f"{received_contigs}/{self.num_contigs} contigs ⚙ | WRITING TO DB 💾 ...")
                    self.store_contigs_buffer()
                    for c in self.contigs:
                        for split in c.splits:
                            del split.coverage
                            del split.auxiliary
                            del split
                        del c.splits[:]
                        del c.coverage
                        del c
                    del self.contigs[:]
                    gc.collect()

            except KeyboardInterrupt:
                self.run.info_single("Anvi'o profiler received SIGINT, terminating all processes...", nl_before=2)
                # Clean up shared memory before breaking
                shared_sequence_store.close()
                shared_sequence_store.unlink()
                shared_nt_positions_store.close()
                shared_nt_positions_store.unlink()
                break

            except Exception as worker_error:
                # An exception was thrown in one of the profile workers. We kill all processes in this case
                self.progress.end()
                for proc in processes:
                    proc.terminate()
                # Clean up shared memory before raising
                shared_sequence_store.close()
                shared_sequence_store.unlink()
                shared_nt_positions_store.close()
                shared_nt_positions_store.unlink()
                raise worker_error

        for proc in processes:
            proc.terminate()

        self.progress.update(f"{received_contigs}/{self.num_contigs} contigs ⚙ | WRITING TO DB 💾 ...")
        self.store_contigs_buffer()
        self.auxiliary_db.close()

        # Clean up shared memory
        shared_sequence_store.close()
        shared_sequence_store.unlink()
        shared_nt_positions_store.close()
        shared_nt_positions_store.unlink()

        self.progress.end(timing_filepath='anvio.debug.timing.txt' if anvio.DEBUG else None)

        overall_mean_coverage = 1
        if self.total_length_of_all_contigs != 0:
            overall_mean_coverage = self.total_coverage_values_for_all_contigs / self.total_length_of_all_contigs

        # FIXME: We know this is ugly. You can keep your opinion to yourself.
        if overall_mean_coverage > 0.0:
            # avoid dividing by zero
            dbops.ProfileDatabase(self.profile_db_path).db._exec("UPDATE abundance_splits SET value = value / " + str(overall_mean_coverage) + " * 1.0;")
            dbops.ProfileDatabase(self.profile_db_path).db._exec("UPDATE abundance_contigs SET value = value / " + str(overall_mean_coverage) + " * 1.0;")

        if not self.skip_SNV_profiling:
            self.layer_additional_data['num_SNVs_reported'] = self.variable_nts_table.num_entries
            self.layer_additional_keys.append('num_SNVs_reported')
            self.run.info("Num SNVs reported", self.layer_additional_data['num_SNVs_reported'], nl_before=1)

        if not self.skip_INDEL_profiling:
            self.layer_additional_data['num_INDELs_reported'] = self.indels_table.num_entries
            self.layer_additional_keys.append('num_INDELs_reported')
            self.run.info("Num INDELs reported", self.layer_additional_data['num_INDELs_reported'])

        if self.profile_SCVs:
            self.layer_additional_data['num_SCVs_reported'] = self.variable_codons_table.num_entries
            self.layer_additional_keys.append('num_SCVs_reported')
            self.run.info("Num SCVs reported", self.layer_additional_data['num_SCVs_reported'])

        if self.total_reads_kept != self.num_reads_mapped:
            # Num reads in profile do not equal num reads in bam
            diff = self.num_reads_mapped - self.total_reads_kept
            perc = (1 - self.total_reads_kept / self.num_reads_mapped) * 100
            self.run.warning("There were %d reads present in the BAM file that did not end up being used "
                             "by anvi'o. That corresponds to about %.2f percent of all reads in the bam file. "
                             "This could be either because you supplied --contigs-of-interest, "
                             "or because pysam encountered reads it could not deal with, e.g. they mapped "
                             "but had no defined sequence, or they had a sequence but did not map. "
                             "Regardless, anvi'o thought you should be aware of this." % (diff, perc))

        self.layer_additional_data['total_reads_kept'] = self.total_reads_kept
        self.layer_additional_keys.append('total_reads_kept')

        self.check_contigs(num_contigs=received_contigs)


    def store_contigs_buffer(self):
        for contig in self.contigs:
            self.total_length_of_all_contigs += contig.length
            self.total_coverage_values_for_all_contigs += contig.coverage.mean * contig.length

            # we will divide every abundance after profiling is done.
            contig.abundance = contig.coverage.mean
            for split in contig.splits:
                split.abundance = split.coverage.mean

            self.total_reads_kept += contig.coverage.num_reads

        self.generate_variable_nts_table()
        self.generate_variable_codons_table()
        self.generate_indels_table()
        self.store_split_coverages()

        # the crux of the profiling
        for atomic_data_field in self.essential_data_fields_for_anvio_profiles:
            view_data_splits, view_data_contigs = contigops.get_atomic_data(self.sample_id, self.contigs, atomic_data_field)

            table_name = '_'.join([atomic_data_field, 'splits'])
            TablesForViews(self.profile_db_path,
                           progress=self.progress) \
                                .create_new_view(view_data=view_data_splits,
                                                 table_name=table_name,
                                                 view_name=atomic_data_field,
                                                 append_mode=True)

            table_name = '_'.join([atomic_data_field, 'contigs'])
            TablesForViews(self.profile_db_path,
                           progress=self.progress) \
                                .create_new_view(view_data=view_data_contigs,
                                                 table_name=table_name,
                                                 view_name=None,
                                                 append_mode=True)


    def check_contigs(self, num_contigs=None):
        if not num_contigs:
            num_contigs = len(self.contigs)

        if not num_contigs:
            raise ConfigError("0 contigs to work with. Bye.")


    def cluster_contigs(self):
        default_clustering_config = constants.blank_default if self.blank else constants.single_default

        dbops.do_hierarchical_clustering_of_items(self.profile_db_path, self.clustering_configs, self.split_names, self.database_paths, \
                                                  input_directory=self.output_directory, default_clustering_config=default_clustering_config, \
                                                  distance=self.distance, linkage=self.linkage, run=self.run, progress=self.progress)


    def check_args(self):
        if self.blank:
            self.run.warning("You are about to generate a blank profile. This is what we do when we have nothing "
                             "but a contigs database to play with. Because anvi'o is lazy, it will not check the "
                             "rest of the parameters you may have declred. Most of them will not matter.")

            if not self.output_directory:
                raise ConfigError("If you want to generate a blank profile, you need to declare an output directory path.")
            if not self.sample_id:
                raise ConfigError("Mock profiles require a sample name to be declared. Because :/")
            return

        if (not self.input_file_path) and (not self.serialized_profile_path):
            raise ConfigError("You didn't declare any input files :/ If you intend to create a blank profile without any,\
                                input file, you should be a bit more explicit about your intention (you know, in the help\
                                there is a flag for it and all). Otherwise you should either provide an input BAM file, or\
                                a serialized anvi'o profile. See '--help' maybe?")
        if self.input_file_path and self.serialized_profile_path:
            raise ConfigError("You can't declare both an input file and a serialized profile.")
        if self.serialized_profile_path and (not self.output_directory):
            raise ConfigError("When loading serialized profiles, you need to declare an output directory.")
        if self.input_file_path and not os.path.exists(self.input_file_path):
            raise ConfigError("No such file: '%s'" % self.input_file_path)
        if self.serialized_profile_path and not os.path.exists(self.serialized_profile_path):
            raise ConfigError("No such file: '%s'" % self.serialized_profile_path)
        if not self.min_coverage_for_variability >= 0:
            raise ConfigError("Minimum coverage for variability must be 0 or larger.")
        if not self.min_contig_length >= 0:
            raise ConfigError("Minimum contig length must be 0 or larger.")
        if not self.max_contig_length >= 100:
            raise ConfigError("Maximum contig length can't be less than 100 base pairs.")
        if self.min_contig_length >= self.max_contig_length:
            raise ConfigError("Maximum contig length (%s) must be larger than the minimum "
                              "contig length (%s). Seriously though." % (pp(self.max_contig_length), pp(self.min_contig_length)))

        if self.num_threads < 1:
            raise ConfigError("Nice try. Obviously, number of threds can not be less than 1.")

        if not self.queue_size:
            self.queue_size = self.num_threads * 2

        if not self.write_buffer_size:
            self.run.warning("You set the write buffer size to 0. Which means, the profiling data will be kept in memory until "
                             "the very end of the processing.")

        if self.write_buffer_size < 0:
            raise ConfigError('No. Write buffer size can not have a negative value.')

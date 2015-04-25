# coding: utf-8
"""Summarizes information stored in a profile database using an annotation database for
   a given list of splits."""

import os
import sys
import numpy
import textwrap

from collections import Counter

import anvio.tables as t
import anvio.utils as utils
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths
import anvio.ccollections as ccollections
import anvio.completeness as completeness

from anvio.errors import ConfigError
from anvio.dbops import DatabasesMetaclass
from anvio.summaryhtml import SummaryHTMLOutput


__author__ = "A. Murat Eren"
__copyright__ = "Copyright 2015, The anvio Project"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = t.profile_db_version
__maintainer__ = "A. Murat Eren"
__email__ = "a.murat.eren@gmail.com"
__status__ = "Development"


pp = terminal.pretty_print
run = terminal.Run()
progress = terminal.Progress()


class Summarizer(DatabasesMetaclass):
    """Creates an über dictionary of 'summary'."""
    def __init__(self, args = None, r = run, p = progress):
        self.summary = {}

        self.debug = False
        self.profile_db_path = None
        self.annotation_db_path = None
        self.output_directory = None

        self.run = run
        self.progress = progress

        DatabasesMetaclass.__init__(self, args, self.run, self.progress)

        self.collections = ccollections.Collections()
        self.collections.populate_sources_dict(self.annotation_db_path, t.annotation_db_version)
        self.collections.populate_sources_dict(self.profile_db_path, t.profile_db_version)

        self.completeness = completeness.Completeness(self.annotation_db_path)

        self.collection_id = None

        if args:
            if args.list_collections:
                self.collections.list_collections()
                sys.exit()

            self.collection_id = args.collection_id
            self.output_directory = args.output_directory
            self.debug = args.debug

        self.sanity_check()

        filesnpaths.gen_output_directory(self.output_directory)

        # set samples
        self.samples = sorted([s.strip() for s in self.p_meta['samples'].split(',')])


    def sanity_check(self):
        if not self.collection_id:
            raise ConfigError, "You must specify a collection id :/"

        if self.collection_id not in self.collections.sources_dict:
            raise ConfigError, "%s is not a valid collection ID. See a list of available ones with '--list-collections' flag" % self.collection_id

        self.output_directory = filesnpaths.check_output_directory(self.output_directory)


    def process(self):
        # learn who you are:
        collection_dict = self.collections.get_collection_dict(self.collection_id)
        collection_colors = self.collections.get_collection_colors(self.collection_id)

        # set up the initial summary dictionary
        self.summary['meta'] = {'output_directory': self.output_directory,
                                'collection': collection_dict.keys(),
                                'num_bins': len(collection_dict.keys()),
                                'collection_id': self.collection_id,
                                'total_length': 0,
                                'num_contigs': 0,
                                'profile': self.p_meta,
                                'annotation': self.a_meta,
                                'samples': self.samples,
                                'percent_described': 0.0}

        self.summary['collection'] = {}

        for bin_id in collection_dict: 
            bin = Bin(self, bin_id, collection_dict[bin_id], self.run, self.progress)
            bin.output_directory = os.path.join(self.output_directory, 'collections', bin_id)

            self.summary['collection'][bin_id] = bin.create()
            self.summary['collection'][bin_id]['color'] = collection_colors[bin_id] or '#212121'
            self.summary['meta']['total_length'] += self.summary['collection'][bin_id]['total_length']
            self.summary['meta']['num_contigs'] += self.summary['collection'][bin_id]['num_contigs']

        # final additions
        self.summary['meta']['percent_described'] = '%.2f' % (self.summary['meta']['total_length'] * 100.0 / int(self.a_meta['total_length']))
        self.summary['meta']['bins'] = self.get_bins_ordered_by_completeness_and_size()

        SummaryHTMLOutput(self.summary, r = self.run, p = self.progress).generate()

        if self.debug:
            import json
            print json.dumps(self.summary, sort_keys=True, indent=4)


    def get_bins_ordered_by_completeness_and_size(self):
        return [t[2] for t in sorted([(self.summary['collection'][bin]['percent_complete'], self.summary['collection'][bin]['total_length'], bin) for bin in self.summary['collection']], reverse=True)]



class Bin:
    def __init__(self, summary, bin_id, split_ids, r = run, p = progress):
        self.summary = summary
        self.bin_id = bin_id
        self.split_ids = split_ids
        self.progress = p
        self.run = r

        self.bin_info_dict = {'files': {}}

        self.output_directory = None
        self.contig_lengths = []

        # make sure all split_ids in the collection is actually in the annotation database.
        # in collections stored in the annotation database, split_ids that are not in the
        # oritinal contigs used to generate annotation database *may* end up in the
        # collections table. we gotta make sure we deal with them properly:
        missing_ids = [split_id for split_id in self.split_ids if split_id not in self.summary.split_sequences]
        if len(missing_ids):
            for missing_id in missing_ids:
                self.split_ids.remove(missing_id)

            self.run.warning('%d split id(s) in bin "%s" reported by collection "%s" is not found in the\
                              annotation database and removed from the bin summary. If this does not make\
                              any sense, you may need make sure everything is in order. The thing is,\
                              sometimes external clustering results that are added to the annotation via\
                              `anvi-populate-collections-table` may include split names that are not used\
                              while the annotation database was generated.'\
                                                % (len(missing_ids), bin_id, self.summary.collection_id))


    def create(self):
        self.progress.new('[Collection "%s"] Creating the output directory' % self.bin_id)
        self.create_bin_dir()

        self.progress.new('[Collection "%s"] Creating the FASTA file' % self.bin_id)
        self.store_contigs_fasta()

        self.progress.new('[Collection "%s"] Accessing completeness scores' % self.bin_id)
        self.access_completeness_scores()

        self.progress.new('[Collection "%s"] Computing basic stats' % self.bin_id)
        self.compute_basic_stats()

        self.progress.new('[Collection "%s"] Filling in taxonomy info' % self.bin_id)
        self.set_taxon_calls()

        return self.bin_info_dict


    def create_bin_dir(self):
        self.progress.update('...')

        if not self.output_directory:
            self.progress.end()
            raise ConfigError, 'You caled Bin.create before setting an output directory. Anvio says "nope, thanks".'

        filesnpaths.gen_output_directory(self.output_directory)

        self.progress.end()


    def get_output_file_handle(self, prefix = 'output.txt', overwrite = False):
        file_path = os.path.join(self.output_directory, '%s-%s' % (self.bin_id, prefix))
        if os.path.exists(file_path) and not overwrite:
            raise ConfigError, 'get_output_file_handle: well, this file already exists: "%s"' % file_path


        key = prefix.split('.')[0].replace('-', '_')
        self.bin_info_dict['files'][key] = file_path[len(self.summary.output_directory):].strip('/')

        return open(file_path, 'w')


    def access_completeness_scores(self):
        self.progress.update('...')

        completeness = self.summary.completeness.get_info_for_splits(set(self.split_ids))

        self.bin_info_dict['completeness'] = completeness

        num_sources = len(completeness)

        # set up for the average completeness / contamination scores:
        for k in ['percent_contamination', 'percent_complete']:
            self.bin_info_dict[k] = 0.0

        # go through all single-copy gene reporting sources
        for c in completeness.values():
            for k in ['percent_contamination', 'percent_complete']:
                self.bin_info_dict[k] += c[k]

        for k in ['percent_contamination', 'percent_complete']:
            self.bin_info_dict[k] /= num_sources

        self.progress.end()


    def store_contigs_fasta(self):
        """Storing contig sequences.
        
           This is not an easy problem. We split contigs into smaller sequences at the beginning. Only
           a part of a given contig may be used during the binning process. On the other hand we can't
           simply store sequences of splits, whenever possible, we must store the entire sequence of
           the contig (only if all splits are selected from a contig in to the same bin). So, this
           function first identifies all splits coming from the same parent, then identifies sequential
           blocks of splits (see `SequentialBlocks` class), then checks whether all splits of a given
           contig is included in the bin. If that is the case, it puts the contig as a single entry,
           witht he identical FASTA id to the original contigs in the assembly file. Otherwise it appends
           `_partial_X_Y` to the FASTA id, X and Y being the start and stop positions.
        """

        fasta_file = self.get_output_file_handle('contigs.fa')

        # some null values:
        self.bin_info_dict['total_length'] = 0
        self.bin_info_dict['num_contigs'] = 0

        # this dict will keep all the contig ids found in this bin:
        contigs_represented = {}

        # go through all splits in this bin, and populate `contigs_represented`
        self.progress.update('Identifying contigs involved ...')
        for split_id in self.split_ids:
            s = self.summary.splits_basic_info[split_id]
            if s['parent'] in contigs_represented:
                contigs_represented[s['parent']][s['order_in_parent']] = split_id
            else:
                contigs_represented[s['parent']] = {s['order_in_parent']: split_id}

        # now it is time to go through each contig found in contigs_represented to
        # figure out how much of the contig is in fact in this bin
        for contig_id in contigs_represented:
            splits_order = contigs_represented[contig_id].keys()

            # this is critical: `sequential_blocks` is a list of one ore more lists,
            # each describes splits that follow each other to represent a coherent
            # chunk of the parent sequence:
            self.progress.update('Identifying sequential blocks ...')
            sequential_blocks = SequentialBlocks(splits_order).process()

            for sequential_block in sequential_blocks:
                self.progress.update('Identifying the portion of contig represented ...')
                first_split = contigs_represented[contig_id][sequential_block[0]]
                last_split = contigs_represented[contig_id][sequential_block[-1]]

                contig_sequence_start_in_splits = self.summary.splits_basic_info[first_split]['start']
                contig_sequence_end_in_splits = self.summary.splits_basic_info[last_split]['end']

                # so this much of the contig is represented by its splits:
                total_contig_length_in_splits = contig_sequence_end_in_splits - contig_sequence_start_in_splits

                # and this is is actual length:
                contig_sequence_length = self.summary.contigs_basic_info[contig_id]['length']

                if contig_sequence_length == total_contig_length_in_splits:
                    # the entireity of the contig is represented!
                    appendix = ''
                else:
                    appendix = '_partial_%d_%d' % (contig_sequence_start_in_splits, contig_sequence_end_in_splits)

                sequence = ''
                self.progress.update('Reconstructing contig sequence from splits ...')
                for split_order in sequential_block:
                    sequence += self.summary.split_sequences[contigs_represented[contig_id][split_order]]

                fasta_id = contig_id + appendix

                self.progress.update('Writing contig sequence into file ...')
                fasta_file.write('>%s\n' % fasta_id)
                fasta_file.write('%s\n' % textwrap.fill(sequence, 80, break_on_hyphens = False))

                # fill in basic info about contigs in bin
                len_seq = len(sequence)
                self.bin_info_dict['total_length'] += len_seq
                self.contig_lengths.append(len_seq)
                self.bin_info_dict['num_contigs'] += 1

        fasta_file.close()
        self.progress.end()


    def set_taxon_calls(self):
        self.progress.update('...')

        taxon_calls_counter = Counter()
        for split_id in self.split_ids:
            taxon_calls_counter[self.summary.genes_in_splits_summary_dict[split_id]['taxonomy']] += 1

        taxon_calls = sorted([list(tc) for tc in taxon_calls_counter.items()], key = lambda x: int(x[1]), reverse = True)

        self.bin_info_dict['taxon_calls'] = taxon_calls

        # taxon_calls = [(None, 129), ('Propionibacterium avidum', 120), ('Propionibacterium acnes', 5)]
        l = [tc for tc in taxon_calls if tc[0]]
        num_calls = sum(taxon_calls_counter.values())
        # l = [('Propionibacterium avidum', 120), ('Propionibacterium acnes', 5)]
        if l and l[0][1] > num_calls / 4.0:
            # if l[0] is associated with more than 25 percent of splits:
            self.bin_info_dict['taxon'] = l[0][0]
        else:
            self.bin_info_dict['taxon'] = 'Unknown'

        # convert to percents..
        for tc in taxon_calls:
            tc[1] = tc[1] * 100.0 / num_calls

        self.progress.end()


    def compute_basic_stats(self):
        self.progress.update('...')

        self.bin_info_dict['N50'] = utils.get_N50(self.contig_lengths)
        self.bin_info_dict['GC_content'] = numpy.mean([self.summary.splits_basic_info[split_id]['gc_content'] for split_id in self.split_ids]) * 100

        self.progress.end()


class SequentialBlocks:
    """Gets a list that goes like this: [1, 2, 3, 5, 6, 9], and returns another list
       that goes like this: [[1, 2, 3], [5, 6], [9]]"""
    def __init__(self, l):
        self.l = sorted(list(set(l)))
        self.blocks = []
        self.current_block = []


    def finalize_block(self):
        self.blocks.append(self.current_block)
        self.current_block = []


    def process(self):
        while 1:
            if not self.l:
                break

            current = self.l.pop(0)

            if not len(self.current_block) or current == self.current_block[-1] + 1:
                self.current_block.append(current)
            else:
                self.finalize_block()
                self.current_block.append(current)

        self.finalize_block()

        return self.blocks


#!/usr/bin/env python
# -*- coding: utf-8

# Copyright (C) 2014, A. Murat Eren
#
# This program is free software; you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free
# Software Foundation; either version 2 of the License, or (at your option)
# any later version.
#
# Please read the COPYING file.

import os
import sys
import json
import numpy
import pysam
import random
import cPickle
import operator
import subprocess
import PaPi.utils as utils
from PaPi.utils import pretty_print as pp
from PaPi.contig_stats_essential import Essential
from PaPi.contig_stats_auxiliary import Auxiliary
from PaPi.contig_composition import Composition


class BAMProfiler:
    """Creates an über class for BAM file operations"""
    def __init__(self, args = None):
        if args:
            self.args = args
            self.input_file_path = args.input_file
            self.serialized_profile_path = args.profile
            self.output_directory = args.output_directory
            self.list_contigs_and_exit = args.list_contigs
            self.min_contig_length = args.min_contig_length
            self.min_mean_coverage = args.min_mean_coverage
            self.number_of_threads = 4 
            self.no_trehading = True

            if args.contigs:
                if os.path.exists(args.contigs):
                    self.contigs_of_interest = [c.strip() for c in open(args.contigs).readlines() if c.strip() and not c.startswith('#')]
                else:
                    self.contigs_of_interest = [c.strip() for c in args.contigs.split(',')] if args.contigs else None

            else:
                self.contigs_of_interest = None

        else:
            self.args = args
            self.input_file_path = None 
            self.serialized_profile_path = None 
            self.output_directory = None 
            self.list_contigs_and_exit = None 
            self.min_contig_length = 10000 
            self.min_mean_coverage = 10
            # FIXME: Parameterize these two:
            self.number_of_threads = 4 
            self.no_trehading = False

        self.bam = None
        self.contigs_dict = {}

        self.progress = utils.Progress()
        self.run = utils.Run()


    def _run(self):
        self.check_args()

        if self.input_file_path:
            self.init_profile_from_BAM()
            self.profile()
            self.store_profile()
        else:
            self.init_serialized_profile()

        self.report()

        runinfo_serialized = self.generate_output_destination('RUNINFO.cPickle')
        self.run.info('runinfo', runinfo_serialized)
        self.run.store_info_dict(runinfo_serialized)
        self.run.quit()


    def init_serialized_profile(self):
        self.progress.new('Init')
        self.progress.update('Reading serialized profile')

        
        self.contigs_dict = cPickle.load(open(self.serialized_profile_path))
        self.progress.end()
        self.run.info('profile_loaded_from', self.serialized_profile_path)

        self.contigs = self.contigs_dict.keys()
        self.lengths = [self.contigs_dict[contig]['essential']['length'] for contig in self.contigs]

        self.run.info('num_contigs', pp(len(self.contigs)))

        if self.list_contigs_and_exit:
            print "\nContigs in the file:\n"
            for (contig, length) in zip(self.contigs, self.lengths):
                print "\t- %s (%s)" % (contig, pp(int(length)))
            print
            sys.exit()

        if self.contigs_of_interest:
            indexes = [self.contigs.index(r) for r in self.contigs_of_interest if r in self.contigs]
            self.contigs = [self.contigs[i] for i in indexes]
            self.lengths = [self.lengths[i] for i in indexes]
            self.run.info('num_contigs_selected_for_analysis', pp(len(self.contigs)))

        contigs_longer_than_M = set()
        for i in range(0, len(self.contigs)):
            if self.lengths[i] > self.min_contig_length:
                contigs_longer_than_M.add(i)
        if not len(contigs_longer_than_M):
            raise utils.ConfigError, "0 contigs larger than %s nts." % pp(self.min_contig_length)
        else:
            self.contigs = [self.contigs[i] for i in contigs_longer_than_M]
            self.contig_lenghts = [self.lengths[i] for i in contigs_longer_than_M]
            self.run.info('contigs_raw_longer_than_M', len(self.contigs))

        self.progress.new('Init')
        self.progress.update('Initializing the output directory ...')
        self.init_output_directory()
        self.progress.end()
        self.run.info('output_dir', self.output_directory)


    def init_profile_from_BAM(self):
        self.progress.new('Init')
        self.progress.update('Reading BAM File')
        self.bam = pysam.Samfile(self.input_file_path, 'rb')
        self.progress.end()
        self.run.info('input_bam', self.input_file_path)

        self.contigs = self.bam.references
        self.contig_lenghts = self.bam.lengths

        try:
            self.num_reads_mapped = self.bam.mapped
        except ValueError:
            raise utils.ConfigError, "It seems the BAM file is not indexed. See 'papi-init-bam' script."

        self.progress.new('Init')
        self.progress.update('Initializing the output directory ...')
        self.init_output_directory()
        self.progress.end()

        runinfo = self.generate_output_destination('RUNINFO')
        self.run.init_info_file_obj(runinfo)
        self.run.info('output_dir', self.output_directory)

        self.run.info('total_reads_mapped', pp(int(self.num_reads_mapped)))
        self.run.info('num_contigs', pp(len(self.contigs)))

        if self.list_contigs_and_exit:
            print "\nContigs in the file:\n"
            for (contig, length) in zip(self.contigs, self.contig_lenghts):
                print "\t- %s (%s)" % (contig, pp(int(length)))
            print
            sys.exit()

        if self.contigs_of_interest:
            indexes = [self.contigs.index(r) for r in self.contigs_of_interest if r in self.contigs]
            self.contigs = [self.contigs[i] for i in indexes]
            self.contig_lenghts = [self.contig_lenghts[i] for i in indexes]
            self.run.info('num_contigs_selected_for_analysis', pp(len(self.contigs)))

        contigs_longer_than_M = set()
        for i in range(0, len(self.contigs)):
            if self.contig_lenghts[i] > self.min_contig_length:
                contigs_longer_than_M.add(i)
        if not len(contigs_longer_than_M):
            raise utils.ConfigError, "0 contigs larger than %s nts." % pp(self.min_contig_length)
        else:
            self.contigs = [self.contigs[i] for i in contigs_longer_than_M]
            self.contig_lenghts = [self.contig_lenghts[i] for i in contigs_longer_than_M]
            self.run.info('contigs_raw_longer_than_M', len(self.contigs))


    def init_output_directory(self):
        Absolute = lambda x: os.path.join(os.getcwd(), x) if not x.startswith('/') else x

        if not self.output_directory:
            self.output_directory = Absolute(self.input_file_path) + '-PaPi-OUTPUT'
        else:
            self.output_directory = Absolute(self.output_directory)

        if not os.path.exists(self.output_directory):
            try:
                os.makedirs(self.output_directory)
            except:
                self.progress.end()
                raise utils.ConfigError, "Output directory does not exist (attempt to create one failed as well): '%s'" % \
                                                                (self.output_directory)
        if not os.access(self.output_directory, os.W_OK):
            self.progress.end()
            raise utils.ConfigError, "You do not have write permission for the output directory: '%s'" % self.output_directory


    def generate_output_destination(self, postfix, directory = False):
        return_path = os.path.join(self.output_directory, postfix)

        if directory == True:
            if os.path.exists(return_path):
                shutil.rmtree(return_path)
            os.makedirs(return_path)

        return return_path


    def profile(self):
        """Big deal function"""

        self.progress.new('Profiling the BAM file for Essential Stats')
        # So we start with essential stats. In the section below, we will simply go through each contig (contig),
        # in the BAM file and populate the contigs dictionary for the first time. There are two major sections,
        # one for no_threading option, and the other with multiple threads.
        if self.no_trehading:
            for i in range(0, len(self.contigs)):
                contig = self.contigs[i]
                self.contigs_dict[contig] = {}

                #fill in basics
                self.progress.update('Essential stats for "%s" (%d of %d) ...' % (contig, i + 1, len(self.contigs)))
                self.contigs_dict[contig]['essential'] = Essential(contig, self.bam.pileup(contig)).report()
                self.contigs_dict[contig]['essential']['length'] = self.contig_lenghts[i]

        else:
            def worker(contig, shared_contigs_dict):
                shared_contigs_dict[contig] = Essential(contig, self.bam.pileup(contig)).report()

            mp = utils.Multiprocessing(worker, self.number_of_threads)
            shared_contigs_dict = mp.get_empty_shared_dict()

            # arrange processes
            processes_to_run = []
            for contig in self.contigs:
                processes_to_run.append((contig, shared_contigs_dict),)

            # start the main loop to run all processes
            mp.run_processes(processes_to_run, self.progress)
            for contig in self.contigs:
                self.contigs_dict[contig] = {}
                self.contigs_dict[contig]['essential'] = shared_contigs_dict.pop(contig)
        self.progress.end()


        # filtering based on C starts.
        self.progress.new('Filtering contigs based on mean coverage')
        # this is also important. here we are going to remove any contig with a mean coverage less than C; mean
        # coverage info is stored in ['essential']['mean_coverage']. the mean coverage does not include areas
        # where zero reads mapped. 
        self.progress.update('Screening coverage for each contig ...')
        contigs_to_discard = set()
        for contig in self.contigs_dict:
            if self.contigs_dict[contig]['essential']['mean_coverage'] < self.min_mean_coverage:
                contigs_to_discard.add(contig)

        if len(contigs_to_discard):
            for contig in contigs_to_discard:
                self.contigs_dict.pop(contig)
            self.contigs = self.contigs_dict.keys()
            self.progress.end()
            self.run.info('contigs_after_C', pp(len(self.contigs)))
        else:
            self.progress.end()

        if not len(self.contigs):
            raise utils.ConfigError, "0 contigs passed minimum mean coverage parameter."

        # QA/QC is done. Now we go into Auxiliary analyses.
        self.progress.new('Computing auxiliary stats')
        if self.no_trehading:
            for i in range(0, len(self.contigs)):
                contig = self.contigs[i]
                # fill in entropy and representatives
                self.progress.update('Auxiliary stats for "%s" (%d of %d) ...' % (contig, i + 1, len(self.contigs)))
                self.contigs_dict[contig]['auxiliary'] = Auxiliary(contig,
                                                                         self.bam.pileup(contig),
                                                                         self.contigs_dict[contig]['essential']\
                                                                         ).report()
        else:
            def worker(contig, shared_contigs_dict):
                shared_contigs_dict[contig] = Auxiliary(contig,
                                                              self.bam.pileup(contig),
                                                              self.contigs_dict[contig]['essential']\
                                                              ).report()

            mp = utils.Multiprocessing(worker, self.number_of_threads)
            shared_contigs_dict = mp.get_empty_shared_dict()

            # arrange processes
            processes_to_run = []
            for contig in self.contigs:
                processes_to_run.append((contig, shared_contigs_dict),)

            # start the main loop to run all processes
            mp.run_processes(processes_to_run, self.progress)
            for contig in self.contigs:
                self.contigs_dict[contig]['auxiliary'] = shared_contigs_dict.pop(contig)
        self.progress.end()


        # it is time to fill in the tetranucleotide frequency info per contig
        self.progress.new('TNF Stats')
        if self.no_trehading:
            for i in range(0, len(self.contigs)):
                contig = self.contigs[i]
                self.progress.update('Computing TNF for "%s" (%d of %d) ...' % (contig, i + 1, len(self.contigs)))
                self.contigs_dict[contig]['composition'] = Composition(contig,
                                                                             self.contigs_dict[contig]['auxiliary']\
                                                                             ).report()
        else:
            def worker(contig, shared_contigs_dict):
                shared_contigs_dict[contig] = Composition(contig,
                                                                self.contigs_dict[contig]['auxiliary']\
                                                                ).report()

            mp = utils.Multiprocessing(worker, self.number_of_threads)
            shared_contigs_dict = mp.get_empty_shared_dict()

            # arrange processes
            processes_to_run = []
            for contig in self.contigs:
                processes_to_run.append((contig, shared_contigs_dict),)

            # start the main loop to run all processes
            mp.run_processes(processes_to_run, self.progress)
            for contig in self.contigs:
                self.contigs_dict[contig]['composition'] = shared_contigs_dict.pop(contig)
        self.progress.end()

        # Profiling is done.


    def store_profile(self):
        output_file = self.generate_output_destination('PROFILE.cPickle')
        self.progress.new('Storing Profile')
        self.progress.update('Serializing information for %s contigs ...' % pp(len(self.contigs_dict)))
        cPickle.dump(self.contigs_dict, open(output_file, 'w'))
        self.progress.end()
        self.run.info('profile_dict', output_file)


    def load_profile(self):
        pass


    def report(self):
        # generate a sorted list of contigs based on length
        self.contigs = [t[1] for t in sorted([(self.contigs_dict[k]['essential']['length'], k)\
                                                for k in self.contigs], reverse = True)]

        self.progress.new('Generating reports')
        self.progress.update('TNF matrix for contigs')
        TNF_matrix_file_path = self.generate_output_destination('TETRANUCLEOTIDE-FREQ-MATRIX.txt')
        output = open(TNF_matrix_file_path, 'w')
        kmers = sorted(self.contigs_dict[self.contigs[0]]['composition']['tnf'].keys())
        output.write('contigs\t%s\n' % ('\t'.join(kmers)))
        for contig in self.contigs:
            output.write('%s\t' % (contig))
            output.write('%s\n' % '\t'.join([str(self.contigs_dict[contig]['composition']['tnf'][kmer]) for kmer in kmers]))
        output.close()
        self.progress.end()
        self.run.info('tnf_matrix', TNF_matrix_file_path)


        self.progress.new('Generating reports')
        self.progress.update('Generating the tree of contigs')
        newick_tree_file_path = self.generate_output_destination('TNF-NEWICK-TREE.txt')
        env = os.environ.copy()
        subprocess.call(['papi-TNF-matrix-to-newick.R', '-o', newick_tree_file_path, TNF_matrix_file_path], env = env)
        self.progress.end()
        self.run.info('tnf_tree', newick_tree_file_path)


        # metadata
        self.progress.new('Generating reports')
        self.progress.update('Metadata for contigs')
        metadata_fields = [('essential', 'length'), ('essential', 'mean_coverage'), ('essential', 'std_coverage'), 
                           ('composition', 'GC_content')]

        metadata_txt = open(self.generate_output_destination('METADATA.txt'), 'w')
        metadata_json = open(self.generate_output_destination('METADATA.json'), 'w')
        metadata_json_buffer = []

        fields = [m[1] for m in metadata_fields]
        metadata_txt.write('contigs\t%s\n' % ('\t'.join(fields)))
        metadata_json_buffer.append([''] + fields)

        for contig in self.contigs:
            l = [self.contigs_dict[contig][major][minor] for major, minor in metadata_fields]
            fields = [contig] + ['%.4f' % x for x in l]
            metadata_txt.write('%s\n' % '\t'.join(fields))
            metadata_json_buffer.append(fields)

        metadata_txt.close()
        metadata_json.write(json.dumps(metadata_json_buffer))
        metadata_json.close()
        self.progress.end()
        self.run.info('metadata_txt', metadata_txt.name)
        self.run.info('metadata_json', metadata_json.name)


        # contigs FASTA
        self.progress.new('Generating reports')
        self.progress.update('Consensus FASTA file for contigs')
        contigs_fasta = open(self.generate_output_destination('CONTIGS-CONSENSUS.fa'), 'w')
        for contig in self.contigs:
            contigs_fasta.write(">%s\n%s\n" % (contig,
                                               self.contigs_dict[contig]['auxiliary']['rep_seq']))
        contigs_fasta.close()
        self.progress.end()
        self.run.info('contigs_fasta', contigs_fasta.name)


    def check_args(self):
        if (not self.input_file_path) and (not self.serialized_profile_path):
            raise utils.ConfigError, "You must declare either an input file, or a serialized profile."
        if self.input_file_path and self.serialized_profile_path:
            raise utils.ConfigError, "You can't declare both an input file and a serialized profile."
        if self.serialized_profile_path and (not self.output_directory):
            raise utils.ConfigError, "When loading serialized profiles, you need to declare an output directory."
        if self.input_file_path and not os.path.exists(self.input_file_path):
            raise utils.ConfigError, "No such file: '%s'" % self.input_file_path
        if self.serialized_profile_path and not os.path.exists(self.serialized_profile_path):
            raise utils.ConfigError, "No such file: '%s'" % self.serialized_profile_path
        if not self.min_mean_coverage > 0:
            raise utils.ConfigError, "Minimum mean coverage must be 1 or larger"
        if not self.min_contig_length > 0:
            raise utils.ConfigError, "Minimum contig length must be 1 or larger (although using anything below 5,000 is kinda silly)."

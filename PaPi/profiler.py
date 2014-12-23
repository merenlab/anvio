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
import numpy
import pysam
import random
import string
import shutil
import operator
import subprocess

import PaPi.db as db
import PaPi.utils as utils
import PaPi.dictio as dictio
import PaPi.terminal as terminal
import PaPi.constants as constants
import PaPi.annotation as annotation
import PaPi.filesnpaths as filesnpaths
import PaPi.clustering as clustering

from PaPi.genes import Genes
from PaPi.contig import Split
from PaPi.contig import Contig, set_contigs_abundance
from PaPi.clusteringconfuguration import ClusteringConfiguration

pp = terminal.pretty_print

__version__ = '0.4.1'


class BAMProfiler:
    """Creates an über class for BAM file operations"""
    def __init__(self, args = None):
        self.args = None
        self.input_file_path = None 
        self.annotation_db_path = None
        self.serialized_profile_path = None 
        self.output_directory = None 
        self.list_contigs_and_exit = None 
        self.min_contig_length = 10000 
        self.min_mean_coverage = 10
        self.split_length = 20000
        self.contig_names_of_interest = None

        if args:
            self.args = args
            self.input_file_path = args.input_file
            self.annotation_db_path = args.annotation_db_path
            self.serialized_profile_path = args.profile
            self.output_directory = args.output_directory
            self.list_contigs_and_exit = args.list_contigs
            self.min_contig_length = args.min_contig_length
            self.min_mean_coverage = args.min_mean_coverage
            self.number_of_threads = 4 
            self.no_trehading = True
            self.split_length = args.split_length
            self.sample_id = args.sample_id

            if args.contigs_of_interest:
                if not os.path.exists(args.contigs_of_interest):
                    raise utils.ConfigError, "Contigs file (%s) is missing..." % (args.contigs_of_interest)

                self.contig_names_of_interest = set([c.strip() for c in open(args.contigs_of_interest).readlines()\
                                                                               if c.strip() and not c.startswith('#')])

        self.bam = None
        self.contigs = {}
        self.contig_ORFs = {}
        self.annotation_db = None

        self.profile_db = None
        self.profile_db_path = None

        self.clustering_configs = constants.clustering_configs['single']

        self.progress = terminal.Progress()
        self.run = terminal.Run(width=35)


    def init_dirs_and_dbs(self):
        self.progress.new('Initializing')

        self.progress.update('Creating the output directory ...')

        Absolute = lambda x: os.path.join(os.getcwd(), x) if not x.startswith('/') else x

        if not self.output_directory:
            self.output_directory = Absolute(self.input_file_path) + '-PAPI_PROFILE'
        else:
            self.output_directory = Absolute(self.output_directory)
        filesnpaths.gen_output_directory(self.output_directory, self.progress)

        self.progress.update('Initializing the profile database ...')

        # init a new db
        self.profile_db_path = self.generate_output_destination('PROFILE.db')
        self.profile_db = db.DB(self.profile_db_path, __version__, new_database = True)

        # put sample id into the meta table
        self.profile_db.set_meta_value('sample_id', self.sample_id)
        self.profile_db.set_meta_value('merged', False)

        if self.annotation_db_path:
            self.progress.update('Initializing the annotation database ...')
            self.annotation_db = annotation.Annotation(self.annotation_db_path)
            self.annotation_db.init_database()

        self.progress.end()


    def _run(self):
        self.check_args()

        self.set_sample_id()

        if not self.list_contigs_and_exit:
            self.init_dirs_and_dbs()

        # we will set up things here so the information in the annotation_db
        # can be utilized directly from within the contigs for loop. contig to
        # gene associations will be stored in self.contig_ORFs dictionary for
        # fast access.
        if self.annotation_db and not self.list_contigs_and_exit:
            self.populate_contig_ORFs()

        self.run.info('profiler_version', __version__)
        self.run.info('sample_id', self.sample_id)
        self.run.info('profile_db', self.profile_db_path)
        self.run.info('annotation_db', True if self.annotation_db_path else False, quiet = True)
        self.run.info('default_clustering', constants.single_default)
        self.run.info('cmd_line', utils.get_cmd_line())
        self.run.info('merged', False)

        # this is kinda important. we do not run full-blown profile function if we are dealing with a summarized
        # profile...
        if self.input_file_path:
            self.init_profile_from_BAM()
            self.profile()
            self.store_profile()
            self.store_summarized_profile()
        else:
            self.init_serialized_profile()
            self.store_summarized_profile()

        if self.annotation_db:
            # generate and store genes table
            self.generate_genes_db()
            self.run.info('genes_table', True, quiet = True)
        else:
            self.run.info('genes_table', False, quiet = True)

        self.report()
        self.cluster_contigs()

        runinfo_serialized = self.generate_output_destination('RUNINFO.cp')
        self.run.info('runinfo', runinfo_serialized)
        self.run.store_info_dict(runinfo_serialized, strip_prefix = self.output_directory)

        self.run.quit()
        self.profile_db.disconnect()


    def populate_contig_ORFs(self):
        self.progress.new('Annotation')
        self.progress.update('Reading annotation table')
        annotation_table = self.annotation_db.db.get_table_as_dict('annotation', annotation.annotation_table_structure)

        self.progress.update('Populating ORFs dictionary for each contig ...')
        for gene in annotation_table:
            e = annotation_table[gene]
            if self.contig_ORFs.has_key(e['contig']):
                self.contig_ORFs[e['contig']].add((gene, e['start'], e['stop']), )
            else:
                self.contig_ORFs[e['contig']] = set([(gene, e['start'], e['stop']), ])

        self.progress.end()
        self.run.info('annotation_db', "%d genes processed successfully." % len(annotation_table), display_only = True)


    def generate_genes_db(self):
        self.genes = Genes(self.progress)
        self.progress.new('Profiling genes')
        num_contigs = len(self.contigs)
        contig_names = self.contigs.keys()
        for i in range(0, num_contigs):
            contig = contig_names[i]
            self.progress.update('Processing contig %d of %d' % (i + 1, num_contigs))
            self.genes.analyze_contig(self.contigs[contig], self.sample_id, self.contig_ORFs[contig])

        self.genes.create_genes_table(self.profile_db)
        self.progress.end()
        


    def set_sample_id(self):
        if self.sample_id:
            utils.check_sample_id(self.sample_id)
        else:
            if self.input_file_path:
                self.input_file_path = utils.ABS(self.input_file_path)
                self.sample_id = os.path.basename(self.input_file_path).upper().split('.BAM')[0]
                utils.check_sample_id(self.sample_id)
            if self.serialized_profile_path:
                self.serialized_profile_path = utils.ABS(self.serialized_profile_path)
                self.sample_id = os.path.basename(os.path.dirname(self.serialized_profile_path))


    def init_serialized_profile(self):
        self.progress.new('Init')
        self.progress.update('Reading serialized profile')
        self.contigs = dictio.read_serialized_object(self.serialized_profile_path)
        self.progress.end()

        self.run.info('profile_loaded_from', self.serialized_profile_path)
        self.run.info('num_contigs', pp(len(self.contigs)))

        # we learn split length from the structure itself (and here we poorly assume that
        # all contigs have the same split length). so command line argument of split
        # length does not mean anything when the profile is loaded from a serialized
        # file.
        self.split_length = self.contigs.values()[0].split_length

        if self.list_contigs_and_exit:
            for tpl in sorted([(int(self.contigs[contig].length), contig) for contig in self.contigs]):
                print '%-40s %s' % (tpl[1], pp(int(tpl[0])))
            sys.exit()

        if self.contig_names_of_interest:
            contigs_to_discard = set()
            for contig in self.contigs:
                if contig not in self.contig_names_of_interest:
                    contigs_to_discard.add(contig)

            if len(contigs_to_discard):
                for contig in contigs_to_discard:
                    self.contigs.pop(contig)
            self.run.info('num_contigs_selected_for_analysis', pp(len(self.contigs)))

        self.check_contigs()

        # lets make sure that each contig name is found in the annotation db
        if self.annotation_db:
            for contig in self.contigs:
                if contig not in self.contig_ORFs:
                    raise utils.ConfigError, "You instructed the profiling to use an annotation database,\
                                              however, at least one contig ('%s') PaPi found in your input file that\
                                              is missing from the database. This would have never happened if you had\
                                              generated the annotation database and done mapping using the same contigs.\
                                              At this point PaPi does not know who to blame, but this is\
                                              not right :/" % contig

        contigs_to_discard = set()
        for contig in self.contigs.values():
            if contig.length < self.min_contig_length:
                contigs_to_discard.add(contig.name)
        if len(contigs_to_discard):
            for contig in contigs_to_discard:
                self.contigs.pop(contig)
            self.run.info('contigs_raw_longer_than_M', len(self.contigs))

        self.check_contigs()


    def init_profile_from_BAM(self):
        self.progress.new('Init')
        self.progress.update('Reading BAM File')
        self.bam = pysam.Samfile(self.input_file_path, 'rb')
        self.progress.end()

        self.contig_names = self.bam.references
        self.contig_lenghts = self.bam.lengths

        if self.list_contigs_and_exit:
            for tpl in sorted(zip(self.contig_lenghts, self.contig_names), reverse = True):
                print '%-40s %s' % (tpl[1], pp(int(tpl[0])))
            sys.exit()

        utils.check_contig_names(self.contig_names)

        try:
            self.num_reads_mapped = self.bam.mapped
        except ValueError:
            raise utils.ConfigError, "It seems the BAM file is not indexed. See 'papi-init-bam' script."

        runinfo = self.generate_output_destination('RUNINFO')
        self.run.init_info_file_obj(runinfo)
        self.run.info('input_bam', self.input_file_path)
        self.run.info('output_dir', self.output_directory)
        self.run.info('total_reads_mapped', pp(int(self.num_reads_mapped)))
        self.run.info('num_contigs', pp(len(self.contig_names)))

        if self.contig_names_of_interest:
            indexes = [self.contig_names.index(r) for r in self.contig_names_of_interest if r in self.contig_names]
            self.contig_names = [self.contig_names[i] for i in indexes]
            self.contig_lenghts = [self.contig_lenghts[i] for i in indexes]
            self.run.info('num_contigs_selected_for_analysis', pp(len(self.contig_names)))


        # lets make sure that each contig name is found in the annotation db
        if self.annotation_db:
            for contig in self.contig_names:
                if contig not in self.contig_ORFs:
                    raise utils.ConfigError, "You instructed the profiling to use an annotation database,\
                                              however, at least one contig ('%s') PaPi found in the BAM file is\
                                              missing from the database.  This would have never happened if you had\
                                              generated the annotation database and done mapping using the same contigs.\
                                              At this point PaPi does not know who to blame, but this is\
                                              not right :/" % contig

        # check for the -M parameter.
        contigs_longer_than_M = set()
        for i in range(0, len(self.contig_names)):
            if self.contig_lenghts[i] > self.min_contig_length:
                contigs_longer_than_M.add(i)
        if not len(contigs_longer_than_M):
            raise utils.ConfigError, "0 contigs larger than %s nts." % pp(self.min_contig_length)
        else:
            self.contig_names = [self.contig_names[i] for i in contigs_longer_than_M]
            self.contig_lenghts = [self.contig_lenghts[i] for i in contigs_longer_than_M]
            self.run.info('contigs_raw_longer_than_M', len(self.contig_names))

        # finally, compute contig splits.
        self.contig_splits = [utils.get_chunks(self.contig_lenghts[i], self.split_length)\
                                                                 for i in range(0, len(self.contig_names))]


    def generate_output_destination(self, postfix, directory = False):
        return_path = os.path.join(self.output_directory, postfix)

        if directory == True:
            if os.path.exists(return_path):
                shutil.rmtree(return_path)
            os.makedirs(return_path)

        return return_path


    def profile(self):
        """Big deal function"""

        # So we start with essential stats. In the section below, we will simply go through each contig
        # in the BAM file and populate the contigs dictionary for the first time.
        for i in range(0, len(self.contig_names)):
        
            contig_name = self.contig_names[i]
            contig_splits = self.contig_splits[i]

            contig = Contig(contig_name)
            contig.length = self.contig_lenghts[i]
            contig.split_length = self.split_length

            self.progress.new('Profiling "%s" (%d of %d) (%s nts)' % (contig.name,
                                                                      i + 1,
                                                                      len(self.contig_names),
                                                                      pp(int(contig.length))))

            # populate contig with empty split objects and 
            for split_order in range(0, len(contig_splits)):
                start, end = contig_splits[split_order]
                split = Split(contig.name, split_order, start, end)
                contig.splits.append(split)

            # analyze coverage for each split
            contig.analyze_coverage(self.bam, self.progress)

            # test the mean coverage of the contig.
            discarded_contigs_due_to_C = set([])
            if contig.coverage.mean < self.min_mean_coverage:
                # discard this contig and continue
                discarded_contigs_due_to_C.add(contig.name)
                self.progress.end()
                continue

            contig.analyze_auxiliary(self.bam, self.progress)

            contig.analyze_composition(self.bam, self.progress)

            contig.analyze_tnf(self.progress)

            self.progress.end()

            # add contig to the dict.
            self.contigs[contig_name] = contig


        if discarded_contigs_due_to_C:
            self.run.info('contigs_after_C', pp(len(self.contigs)))

        # set contig abundance
        set_contigs_abundance(self.contigs)

        self.check_contigs()


    def store_profile(self):
        output_file = self.generate_output_destination('PROFILE.cp')
        self.progress.new('Storing Profile')
        self.progress.update('Serializing information for %s contigs ...' % pp(len(self.contigs)))
        dictio.write_serialized_object(self.contigs, output_file)
        self.progress.end()
        self.run.info('profile_dict', output_file)


    def store_summarized_profile(self):
        summary_index = {}
        summary_index_output_path = self.generate_output_destination('SUMMARY.cp')
        summary_dir = self.generate_output_destination('SUMMARY', directory=True)
        self.progress.new('Storing summary files')

        counter = 1

        for contig in self.contigs:
            self.progress.update('working on contig %s of %s' % (pp(counter), pp(len(self.contigs))))
            for split in self.contigs[contig].splits:
                split_summary_path = self.generate_output_destination(os.path.join(summary_dir, '%.6d.cp' % counter))
                dictio.write_serialized_object({self.sample_id: {'coverage': split.coverage.c,
                                                                 'variability': split.auxiliary.v,
                                                                 'competing_nucleotides': split.auxiliary.competing_nucleotides}},
                                                                 split_summary_path)
                summary_index[split.name] = split_summary_path
                counter += 1

        self.progress.end()
        self.run.info('profile_summary_dir', summary_dir)
        dictio.write_serialized_object(dictio.strip_prefix_from_dict_values(summary_index, self.output_directory), summary_index_output_path)
        self.run.info('profile_summary_index', summary_index_output_path)


    def check_contigs(self):
        if not len(self.contigs):
            raise utils.ConfigError, "0 contigs to work with. Bye."


    def report(self):
        self.run.info('split_length', self.split_length)
        self.run.info('min_contig_length', self.min_contig_length)
        self.run.info('min_mean_coverage', self.min_mean_coverage)

        # generate a sorted list of contigs based on length
        self.contig_names = [t[1] for t in sorted([(self.contigs[k].length, k)\
                                                for k in self.contigs], reverse = True)]

        for target in ["contigs", "splits"]:
            self.progress.new('TNF Matrix')
            self.progress.update('Being generated for %s' % target)
            matrix_path = self.generate_output_destination('TNF-MATRIX-%s.txt' % target.upper())
            output = open(matrix_path, 'w')
            kmers = sorted(self.contigs[self.contigs.keys()[0]].tnf.keys())
            output.write('contigs\t%s\n' % ('\t'.join(kmers)))
            for contig in self.contigs:
                for split in self.contigs[contig].splits:
                    output.write('%s\t' % (split.name))
                    if target == "contigs":
                        output.write('%s\n' % '\t'.join([str(self.contigs[contig].tnf[kmer]) for kmer in kmers]))
                    else:
                        output.write('%s\n' % '\t'.join([str(split.tnf[kmer]) for kmer in kmers]))
            output.close()
            self.progress.end()
            self.run.info('tnf_matrix_%s' % target, matrix_path)

        F = lambda x: '%.4f' % x
        I = lambda x: '%d' % x

        # metadata
        for target in ["contigs", "splits"]:
            self.progress.new('Metadata')
            self.progress.update('Being generated for %s' % target)
            metadata_fields_to_report = ['contigs', 'length', 'GC_content', 'std_coverage', 'mean_coverage', 'noramlized_coverage', 'portion_covered', 'abundance', 'variability', '__parent__']
            metadata_txt = open(self.generate_output_destination('METADATA-%s.txt' % target.upper()), 'w')
            metadata_txt.write('%s\n' % ('\t'.join(metadata_fields_to_report)))
            for contig_name in self.contigs:
                contig = self.contigs[contig_name]
                for split in contig.splits:
                    if target == "contigs":
                        obj = contig
                    else:
                        obj = split

                    # FIXME: except for the split.name, all should be obj's. but the poor design in contig.py
                    # prevents that. 
                    fields = [split.name,
                              I(split.length),
                              F(split.composition.GC_content),
                              F(obj.coverage.std),
                              F(obj.coverage.mean),
                              F(obj.coverage.normalized),
                              F(obj.coverage.portion_covered),
                              F(obj.abundance),
                              F(split.auxiliary.variability_score),
                              contig_name] 
                    metadata_txt.write('%s\n' % '\t'.join(fields))
            metadata_txt.close()
            self.progress.end()
            self.run.info('metadata_%s' % target, metadata_txt.name)


        # splits FASTA
        self.progress.new('Consensus FASTA')
        splits_fasta = open(self.generate_output_destination('SPLITS-CONSENSUS.fa'), 'w')
        for contig in self.contig_names:
            for split in self.contigs[contig].splits:
                splits_fasta.write(">%s\n%s\n" % (split.name,
                                                  split.auxiliary.rep_seq))
        splits_fasta.close()
        self.progress.end()
        self.run.info('splits_fasta', splits_fasta.name)


    def cluster_contigs(self):
        clusterings = {}
        for config_name in self.clustering_configs:
            config_path = self.clustering_configs[config_name]
            config = ClusteringConfiguration(config_path, self.output_directory)
            newick_path = clustering.order_contigs_simple(config, progress = self.progress)
            clusterings[config_name] = os.path.basename(newick_path)
        self.run.info('clusterings', clusterings)


    def check_args(self):
        if (not self.input_file_path) and (not self.serialized_profile_path):
            raise utils.ConfigError, "You must declare either an input file, or a serialized profile. Use '--help'\
                                      to learn more about the command line parameters."
        if self.input_file_path and self.serialized_profile_path:
            raise utils.ConfigError, "You can't declare both an input file and a serialized profile."
        if self.serialized_profile_path and ('-L' in sys.argv or '--split-length' in sys.argv):
            raise utils.ConfigError, "You can't change split size when you use a serialized profile as an input.\
                                      Unfortunately, the only way to change the split size is to run the profiler\
                                      on the BAM file from scratch."
        if self.serialized_profile_path and (not self.output_directory):
            raise utils.ConfigError, "When loading serialized profiles, you need to declare an output directory."
        if self.input_file_path and not os.path.exists(self.input_file_path):
            raise utils.ConfigError, "No such file: '%s'" % self.input_file_path
        if self.serialized_profile_path and not os.path.exists(self.serialized_profile_path):
            raise utils.ConfigError, "No such file: '%s'" % self.serialized_profile_path
        if not self.min_mean_coverage >= 0:
            raise utils.ConfigError, "Minimum mean coverage must be 0 or larger."
        if not self.min_contig_length >= 0:
            raise utils.ConfigError, "Minimum contig length must be 0 or larger (although using anything\
                                      below 5,000 is kinda silly, UNLESS you are working with mappings of\
                                      multiple samples to a single assembly)."

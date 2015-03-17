# coding: utf-8
"""Provides the necessary class to profile BAM files."""

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
import PaPi.clustering as clustering
import PaPi.annotation as annotation
import PaPi.filesnpaths as filesnpaths
import PaPi.ccollections as ccollections

from PaPi.genes import Genes
from PaPi.metadata import Metadata
from PaPi.contig import Split, Contig, set_contigs_abundance
from PaPi.clusteringconfuguration import ClusteringConfiguration


__author__ = "A. Murat Eren"
__copyright__ = "Copyright 2015, The PaPi Project"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = "0.9.0"
__maintainer__ = "A. Murat Eren"
__email__ = "a.murat.eren@gmail.com"
__status__ = "Development"


pp = terminal.pretty_print


class BAMProfiler:
    """Creates an über class for BAM file operations"""
    def __init__(self, args = None):
        self.args = None
        self.input_file_path = None 
        self.annotation_db_path = None
        self.serialized_profile_path = None 
        self.split_length = 20000
        self.output_directory = None 
        self.list_contigs_and_exit = None 
        self.min_contig_length = 10000 
        self.min_mean_coverage = 10
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
            self.sample_id = args.sample_id

            if args.contigs_of_interest:
                if not os.path.exists(args.contigs_of_interest):
                    raise utils.ConfigError, "Contigs file (%s) is missing..." % (args.contigs_of_interest)

                self.contig_names_of_interest = set([c.strip() for c in open(args.contigs_of_interest).readlines()\
                                                                               if c.strip() and not c.startswith('#')])

        self.bam = None
        self.contigs = {}
        self.genes_in_contigs = {}
        self.annotation_db = None

        self.profile_db = None
        self.profile_db_path = None

        self.clustering_configs = constants.clustering_configs['single']

        self.progress = terminal.Progress()
        self.run = terminal.Run(width=35)

        self.metadata = Metadata(self.progress)


    def init_dirs_and_dbs(self):
        if not self.annotation_db_path:
            raise utils.ConfigError, "You can not run profiling without an annotation database. You can create\
                                      one using 'papi-gen-annotation-database'. Not sure how? Please see the\
                                      user manual."
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

        # know thyself
        self.profile_db.set_meta_value('db_type', 'profile')

        # this will be the unique information that will be passed downstream whenever this db is used:
        # create empty collections tables in newly generated profile database:
        ccollections.create_blank_collections_tables(self.profile_db)

        # put sample id into the meta table
        self.profile_db.set_meta_value('sample_id', self.sample_id)
        # the database design also requires to have a metadata variable that holds all the "samples"
        # described in the database. it doesn't make much sense for an individual profile, but it is
        # well used for merged studies.
        self.profile_db.set_meta_value('samples', self.sample_id)
        self.profile_db.set_meta_value('merged', False)

        self.progress.update('Initializing the annotation database ...')
        self.annotation_db = annotation.AnnotationDatabase(self.annotation_db_path)
        self.split_length = int(self.annotation_db.db.get_meta_value('split_length'))

        self.progress.end()


    def _run(self):
        self.check_args()

        self.set_sample_id()

        if not self.list_contigs_and_exit:
            self.init_dirs_and_dbs()

        # we will set up things here so the information in the annotation_db
        # can be utilized directly from within the contigs for loop. contig to
        # gene associations will be stored in self.genes_in_contigs dictionary for
        # fast access.
        if self.annotation_db and not self.list_contigs_and_exit:
            self.populate_genes_in_contigs_dict()

        self.run.info('profiler_version', __version__)
        self.run.info('sample_id', self.sample_id)
        self.run.info('profile_db', self.profile_db_path)
        self.run.info('annotation_db', True if self.annotation_db_path else False)
        self.run.info('default_clustering', constants.single_default)
        self.run.info('cmd_line', utils.get_cmd_line())
        self.run.info('merged', False)
        self.run.info('split_length', self.split_length)
        self.run.info('min_contig_length', self.min_contig_length)
        self.run.info('min_mean_coverage', self.min_mean_coverage)

        # this is kinda important. we do not run full-blown profile function if we are dealing with a summarized
        # profile...
        if self.input_file_path:
            self.init_profile_from_BAM()
            self.profile()
            self.store_profile()
            self.store_summarized_profile_for_each_split()
        else:
            self.init_serialized_profile()
            self.store_summarized_profile_for_each_split()

        annotation_hash = self.annotation_db.db.get_meta_value('annotation_hash')
        self.profile_db.set_meta_value('annotation_hash', annotation_hash)
        self.run.info('annotation_hash', annotation_hash)
        self.generate_genes_table()
        self.run.info('genes_table', True, quiet = True)

        # here we store both metadata and TNF information into the database:
        self.metadata.store_metadata_for_contigs_and_splits(self.sample_id, self.contigs, self.profile_db)

        # so this is a little sloppy. views variable holds all the table names that are appropriate for visualization.
        # for single runs there is only one table in the PROFILE.db that is relevant for visualization; which is the
        # 'metadata_splits' table. so we set it up here in such a way that it will be seamless to visualize both single
        # and merged runs:
        self.run.info('default_view', 'single', quiet = True)
        self.run.info('views', {'single': 'metadata_splits'}, quiet = True)

        self.store_consenus_FASTA_files_for_splits_and_contigs()

        self.cluster_contigs()

        runinfo_serialized = self.generate_output_destination('RUNINFO.cp')
        self.run.info('runinfo', runinfo_serialized)
        self.run.store_info_dict(runinfo_serialized, strip_prefix = self.output_directory)

        self.run.quit()
        self.profile_db.disconnect()


    def populate_genes_in_contigs_dict(self):
        self.progress.new('Annotation')
        self.progress.update('Reading genes in contigs table')
        genes_in_contigs_table = self.annotation_db.db.get_table_as_dict(annotation.genes_contigs_table_name, annotation.genes_contigs_table_structure)

        self.progress.update('Populating ORFs dictionary for each contig ...')
        for gene in genes_in_contigs_table:
            e = genes_in_contigs_table[gene]
            if self.genes_in_contigs.has_key(e['contig']):
                self.genes_in_contigs[e['contig']].add((gene, e['start'], e['stop']), )
            else:
                self.genes_in_contigs[e['contig']] = set([(gene, e['start'], e['stop']), ])

        self.progress.end()
        self.run.info('annotation_db', "%d genes processed successfully." % len(genes_in_contigs_table), display_only = True)


    def generate_genes_table(self):
        self.genes = Genes(self.progress)
        self.progress.new('Profiling genes')
        num_contigs = len(self.contigs)
        contig_names = self.contigs.keys()
        for i in range(0, num_contigs):
            contig = contig_names[i]
            self.progress.update('Processing contig %d of %d' % (i + 1, num_contigs))
            # if no open reading frames were found in a contig, it wouldn't have an entry in the annotation table,
            # therefore there wouldn't be any record of it in contig_ORFs; so we better check ourselves before
            # we wreck ourselves and the ultimately the analysis of this poor user:
            if self.genes_in_contigs.has_key(contig):
                self.genes.analyze_contig(self.contigs[contig], self.sample_id, self.genes_in_contigs[contig])

        self.genes.create_genes_table(self.profile_db)
        self.progress.end()


    def set_sample_id(self):
        if self.sample_id:
            utils.check_sample_id(self.sample_id)
        else:
            if self.input_file_path:
                self.input_file_path = os.path.abspath(self.input_file_path)
                self.sample_id = os.path.basename(self.input_file_path).upper().split('.BAM')[0]
                self.sample_id = self.sample_id.replace('-', '_')
                self.sample_id = self.sample_id.replace('.', '_')
                if self.sample_id[0] in constants.digits:
                    self.sample_id = 's' + self.sample_id
                utils.check_sample_id(self.sample_id)
            if self.serialized_profile_path:
                self.serialized_profile_path = os.path.abspath(self.serialized_profile_path)
                self.sample_id = os.path.basename(os.path.dirname(self.serialized_profile_path))


    def check_contigs_without_any_ORFs(self, contig_names):
        if not self.annotation_db:
            return
        contig_names = set(contig_names)
        contigs_without_annotation = [c for c in contig_names if c not in self.genes_in_contigs]

        if len(contigs_without_annotation):
            import random
            P = lambda x: 'are %d contigs' % (x) if x > 1 else 'there is one contig'
            self.run.warning('You have instructed profiling to use an annotation database,\
                              however, there %s in your BAM file that did not get annotated. Which means\
                              whatever method you used to identify open reading frames in these contigs\
                              failed to find any open reading frames in those. Which may be normal\
                              (a) if your contigs are very short, or (b) if your gene finder is not\
                              capable of dealing with your stuff. If you know what you are doing, that\
                              is fine. Otherwise please double check. Here is one contig missing\
                              annotation if you would like to play: %s"' %\
                                      (P(len(contigs_without_annotation)), random.choice(contigs_without_annotation)))


    def init_serialized_profile(self):
        self.progress.new('Init')
        self.progress.update('Reading serialized profile')
        self.contigs = dictio.read_serialized_object(self.serialized_profile_path)
        self.progress.end()

        self.run.info('profile_loaded_from', self.serialized_profile_path)
        self.run.info('num_contigs', pp(len(self.contigs)))

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

        # it brings good karma to let the user know what the hell is wrong with their data:
        self.check_contigs_without_any_ORFs(self.contigs.keys())

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

        # store num reads mapped for later use.
        self.profile_db.set_meta_value('total_reads_mapped', int(self.num_reads_mapped))

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


        # it brings good karma to let the user know what the hell is wrong with their data:
        self.check_contigs_without_any_ORFs(self.contig_names)

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
        self.splits = self.annotation_db.db.get_table_as_dict(annotation.splits_info_table_name)
        self.contig_name_to_splits = {}
        for split_name in self.splits:
            parent = self.splits[split_name]['parent']
            if self.contig_name_to_splits.has_key(parent):
                self.contig_name_to_splits[parent].append(split_name)
            else:
                self.contig_name_to_splits[parent] = [split_name]


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

            contig = Contig(contig_name)
            contig.length = self.contig_lenghts[i]
            contig.split_length = self.split_length

            self.progress.new('Profiling "%s" (%d of %d) (%s nts)' % (contig.name,
                                                                      i + 1,
                                                                      len(self.contig_names),
                                                                      pp(int(contig.length))))

            # populate contig with empty split objects and 
            for split_name in self.contig_name_to_splits[contig_name]:
                s = self.splits[split_name]
                split = Split(split_name, contig_name, s['order_in_parent'], s['start'], s['end'])
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


    def store_summarized_profile_for_each_split(self):
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


    def store_consenus_FASTA_files_for_splits_and_contigs(self):
        # generate a sorted list of contigs based on length
        self.contig_names = [t[1] for t in sorted([(self.contigs[k].length, k)\
                                                for k in self.contigs], reverse = True)]

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
        # FIXME: need a profiledb.py to do all these stuff and the things that are done in the metadata.py:
        clusterings_table_name      = 'clusterings'
        clusterings_table_structure = ['clustering', 'newick' ]
        clusterings_table_types     = [   'str'    ,  'str'   ]
        self.profile_db.create_table(clusterings_table_name, clusterings_table_structure, clusterings_table_types)

        clusterings = []

        for config_name in self.clustering_configs:
            config_path = self.clustering_configs[config_name]

            config = ClusteringConfiguration(config_path, self.output_directory, version = __version__)
            newick = clustering.order_contigs_simple(config, progress = self.progress)

            clusterings.append(config_name)
            db_entries = tuple([config_name, newick])
            self.profile_db._exec('''INSERT INTO %s VALUES (?,?)''' % clusterings_table_name, db_entries)

        self.run.info('available_clusterings', clusterings)


    def check_args(self):
        if (not self.input_file_path) and (not self.serialized_profile_path):
            raise utils.ConfigError, "You must declare either an input file, or a serialized profile. Use '--help'\
                                      to learn more about the command line parameters."
        if self.input_file_path and self.serialized_profile_path:
            raise utils.ConfigError, "You can't declare both an input file and a serialized profile."
        if self.serialized_profile_path and (not self.output_directory):
            raise utils.ConfigError, "When loading serialized profiles, you need to declare an output directory."
        if self.input_file_path and not os.path.exists(self.input_file_path):
            raise utils.ConfigError, "No such file: '%s'" % self.input_file_path
        if self.serialized_profile_path and not os.path.exists(self.serialized_profile_path):
            raise utils.ConfigError, "No such file: '%s'" % self.serialized_profile_path
        if not self.min_mean_coverage >= 0:
            raise utils.ConfigError, "Minimum mean coverage must be 0 or larger."
        if not self.min_contig_length >= 0:
            raise utils.ConfigError, "Minimum contig length must be 0 or larger."

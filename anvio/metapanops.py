# -*- coding: utf-8
# pylint: disable=line-too-long
"""
    The anvi'o metapangenome module.

    anvi-meta-pan-genome is the default client using this.
"""

import numpy as numpy

import anvio
import anvio.dbops as dbops
import anvio.utils as utils
import anvio.terminal as terminal
import anvio.summarizer as summarizer
import anvio.genomedescriptions as genomedescriptions

from anvio.errors import ConfigError
from anvio.tables.miscdata import TableForLayerAdditionalData, TableForItemAdditionalData


__copyright__ = "Copyleft 2015-2018, the Meren Lab (http://merenlab.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "A. Murat Eren"
__email__ = "a.murat.eren@gmail.com"



class MetaPangenome(object):
    def __init__(self, args=None, run=terminal.Run(), progress=terminal.Progress()):
        self.args = args
        self.run = run
        self.progress = progress

        A = lambda x: args.__dict__[x] if x in args.__dict__ else None

        if len([p for p in [A('pan_db'), A('genomes_storage')] if not p]):
            raise ConfigError("MetaPangenome class should be inherited with an `args` object that contains "
                              "an anvi'o pan database (`pan_db`), genomes storage (`genomes_storage`) :/")

        self.pan_db_path = A('pan_db')
        self.genomes_storage_path = A('genomes_storage')
        self.num_threads = A('num_threads')

        self.fraction_of_median_coverage = A('fraction_of_median_coverage') or 0.25
        self.min_detection = A('min_detection') or 0.50

        # This object will be populated to give access to pan summary:
        self.pan_summary = None

        # This object will give access to genome descriptions, and will know everything about
        # our internal genomes. during the initialization of genomes, we will also recover the
        # sample names stored in profile databases.
        self.descriptions = None
        self.init_genome_descriptions()

        # go through each profile database and set sample names.
        self.sample_names = None
        self.set_sample_names()

        # This dict makes sure various operations in this class do not initialize
        # the summary object multiple times for genomes described in the
        # same profile database. ALTHOUGH, we currently allow only one profile database and
        # one collection (see init_genome_descriptions function). There is no theoretical
        # reason to not change this. But for now, we will keep it simple simply becasuse we
        # don't have enough test material.
        self.unique_profile_db_path_to_internal_genome_name = None
        self.set_unique_profile_db_path_to_internal_genome_name_dict()


    def init_pan_summary(self):
        self.progress.new('Pan summary')
        self.progress.update('Initializing...')

        args = summarizer.ArgsTemplateForSummarizerClass()
        args.pan_db = self.pan_db_path
        args.genomes_storage = self.genomes_storage_path
        args.skip_check_collection_name = True
        args.skip_init_functions = True
        args.output_dir = None

        self.pan_summary = summarizer.PanSummarizer(args)

        self.progress.end()


    def init_genome_descriptions(self):
        self.progress.new('Genome descriptions')
        self.progress.update('Initializing')

        self.descriptions = genomedescriptions.GenomeDescriptions(self.args)

        if len(self.descriptions.external_genomes_dict):
            raise ConfigError("Anvi'o doesn't know how you did it, but you managed to inherit this class with an "
                              "`args` object that describes external genomes. Unfortunately anvi'o metapangneomic "
                              "workflow only works with internal genomes since it is all about making sense of "
                              "pangenomes in the context of metageomes. So this is not really working for us :(")

        if len(set([g['profile_db_path'] for g in list(self.descriptions.internal_genomes_dict.values())])) > 1:
            raise ConfigError("There are multiple profile databases in your internal genomes file. We are simply "
                              "not ready to deal with this complexity. If you think this is a mistake, let us know "
                              "and we will work with you to make anvi'o work with multiple profile databases (in "
                              "fact anvi'o is able to make sense of internal genomes across multiple profile "
                              "databases, but we haven't tested it to understand potential caveats associated "
                              "with that level of complexity).")

        if len(set([g['collection_id'] for g in list(self.descriptions.internal_genomes_dict.values())])) > 1:
            raise ConfigError("For the sake of simplicity, we expect collection names to be identical in a given "
                              "internal genomes file. Anvi'o is kind of a paranoid, and it apologizes for it.")

        self.descriptions.load_genomes_descriptions(skip_functions=True)

        self.progress.end()

        self.run.info("Internal genomes found", "%d (%s)" % (len(self.descriptions.internal_genome_names), ', '.join(self.descriptions.internal_genome_names)))


    def set_sample_names(self):
        """Go through all profile databases involved, and learn all sample names"""

        self.sample_names = []

        for profile_db_path in set([g['profile_db_path'] for g in list(self.descriptions.internal_genomes_dict.values())]):
            self.sample_names.extend(sorted(list(dbops.ProfileDatabase(profile_db_path).samples)))

        self.run.info("Samples found", "%d (%s)" % (len(self.sample_names), ', '.join(self.sample_names)), nl_after=1)


    def get_summary_object_for_profile_db(self, profile_db_path, init_gene_coverages=True):
        collection_name = self.descriptions.genomes[self.unique_profile_db_path_to_internal_genome_name[profile_db_path][0]]['collection_id']
        profile_db_path = self.descriptions.genomes[self.unique_profile_db_path_to_internal_genome_name[profile_db_path][0]]['profile_db_path']
        contigs_db_path = self.descriptions.genomes[self.unique_profile_db_path_to_internal_genome_name[profile_db_path][0]]['contigs_db_path']

        # poor-man's whatever
        bin_names_list = [self.descriptions.genomes[g]['bin_id'] for g in self.unique_profile_db_path_to_internal_genome_name[profile_db_path]]

        ARGS = summarizer.ArgsTemplateForSummarizerClass()
        ARGS.profile_db = profile_db_path
        ARGS.contigs_db = contigs_db_path
        ARGS.skip_init_functions = True
        ARGS.init_gene_coverages = init_gene_coverages
        ARGS.collection_name = collection_name
        ARGS.bin_names_list = bin_names_list
        ARGS.output_dir = None

        summary = summarizer.ProfileSummarizer(ARGS)
        summary.init()
        summary.init_collection_profile(collection_name)

        return summary


    def get_genomes_across_metagenomes_dict(self, data_key='mean_coverage'):
        self.progress.new('Recovering data for genomes across metagenomes')
        self.progress.update('...')

        genomes_across_metagenomes = {}
        for internal_genome_name in self.descriptions.internal_genome_names:
            genome_name = self.descriptions.genomes[internal_genome_name]['bin_id']
            genomes_across_metagenomes[internal_genome_name] = {}
            for sample_name in self.sample_names:
                genomes_across_metagenomes[internal_genome_name][sample_name] = 0.0

        D = lambda: summary.collection_profile[genome_name][data_key]
        for profile_db_path in self.unique_profile_db_path_to_internal_genome_name:
            self.progress.update('"%s" from profile db at %s ...' % (data_key, profile_db_path))
            summary = self.get_summary_object_for_profile_db(profile_db_path, init_gene_coverages=False)

            for internal_genome_name in self.unique_profile_db_path_to_internal_genome_name[profile_db_path]:
                genome_name = self.descriptions.genomes[internal_genome_name]['bin_id']
                for sample_name in D():
                    genomes_across_metagenomes[internal_genome_name][sample_name] = D()[sample_name]

        self.progress.end()

        return genomes_across_metagenomes


    def set_unique_profile_db_path_to_internal_genome_name_dict(self):
        self.unique_profile_db_path_to_internal_genome_name = self.descriptions.get_unique_profile_db_path_to_internal_genome_name_dict()

        for profile_db_path in self.unique_profile_db_path_to_internal_genome_name:
            collection_names = set([self.descriptions.genomes[genome_name]['collection_id'] for genome_name in self.unique_profile_db_path_to_internal_genome_name[profile_db_path]])
            if len(collection_names) != 1:
                self.progress.end()
                raise ConfigError("You have to have the same collection for each bin originate from the same profile db.")


    def get_gene_presence_in_the_environment_dict(self):
        if not isinstance(self.fraction_of_median_coverage, float):
            raise ConfigError("Fraction of median coverage must of type `float`.")

        if not isinstance(self.min_detection, float):
            raise ConfigError("Minimum detection must be of type `float`")

        self.run.info('Fraction of median coverage for core genes', self.fraction_of_median_coverage)
        self.run.info('Min detection of a genome in at last one metagenome', self.min_detection)

        self.progress.new('Working on gene presence/absence')
        self.progress.update('...')

        gene_presence_in_the_environment_dict = {}
        for profile_db_path in self.unique_profile_db_path_to_internal_genome_name:
            self.progress.update('Collection info from profile db at %s ...' % (profile_db_path))
            summary = self.get_summary_object_for_profile_db(profile_db_path)

            for internal_genome_name in self.unique_profile_db_path_to_internal_genome_name[profile_db_path]:
                genome_name = self.descriptions.genomes[internal_genome_name]['bin_id']

                self.progress.update('Working on genome %s in profile db %s ...' % (internal_genome_name, profile_db_path))

                # for each genome, first we will see whether it is detected in at least one metagenome
                detection_across_metagenomes = summary.collection_profile[genome_name]['detection']
                num_metagenomes_above_min_detection = [m for m in detection_across_metagenomes if detection_across_metagenomes[m] > self.min_detection]
                not_enough_detection = False if len(num_metagenomes_above_min_detection) else True

                gene_presence_in_the_environment_dict[genome_name] = {}
                split_names_of_interest = self.descriptions.get_split_names_of_interest_for_internal_genome(self.descriptions.genomes[internal_genome_name])

                genome_bin_summary = summarizer.Bin(summary, genome_name, split_names_of_interest)
                gene_coverages_across_samples = utils.get_values_of_gene_level_coverage_stats_as_dict(genome_bin_summary.gene_level_coverage_stats_dict, "mean_coverage")

                # at this point we have all the genes in the genome bin. what we need is to characterize their detection. first,
                # summarize the coverage of each gene in all samples:
                sum_gene_coverages_across_samples = dict([(gene_callers_id, sum(gene_coverages_across_samples[gene_callers_id].values())) for gene_callers_id in gene_coverages_across_samples])

                # now we will identify the median coverage
                median_coverage_across_samples = numpy.median(list(sum_gene_coverages_across_samples.values()))

                # now we will store decide whether a gene found in this genome is also found in the environment, and store that
                # information into `gene_presence_in_the_environment_dict`, and move on to the next stage.
                for gene_caller_id in sum_gene_coverages_across_samples:
                    if not_enough_detection:
                        _class = 'NA'
                    elif sum_gene_coverages_across_samples[gene_caller_id] < median_coverage_across_samples * self.fraction_of_median_coverage:
                        _class = 'EAG'
                    else:
                        _class = 'ECG'

                    gene_presence_in_the_environment_dict[genome_name][gene_caller_id] = _class

        self.progress.end()

        return gene_presence_in_the_environment_dict


    def add_genomes_across_metagenomes_dict_into_pan_database(self):
        genomes_across_metagenomes_dict = self.get_genomes_across_metagenomes_dict()

        self.args.just_do_it = True
        TableForLayerAdditionalData(self.args).add(genomes_across_metagenomes_dict, self.sample_names)


    def add_ECG_EAG_ratio_per_gene_cluster_into_pan_database(self):
        if not self.pan_summary:
            self.init_pan_summary()

        gene_presence_in_the_environment_dict = self.get_gene_presence_in_the_environment_dict()

        self.progress.new('Working on ECG/EAG ratio per gene cluster')
        self.progress.update('...')

        gene_status_frequencies_in_gene_cluster = {}

        gene_cluster_names = list(self.pan_summary.gene_clusters.keys())
        num_gene_clusters = len(gene_cluster_names)
        for i in range(0, num_gene_clusters):
            self.progress.update('%.2f' % ((i + 1) * 100 / num_gene_clusters))
            gene_cluster_name = gene_cluster_names[i]

            status = {'EAG': 0, 'ECG': 0, 'NA': 0}
            for internal_genome_name in self.pan_summary.gene_clusters[gene_cluster_name]:
                genome_name = self.descriptions.genomes[internal_genome_name]['bin_id']

                for gene_caller_id in self.pan_summary.gene_clusters[gene_cluster_name][internal_genome_name]:
                    if genome_name not in gene_presence_in_the_environment_dict:
                        self.progress.end()
                        raise ConfigError("Something is wrong... It seems you generated a pangenome with an internal genomes file "
                                          "that is not identical to the internal genomes file you are using to run this program.")

                    status[gene_presence_in_the_environment_dict[genome_name][gene_caller_id]] += 1
            gene_status_frequencies_in_gene_cluster[gene_cluster_name] = status

        # setup some boring variable names.
        items_additional_data_dict = {}
        key_ECG_EAG_ratio = 'EAG_ECG_ratio'
        key_ECGs_and_EAGs = 'ECGs_and_EAGs'
        list_ECG_EAG_keys = ['EAG', 'ECG', 'NA']

        self.progress.update('Setting up the items data dictionary ..')
        for gene_cluster_name in gene_status_frequencies_in_gene_cluster:
            r = gene_status_frequencies_in_gene_cluster[gene_cluster_name]

            # add ECG and EAG frequencies for the gene cluster
            items_additional_data_dict[gene_cluster_name] = dict([('%s!%s' % (key_ECGs_and_EAGs, status), r[status]) for status in list_ECG_EAG_keys])

            # add ECG / EAG ratio
            items_additional_data_dict[gene_cluster_name][key_ECG_EAG_ratio] = (r['EAG'] / (r['EAG'] + r['ECG']) if (r['EAG'] + r['ECG']) else 0)

        self.progress.end()

        # add that bad boy to the database
        self.args.just_do_it = True
        items_additional_data_keys = [('%s!%s' % (key_ECGs_and_EAGs, status)) for status in list_ECG_EAG_keys] + [key_ECG_EAG_ratio]
        TableForItemAdditionalData(self.args).add(items_additional_data_dict, items_additional_data_keys)


    def process(self):
        """Annotates the pan database with metapangenomic information"""

        self.add_genomes_across_metagenomes_dict_into_pan_database()

        self.add_ECG_EAG_ratio_per_gene_cluster_into_pan_database()


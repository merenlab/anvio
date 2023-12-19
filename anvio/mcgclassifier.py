# -*- coding: utf-8
# pylint: disable=line-too-long
"""
    Classes to classify genes based on coverages across metagenomes.

    anvi-mcg-classifier is the default client using this module
"""


import os
import anvio
import numpy as np
import pandas as pd
import matplotlib
# TODO: according to the warning, this call to set the back-hand is meaningless
# I need to experiment to see what happens if I delete it.
matplotlib.use('pdf')
import anvio.utils as utils
import matplotlib.pyplot as plt
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths

from scipy import odr as odr
from anvio.mcgops import MCGPlots
from anvio.errors import ConfigError, FilesNPathsError
from matplotlib.backends.backend_pdf import PdfPages


__author__ = "Developers of anvi'o (see AUTHORS.txt)"
__copyright__ = "Copyleft 2015-2018, the Meren Lab (http://merenlab.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "Alon Shaiber"
__email__ = "alon.shaiber@gmail.com"


run = terminal.Run()
progress = terminal.Progress()
pp = terminal.pretty_print


columns_for_samples_coverage_stats_dict = ['non_outlier_mean_coverage', 'non_outlier_coverage_std']


class MetagenomeCentricGeneClassifier:
    def __init__(self, args, run=run, progress=progress):
        self.run = run
        self.progress = progress

        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        self.output_file_prefix = A('output_file_prefix')
        self.alpha = A('alpha')
        self.collection_name = A('collection_name')
        self.bin_id = A('bin_id')
        self.bin_ids_file_path = A('bin_ids_file')
        self.exclude_samples = A('exclude_samples')
        self.include_samples = A('include_samples')
        self.outliers_threshold = A('outliers_threshold')
        self.zeros_are_outliers = A('zeros_are_outliers')
        self.gen_figures = A('gen_figures')
        self.overwrite_output_destinations = A('overwrite_output_destinations')
        self.split_coverage_values_per_nt_dict = None
        self.gene_level_coverage_stats_dict = None
        self.gene_level_coverage_stats_dict_of_dataframes = None


        self.profile_db = {}
        self.coverage_values_per_nt = None
        self.gene_coverages = {}
        self.gene_detections = None
        self.samples = None
        self.positive_samples = []
        self.number_of_positive_samples = None
        self.negative_samples = {}
        self.number_of_negative_samples = None
        self.gene_class_df = {}
        self.samples_detection_information = {}
        self.gene_presence_absence_in_samples_initiated = False
        self.gene_presence_absence_in_samples = None
        self.additional_description = ''
        self.total_length = None
        self.samples_coverage_stats_dicts_was_initiated = False
        self.samples_coverage_stats_dicts = {}
        self.non_outlier_indices = {}
        self.gene_coverage_consistency_dict = {}
        self.gene_coverage_consistency_dict_initiated = False
        self.samples_to_exclude = set([])
        self.samples_to_include = set([])
        self.write_output_to_files = None

        if self.exclude_samples:
            # check that there is a file like this
            filesnpaths.is_file_exists(self.exclude_samples)
            self.samples_to_exclude = set([l.split('\t')[0].strip() for l in open(self.exclude_samples, 'r').readlines()])

            if not self.samples_to_exclude:
                raise ConfigError("You asked to exclude samples, but provided an empty list.")

            run.info('Excluding Samples', 'The following samples will be excluded: %s' % self.samples_to_exclude,)

        if self.include_samples:
            # check that there is a file like this
            filesnpaths.is_file_exists(self.include_samples)
            self.samples_to_include = set([l.split('\t')[0].strip() for l in open(self.include_samples, 'r').readlines()])

            if not self.samples_to_include:
                raise ConfigError("You provided an empty list of samples to include.")

            run.info('Including Samples', 'The following samples will be included: %s' % self.samples_to_include,)

        # run sanity check on all input arguments
        self.sanity_check()


    def init(self, gene_level_coverage_stats_dict=None, split_coverage_values_per_nt_dict=None, additional_description=None):
        """ setting the dictionaries for gene coverage stats and for split coverage per nucleotide"""

        if gene_level_coverage_stats_dict is None and split_coverage_values_per_nt_dict is None:
            raise ConfigError("MCGC needs at least one of the following in order to work: "
                               "gene_level_coverage_stats_dict or/and split_coverage_values_per_nt_dict")

        # We want to make sure these are empty in case we use "init" multiple times for different bins
        self.coverage_values_per_nt = None
        self.gene_class_df = {}
        self.samples_detection_information = {}
        self.gene_presence_absence_in_samples_initiated = False
        self.gene_presence_absence_in_samples = None
        self.samples_coverage_stats_dicts_was_initiated = False
        self.samples_coverage_stats_dicts = {}
        self.non_outlier_indices = {}
        self.gene_coverage_consistency_dict = {}
        self.gene_coverage_consistency_dict_initiated = False

        self.gene_level_coverage_stats_dict = gene_level_coverage_stats_dict
        self.split_coverage_values_per_nt_dict = split_coverage_values_per_nt_dict
        if additional_description:
            self.additional_description = '-' + additional_description

        try:
            samples = next(iter(self.gene_level_coverage_stats_dict.values())).keys()
        except:
            samples = next(iter(self.split_coverage_values_per_nt_dict.values())).keys()
        self.init_samples(samples)


    def sanity_check(self):
        """Basic sanity check for class inputs"""

        if self.output_file_prefix:
            filesnpaths.is_output_file_writable(self.output_file_prefix + '-additional-layers.txt', ok_if_exists=self.overwrite_output_destinations)

        try:
            if self.gen_figures:
                plot_dir = self.output_file_prefix + '-nucleotide-coverage-distribution-plots'
                os.makedirs(plot_dir, exist_ok=self.overwrite_output_destinations)
        except FileExistsError as e:
            raise FilesNPathsError("%s already exists, if you would like to overwrite it, then use -W (see help menu)." % plot_dir)

        # checking alpha
        if not isinstance(self.alpha, float):
            raise ConfigError("alpha value must be a type float.")
        # alpha must be a min of 0 and smaller than 0.5
        if self.alpha < 0 or self.alpha >= 0.5:
            raise ConfigError("alpha must be a minimum of 0 and smaller than 0.5")

        if self.exclude_samples and self.include_samples:
            raise ConfigError("You cannot use both --include-samples and --exclude-samples! Please choose one.")


    def init_samples(self, samples_list):
        """ Create the set of samples according to user input and store it in self.samples"""
        # remove the samples that should be excluded
        samples = set(samples_list) - self.samples_to_exclude

        if self.include_samples:
            samples_to_include_that_are_not_there = self.samples_to_include - samples
            if samples_to_include_that_are_not_there:
                raise ConfigError("You requested to include some samples that are not in the profile database. Here are the samples in the profile database: %s. "
                               "And here are the samples you requested, and that are not there: %s" % (samples, samples_to_include_that_are_not_there))
            samples = self.samples_to_include

        self.samples = samples


    def init_gene_level_coverage_stats_dict_of_dataframes(self):
        """ converts the dictionaries of gene_level_coverage_stats_dict to dataframes"""
        self.gene_level_coverage_stats_dict_of_dataframes = {}
        for key in ['mean_coverage', 'detection', 'non_outlier_mean_coverage', 'non_outlier_coverage_std']:
            # Only include samples that the user want
            gene_stat = utils.get_values_of_gene_level_coverage_stats_as_dict(self.gene_level_coverage_stats_dict, key, as_pandas=True, samples_of_interest=self.samples)
            self.gene_level_coverage_stats_dict_of_dataframes[key] = gene_stat

        for key in ['gene_coverage_values_per_nt', 'non_outlier_positions']:
            gene_stat = utils.get_values_of_gene_level_coverage_stats_as_dict(self.gene_level_coverage_stats_dict, key, as_pandas=False, samples_of_interest=self.samples)
            self.gene_level_coverage_stats_dict_of_dataframes[key] = gene_stat


    def init_samples_coverage_stats_dict(self):
        """ populate the samples_coverage_stats_dict, and determine positive, negative, and ambiguous samples with the genome detection information
            (--alpha, --genome-detection-uncertainty)

            The samples_coverage_stats_dict dataframe is used to calculate the gene consistency information.
            It is also used for plotting purposes (both for the nucleotide-coverage-distribution plots and the gene-consistency plots).

            The coverage_values_per_nt is used to calculate the detection value (portion of nucleotides
            covered) for a sample. Then, a cutoff for detection values is used to determine the presence
            or absence of the genome in each sample.
        """

        if self.coverage_values_per_nt is None:
            self.coverage_values_per_nt = get_coverage_values_per_nucleotide(self.split_coverage_values_per_nt_dict, samples=self.samples)

        total_length = len(next(iter(self.coverage_values_per_nt.values())))
        MCG_samples_information_table_structure = ['samples', 'presence', 'detection', 'number_of_taxon_specific_core_detected']

        # create an empty dataframe
        samples_information = pd.DataFrame(index=self.samples, columns=MCG_samples_information_table_structure[1:])
        positive_samples = []
        negative_samples = []

        self.progress.new("Finding nucleotide positions in samples with outlier coverage values")
        progress.update('...')
        num_samples, counter = len(self.samples), 1
        detection = {}
        total_length = len(next(iter(self.coverage_values_per_nt.values())))

        self.samples_coverage_stats_dicts = pd.DataFrame(index=self.samples, columns=columns_for_samples_coverage_stats_dict)
        for sample in self.samples:
            if num_samples > 100 and counter % 100 == 0:
                self.progress.update('%d of %d samples...' % (counter, num_samples))
            # get the non-outlier information
            non_outlier_indices, self.samples_coverage_stats_dicts.loc[sample,] = get_non_outliers_information(self.coverage_values_per_nt[sample], MAD_threshold=self.outliers_threshold, zeros_are_outliers=self.zeros_are_outliers)
            self.non_outlier_indices[sample] = non_outlier_indices
            number_of_non_outliers = len(self.non_outlier_indices[sample])
            if anvio.DEBUG:
                self.run.info_single('The mean and std of non-outliers in sample %s are: %s, %s respectively' % (sample, self.samples_coverage_stats_dicts['non_outlier_mean_coverage'][sample], self.samples_coverage_stats_dicts['non_outlier_coverage_std'][sample]))
                self.run.info_single('The number of non-outliers is %s of %s (%.2f%%)' % (number_of_non_outliers, total_length, 100.0 * number_of_non_outliers / total_length))
            detection[sample] = np.count_nonzero(self.coverage_values_per_nt[sample]) / total_length
            samples_information['presence'][sample] = get_presence_absence_information(number_of_non_outliers/total_length, self.alpha)
            if detection[sample] <= 0.5:
                samples_information['presence'][sample] = False
            if samples_information['presence'][sample]:
                positive_samples.append(sample)
            elif samples_information['presence'][sample] == False:
                negative_samples.append(sample)

            samples_information['detection'][sample] = detection[sample]
            counter += 1

        self.positive_samples = positive_samples
        self.number_of_positive_samples = len(self.positive_samples)
        self.negative_samples = negative_samples
        self.samples_detection_information = samples_information
        self.run.warning('The number of positive samples is %s' % self.number_of_positive_samples)
        self.run.warning('The number of negative samples is %s' % len(self.negative_samples))


        self.samples_coverage_stats_dicts_was_initiated = True
        self.progress.end()


    def plot_nucleotide_coverage_distribution(self):
        """ Creates a pdf file with the following plots for each sample the sorted nucleotide coverages \
        (with the outliers in red and non-outliers in blue), and a histogram of coverages for the non-outliers"""
        # Creating a dircetory for the plots. If running on bins, each bin would be in a separate sub-directory

        if not self.samples_coverage_stats_dicts_was_initiated:
            self.init_samples_coverage_stats_dict()

        plot_dir = self.output_file_prefix + '-nucleotide-coverage-distribution-plots' + '/'
        self.progress.new('Saving figures of taxon specific distributions to pdf')
        progress.update('...')
        number_of_fininshed = 0
        for sample in self.positive_samples:
            coverages_pdf_output = plot_dir + sample + self.additional_description + '-coverages.pdf'
            pdf_output_file = PdfPages(coverages_pdf_output)
            v = self.coverage_values_per_nt[sample]
            # Using argsort so we can use the non_oulier indices
            sorting_indices = np.argsort(v)
            # we would need the reverse of the sorting of the indices to create the x axis for the non-outliers
            reverse_sorted_indices = np.zeros(len(sorting_indices))
            reverse_sorted_indices[sorting_indices] = range(len(reverse_sorted_indices))

            # plotting the ordered coverage values (per nucleotide)
            # the non-outliers are plotted in blue
            # the outlier values are plotted in red
            fig = plt.figure()
            ax = fig.add_subplot(111, rasterized=True)
            ax.set_xlabel = 'Nucleotide Number (ordered)'
            ax.set_ylabel = r'$Nucleotide Coverage^2$'
            x1 = range(len(v)) # FIXME: this shouldn't be in the loop (only here because I need to fix the mock data)
            x2 = reverse_sorted_indices[self.non_outlier_indices[sample]]
            #y2 = v[self.non_outlier_indices[sample]]
            # plot all in red
            ax.semilogy(x1,v[sorting_indices],'r.', rasterized=True)
            # plot on top the non-outliers in blue
            ax.semilogy(x2,v[self.non_outlier_indices[sample]],'b.', rasterized=True)
            fig.suptitle("%s - sorted coverage values with outliers" % sample)
            plt.savefig(pdf_output_file, format='pdf')
            plt.close()

            # plotting a histogram of the non-outliers
            # This would allow to see if they resemble a normal distribution
            hist_range = (min(v[self.non_outlier_indices[sample]]),max(v[self.non_outlier_indices[sample]]))
            # computing the number of bins so that the width of a bin is ~1/4 of the standard deviation
            # FIXME: need to make it so the bins are only of integers (so the smallest bin is of width 1
            # and that bins are integers)
            number_of_hist_bins = np.ceil((hist_range[1] - hist_range[0]) / (self.samples_coverage_stats_dicts['non_outlier_coverage_std'][sample]/4)).astype(int) # setting the histogram bins to be of the width of a quarter of std
            fig = plt.figure()
            ax = fig.add_subplot(111, rasterized=True)
            ax.set_xlabel = 'Coverage'
            ax.hist(v[self.non_outlier_indices[sample]], number_of_hist_bins,hist_range, rasterized=True)
            fig.suptitle("%s - histogram of non-outliers" % sample)
            # adding the mean and std of the non-outliers as text to the plot
            text_for_hist = u'$\mu = %d$\n $\sigma = %d$' %\
                                (self.samples_coverage_stats_dicts['non_outlier_mean_coverage'][sample],\
                                 self.samples_coverage_stats_dicts['non_outlier_coverage_std'][sample])
            ax.text(0.8, 0.9, text_for_hist, ha='center', va='center', transform=ax.transAxes)
            plt.savefig(pdf_output_file, format='pdf')
            plt.close()
            # close the pdf file
            pdf_output_file.close()
            number_of_fininshed += 1
            self.progress.update("Finished %d of %d" % (number_of_fininshed, self.number_of_positive_samples))
        self.progress.end()


    def init_gene_presence_absence_in_samples(self):
        """ Determining presence and absence of genes in samples according to gene detection values."""
        if not self.gene_level_coverage_stats_dict:
            raise ConfigError("gene presence/absence in samples cannot be determined without a gene_level_coverage_stats_dict,\
                                but it seems that you don't have one. maybe you should run init()?")

        if self.gene_level_coverage_stats_dict_of_dataframes is None:
            self.init_gene_level_coverage_stats_dict_of_dataframes()

        gene_callers_id = self.gene_level_coverage_stats_dict_of_dataframes['detection'].index
        self.gene_presence_absence_in_samples = pd.DataFrame(index=gene_callers_id, columns=self.samples)

        T = lambda x: get_presence_absence_information(sum(x)/len(x), self.alpha)
        self.progress.new('Computing gene presence/absence in samples')
        progress.update('...')
        genes_above_outlier_threshold = pd.DataFrame.from_dict(self.gene_level_coverage_stats_dict_of_dataframes['non_outlier_positions'], orient='index').applymap(T)
        genes_with_detection_above_half = self.gene_level_coverage_stats_dict_of_dataframes['detection'].applymap(lambda x: x > 0.5)
        self.gene_presence_absence_in_samples = genes_above_outlier_threshold & genes_with_detection_above_half
        self.gene_presence_absence_in_samples_initiated = True
        self.progress.end()


    def init_gene_coverage_consistency_information(self):
        """ Perform orthogonal distance regression for each gene to determine coverage consistency.

            The question that we are trying to ask is:
                Do the non-outlier nt coverage of the gene in samlpes correlates to the non-outlier
                nt coverage of the genome in samples?

            The regression is performed only for positive samples.
            For each gene, the regression is performed only according to samples in which
            the gene is present (according to the detection critrea).
        """
        if not self.samples_coverage_stats_dicts_was_initiated:
            self.init_samples_coverage_stats_dict()

        if not self.gene_presence_absence_in_samples_initiated:
            self.init_gene_presence_absence_in_samples()

        self.progress.new("Computing coverage consistency for all genes.")
        progress.update('...')
        gene_ids = self.gene_level_coverage_stats_dict_of_dataframes['mean_coverage'].index
        num_genes, counter = len(gene_ids), 1
        for gene_id in gene_ids:
            if num_genes > 100 and counter % 100 == 0:
                self.progress.update('%d of %d genes...' % (counter, num_genes))

            # samples in which the gene is present
            _samples = self.gene_presence_absence_in_samples.loc[gene_id,self.gene_presence_absence_in_samples.loc[gene_id,]==True].index
            # mean and std of non-outlier nt in each sample
            x = list(self.samples_coverage_stats_dicts.loc[_samples,'non_outlier_mean_coverage'].values)
            if "non_outlier_coverage_std" in self.samples_coverage_stats_dicts:
                # we only expect to have the sample coverage std in "full" mode
                std_x = list(self.samples_coverage_stats_dicts.loc[_samples,'non_outlier_coverage_std'].values)
            else:
                std_x = None

            if len(_samples) > 1:
                # mean and std of non-outlier nt in the gene (in each sample)
                y = self.gene_level_coverage_stats_dict_of_dataframes['non_outlier_mean_coverage'].loc[gene_id, _samples].values
                std_y = self.gene_level_coverage_stats_dict_of_dataframes['non_outlier_coverage_std'].loc[gene_id, _samples].values

                # performing the regression using ODR
                _data = odr.RealData(x, y, std_x, std_y)
                _model = lambda B, c: B[0] * c
                _odr = odr.ODR(_data, odr.Model(_model), beta0=[3])
                odr_output = _odr.run()

                # store results
                self.gene_coverage_consistency_dict[gene_id] = {}
                self.gene_coverage_consistency_dict[gene_id]['slope'] = odr_output.beta[0]
                self.gene_coverage_consistency_dict[gene_id]['slope_std'] = odr_output.sd_beta[0]
                self.gene_coverage_consistency_dict[gene_id]['slope_precision'] = odr_output.sd_beta[0] / odr_output.beta[0]

                # compute R squered
                f = lambda b: lambda _x: b*_x
                R_squered = 1 - sum((np.apply_along_axis(f(odr_output.beta[0]),0,x)-y)**2) / sum((y-np.mean(y))**2)

                # Check if converged
                self.gene_coverage_consistency_dict[gene_id]['R_squered'] = R_squered
                if odr_output.stopreason[0] == 'Sum of squares convergence':
                    self.gene_coverage_consistency_dict[gene_id]['converged'] = True
                else:
                    self.gene_coverage_consistency_dict[gene_id]['converged'] = False

        self.gene_coverage_consistency_dict_initiated = True
        self.progress.end()


    def get_gene_specificity(self, gene_id):
        """ return True for gene if it occurs in positive samples and doesn't occur in negative samples.

            Ambiguous occurences are not counted as anything. This means that if a gene is ambiguously
            occuring in a negative sample it could still be counted as "specific". It also means that
            if a gene is only ambiguously occuring in positive samples then it would be considered
            as "non-specific".
        """

        if self.gene_class_df.loc[gene_id, 'occurence_in_positive_samples'] > 1 and self.gene_class_df.loc[gene_id, 'occurence_in_negative_samples'] == 0:
            return True
        else:
            return False
        # TODO: if there are no occurences of the gene at all, then we should maybe return None instead of False


    def get_gene_coverage_consistency(self, gene_id):
        """ return true if the gene's coverage is consistent accross positive samples, False otherwise."""

        # TODO: make sure coverage_consistency_dict has been initiated
        if self.gene_class_df.loc[gene_id, 'occurence_in_positive_samples'] == 0:
            # if the gene doesn't occur in positive samlpes then there is no classification
            return None
        elif self.gene_class_df.loc[gene_id, 'occurence_in_positive_samples'] == 1:
            # if the gene occurs only in one positive sample then return True.
            # XXX: we might prefer to return None, we should consider this in the future.
            return True
        elif self.gene_coverage_consistency_dict[gene_id]['converged']:
                # FIXME: this is where we use an arbitrary threshold again :-(
                # if the slope precision is smaller than the threshold then the regression
                # fit is considered accurate enough and the gene coverage is considered consistent.
            return self.gene_coverage_consistency_dict[gene_id]['slope_precision'] < 0.5
        else:
            # The regression didn't converege so the coverage is probably not consistent.
            return False


    def determine_if_gene_is_core(self, gene_id, gene_specificity):
        """ return True for core gene, False for accessory gene

            If the gene is specific to positive samples, then core would be considered if it
            occurs in all positive samples. Otherwise it would be considered core if it
            occurs in all positive AND all negative samples.
            Ambiguous occurences of a gene are not considered (i.e. they are the same as absence).
        """

        if gene_specificity:
            # return True if the gene occurs in all positive samples.
            return self.gene_class_df.loc[gene_id, 'occurence_in_positive_samples'] == len(self.positive_samples)
        else:
            # return True if the gene occurs in all positive AND all negative samples
            return self.gene_class_df.loc[gene_id, 'occurence_in_positive_and_negative_samples'] == len(self.positive_samples) + len(self.negative_samples)


    def init_gene_class_df(self):
        """ generate dictionary with the class information per gene.

            This dictionary could be later use to produce an additional-layer
            text file for vizualization.
        """

        # TODO: make sure gene presence absence was calculated
        if not self.gene_coverage_consistency_dict_initiated:
            self.init_gene_coverage_consistency_information()
        # XXX: only negative and positive samples are used here
        # ambiguous samples are ignored as if they were never
        # there. This is not ideal, but is easy to do.
        gene_ids = self.gene_level_coverage_stats_dict_of_dataframes['mean_coverage'].index
        self.gene_class_df = pd.DataFrame(index=gene_ids)
        for gene_id in gene_ids:
            # determine the number of occurences in positive samples
            self.gene_class_df.loc[gene_id, 'occurence_in_positive_samples'] = len([s for s in self.positive_samples if self.gene_presence_absence_in_samples.loc[gene_id,s] == True])
            # determine the number of occurences in negative samples
            self.gene_class_df.loc[gene_id, 'occurence_in_negative_samples'] = len([s for s in self.negative_samples if self.gene_presence_absence_in_samples.loc[gene_id,s] == True])
            # set the occurence_in_positive_and_negative_samples
            self.gene_class_df.loc[gene_id, 'occurence_in_positive_and_negative_samples'] = self.gene_class_df.loc[gene_id, 'occurence_in_positive_samples'] + self.gene_class_df.loc[gene_id, 'occurence_in_negative_samples']

            gene_specificity = self.get_gene_specificity(gene_id)
            gene_coverage_consistency = self.get_gene_coverage_consistency(gene_id)
            # determine core accessory
            gene_is_core = self.determine_if_gene_is_core(gene_id, gene_specificity)

            self.gene_class_df.loc[gene_id, 'specificity'] = gene_specificity
            self.gene_class_df.loc[gene_id, 'coverage_consistency'] =gene_coverage_consistency
            self.gene_class_df.loc[gene_id, 'core'] = gene_is_core
            self.gene_class_df.loc[gene_id, 'MCG_class'] = get_class_string(gene_specificity, gene_coverage_consistency, gene_is_core)


    def update_samples_information_from_gene_class_df(self):
        # after running classification we sum up some information regarding
        # the results of the classifier per sample

        for sample in self.samples_detection_information:
            TSC = [g for g in self.gene_class_df.index if (self.gene_class_df.loc[g,'coverage_consistency'] and \
                                                            self.gene_class_df.loc[g,'core'])]
            self.samples_detection_information['number_of_taxon_specific_core_detected'] = len(TSC)


    def gen_gene_consistency_plots(self):
        """ generate and save the gene consistency plots for each gene."""

        if not self.gene_coverage_consistency_dict_initiated:
            self.init_gene_coverage_consistency_information()

        gene_ids = self.gene_level_coverage_stats_dict_of_dataframes['mean_coverage'].index
        num_genes, counter = len(gene_ids), 1
        progress.new('Plotting gene consistency information')
        progress.update('...')
        for gene_id in gene_ids:
            if num_genes > 100 and counter % 100 == 0:
                self.progress.update('%d of %d genes...' % (counter, num_genes))
            p = MCGPlots(self, gene_id, run=run, progress=progress)
            p.plot()

        progress.end()


    def save_gene_class_information_in_additional_layers(self):
        output_file_path = self.output_file_prefix + self.additional_description + '-additional-layers.txt'
        self.gene_class_df.to_csv(output_file_path, sep='\t', index_label='gene_callers_id')


    def save_samples_information(self):
        samples_information_file_name = self.output_file_prefix + self.additional_description + '-samples-information.txt'
        samples_information = pd.concat([self.samples_detection_information, self.samples_coverage_stats_dicts], axis=1, sort=True)
        samples_information.to_csv(samples_information_file_name, sep='\t', index_label='samples')


    def classify(self):
        self.init_gene_class_df()

        self.update_samples_information_from_gene_class_df()

        if self.write_output_to_files:
            self.save_gene_class_information_in_additional_layers()
            self.save_samples_information()

        if self.gen_figures:
            # Create the plots for nucleotide-level coverage data per sample.
            self.plot_nucleotide_coverage_distribution()
            # generate plots for coverage consistency information for each gene.
            self.gen_gene_consistency_plots()


def get_coverage_values_per_nucleotide(split_coverage_values_per_nt_dict, samples=None):
    """ Helper function that accepts a split_coverage_values_per_nt_dict and returns a dictionary with
    samples as keys and the concatenated coverage values for all splits as one array
    """

    if not split_coverage_values_per_nt_dict:
        raise ConfigError("You did not provide a split_coverage_values_per_nt_dict, and we need it...")

    progress.new('Merging coverage values accross splits')
    progress.update('...')

    d = {}
    if samples is None:
        samples = next(iter(split_coverage_values_per_nt_dict.values())).keys()

    number_of_samples = len(samples)
    number_of_finished = 0

    # find the combined legnth of all contigs first
    total_length = 0
    for split in split_coverage_values_per_nt_dict:
        total_length += len(split_coverage_values_per_nt_dict[split][next(iter(samples))])

    for sample in samples:
        # create an array of zero with the total length
        # this is much faster than appending the vectors of splits
        d[sample] = np.zeros(total_length)
        pos = 0
        for split in split_coverage_values_per_nt_dict:
            split_values = split_coverage_values_per_nt_dict[split][sample]
            split_len = len(split_values)
            d[sample][pos:pos+split_len] = split_values
            pos += split_len
        #d[sample] = np.array(d[sample])
        number_of_finished += 1
        progress.update("Finished sample %d of %d" % (number_of_finished,number_of_samples))

    progress.end()

    return d


def get_non_outliers_information(v, MAD_threshold=2.5, zeros_are_outliers=False):
    """ returns the non-outliers for the input pandas series using MAD"""

    d = pd.Series(index=columns_for_samples_coverage_stats_dict)
    outliers = utils.get_list_of_outliers(v, threshold=MAD_threshold, zeros_are_outliers=zeros_are_outliers)
    non_outliers = np.logical_not(outliers)
    non_outlier_indices = np.where(non_outliers)[0]

    if not(len(non_outlier_indices)):
        non_outlier_indices = np.array([])
        d['non_outlier_mean_coverage'] = 0.0
        d['non_outlier_coverage_std'] = 0.0

    else:
        d['non_outlier_mean_coverage'] = np.mean(v[non_outlier_indices])
        d['non_outlier_coverage_std'] = np.std(v[non_outlier_indices])

    return non_outlier_indices, d

# The order of the strings is very important since it is used in get_class_string
class_short_names = ['NNA', 'SNA', 'NCA',\
                     'SCA', 'NNC', 'SNC',\
                     'NCC', 'SCC']
class_long_names = ['Non-specific_Non-consistent_Accessory', 'Specific_Non-consistent_Accessory', 'Non-specific_Consistent_Accessory',\
                    'Specific_Consistent_Accessory', 'Non-specific_Non-consistent_Core', 'Specific_Non-consistent_Core',\
                    'Non-specific_Consistent_Core', 'Specific_Consistent_Core']
class_short_name_long_name_dict = dict(zip(class_short_names,class_long_names))


def get_class_long_name_from_short_name(short_name):
    return class_short_name_long_name_dict[short_name]


def get_class_string(gene_specificity, gene_coverage_consistency, gene_is_core):
    """ Takes the values of the three categories and returns a string to represent the class."""
    value_list = [gene_specificity, gene_coverage_consistency, gene_is_core]
    if None in value_list:
        return 'NA'
    # converting the list of booleans to a number
    # this solution was takes from here: https://stackoverflow.com/a/4066807/7115450
    index = sum(1<<i for i, b in enumerate(value_list) if b)
    return class_short_names[index]


def get_presence_absence_information(number_of_non_outliers, alpha):
    """ Helper function to determine presence/absence according to a threshold."""
    ##### WHAT WE SHOULD DO IN THE FUTURE #####
    # Arbitrary cut-offs are terrible.
    # If we assume there are no accessory genes (we will get back to this later),
    # then if the gnomes is present, then we expect ALL of it to be present. Thus,
    # if we had an unlimited number of reads, then we expect detection to be 1.
    # as the number of reads gets smaller, the expected detection value is smaller.
    # for a given genome size, a given read length, and the number of reads mapped to
    # the genome, we can compute the following value: "what is the probability that
    # the detection value will be greater than the actual detection value", if that
    # probability is high, then that is a good sign that the genome is not present
    # in the sample, and that any reads that we got are due to non-specific coverage.
    # the same thing could be calculated for a given gene.
    # we can create a measure for agreement between the mean coverage of a gene
    # and the detection of the gene. It would simply be the probability that the
    # coverage of the gene would exist with a detection that is higher than the
    # actual detection of the gene. All we need for that is the read length,
    # gene/genome length, and the expected genomic portion shared by two genomes that
    # belong to the population in question.
    if number_of_non_outliers >= 0.5 + alpha:
        return True
    elif np.sum(number_of_non_outliers) <= 0.5 - alpha:
        return False
    else:
        return None

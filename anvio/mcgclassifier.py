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
import matplotlib.pyplot as plt
import anvio.terminal as terminal
import anvio.summarizer as summarizer
import anvio.filesnpaths as filesnpaths

from scipy import odr as odr
from anvio.mcgops import MCGPlots
from anvio.errors import ConfigError
from anvio.dbops import ProfileSuperclass
from anvio.sequence import get_list_of_outliers
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
        self.gen_figures = A('gen_figures')
        self.split_coverage_values_per_nt_dict = {}
        self.gene_level_coverage_stats_dict = {}
        self.profile_db = {}
        self.coverage_values_per_nt = {}
        self.gene_coverages = {}
        self.gene_detections = None
        self.gene_coverage_values_per_nt = {}
        self.gene_non_outlier_coverage_std = None
        self.gene_non_outlier_positions = {}
        self.samples = {}
        self.sample_detection_information_was_initiated = False
        self.positive_samples = []
        self.number_of_positive_samples = None
        self.negative_samples = {}
        self.number_of_negative_samples = None
        self.gene_class_df = {}
        self.samples_detection_information = {}
        self.gene_presence_absence_in_samples_initiated = False
        self.gene_presence_absence_in_samples = {}
        self.gene_coverages_filtered = {}
        self.additional_description = ''
        self.total_length = None
        self.samples_coverage_stats_dicts_was_initiated = False
        self.samples_coverage_stats_dicts = {}
        self.non_outlier_indices = {}
        self.gene_coverage_consistency_dict = {}
        self.gene_coverage_consistency_dict_initiated = False
        self.samples_to_exclude = set([])
        self.samples_to_include = set([])

        if self.exclude_samples:
            # check that there is a file like this
            filesnpaths.is_file_exists(self.exclude_samples)
            self.samples_to_exclude = set([l.split('\t')[0].strip() for l in open(self.exclude_samples, 'rU').readlines()])

            if not self.samples_to_exclude:
                raise ConfigError("You asked to exclude samples, but provided an empty list.")

            run.info('Excluding Samples', 'The following samples will be excluded: %s' % self.samples_to_exclude,)

        if self.include_samples:
            # check that there is a file like this
            filesnpaths.is_file_exists(self.include_samples)
            self.samples_to_include = set([l.split('\t')[0].strip() for l in open(self.include_samples, 'rU').readlines()])

            if not self.samples_to_include:
                raise ConfigError("You provided an empty list of samples to include.")

            run.info('Including Samples', 'The following samples will be included: %s' % self.samples_to_include,)

        # run sanity check on all input arguments
        self.sanity_check()


    def init(self, gene_level_coverage_stats_dict=None, split_coverage_values_per_nt_dict=None, additional_description=None):
        """ setting the dictionaries for gene coverage stats and for split coverage per nucleotide"""
        self.gene_level_coverage_stats_dict = gene_level_coverage_stats_dict
        self.split_coverage_values_per_nt_dict = split_coverage_values_per_nt_dict
        self.additional_description = additional_description


    def get_gene_coverages_and_gene_detection_dicts(self):
        gene_coverages = {}
        gene_detection = {}
        gene_non_outlier_mean_coverage = {}
        gene_non_outlier_coverage_stds = {}

        A = lambda x: self.profile_db.gene_level_coverage_stats_dict[gene_callers_id][sample_name][x]

        gene_caller_ids = list(self.profile_db.gene_level_coverage_stats_dict.keys())

        # populate gene coverage and detection dictionaries
        if self.profile_db.gene_level_coverage_stats_dict:
            for gene_callers_id in gene_caller_ids:
                gene_coverages[gene_callers_id], gene_non_outlier_coverage_stds[gene_callers_id], gene_non_outlier_mean_coverage[gene_callers_id], gene_detection[gene_callers_id] = {}, {}, {}, {}
                self.gene_coverage_values_per_nt[gene_callers_id], self.gene_non_outlier_positions[gene_callers_id] = {}, {}

                for sample_name in self.profile_db.p_meta['samples']:
                    gene_coverages[gene_callers_id][sample_name] = A('mean_coverage')
                    gene_detection[gene_callers_id][sample_name] = A('detection')
                    gene_non_outlier_mean_coverage[gene_callers_id][sample_name] = A('non_outlier_mean_coverage')
                    gene_non_outlier_coverage_stds[gene_callers_id][sample_name] = A('non_outlier_coverage_std')
                    self.gene_coverage_values_per_nt[gene_callers_id][sample_name] = A('gene_coverage_values_per_nt')
                    self.gene_non_outlier_positions[gene_callers_id][sample_name] = A('non_outlier_positions')
                    

        return gene_coverages, gene_detection, gene_non_outlier_mean_coverage, gene_non_outlier_coverage_stds


    def check_if_valid_portion_value(self, arg_name,arg_value):
        """ Helper function to verify that an argument has a valid value for a non-zero portion (i.e. greater than zero and a max of 1)"""
        if arg_value <= 0 or arg_value > 1:
            raise ConfigError("%s value must be greater than zero and a max of 1, the value you supplied %s" % (arg_name,arg_value))


    def sanity_check(self):
        """Basic sanity check for class inputs"""

        if self.output_file_prefix:
            # checking output file
            filesnpaths.is_output_file_writable(self.output_file_prefix + '-additional-layers.txt', ok_if_exists=False)

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
        samples = set(samples_list) - self.samples_to_exclude
        if self.include_samples:
            samples_to_include_that_are_not_there = self.samples_to_include - samples
            if samples_to_include_that_are_not_there:
                raise ConfigError("You requested to include some samples that are not in the profile database. Here are the samples in the profile database: %s. \
                                And here are the samples you requested, and that are not there: %s" % (samples, samples_to_include_that_are_not_there))
            samples = self.samples_to_include
        self.samples = samples


    def init_coverage_and_detection_dataframes(self, gene_coverages_dict, gene_detection_dict, gene_non_outlier_coverages_dict, gene_non_outlier_coverage_stds_dict):
        """ Populate the following: self.gene_coverages, self.Ng, self.gene_detections.

            Notice that this function could get as input either an object of ProfileSuperclass or of summarizer.Bin
        """
        self.gene_coverages = pd.DataFrame.from_dict(gene_coverages_dict, orient='index', dtype=float)
        self.Ng = len(self.gene_coverages.index)
        self.gene_detections = pd.DataFrame.from_dict(gene_detection_dict, orient='index', dtype=float)
        self.gene_non_outlier_coverages = pd.DataFrame.from_dict(gene_non_outlier_coverages_dict, orient='index', dtype=float)
        self.gene_non_outlier_coverage_stds = pd.DataFrame.from_dict(gene_non_outlier_coverage_stds_dict, orient='index', dtype=float)


        if self.include_samples or self.exclude_samples:
            # Only include samples that the user want
            self.gene_coverages = self.gene_coverages[list(self.samples)]
            self.gene_detections = self.gene_detections[list(self.samples)]
            self.gene_non_outlier_coverages = self.gene_non_outlier_coverages[list(self.samples)]
            self.gene_non_outlier_coverage_stds = self.gene_non_outlier_coverage_stds[list(self.samples)]


    def init_sample_detection_information(self):
        """ Determine  positive, negative, and ambiguous samples with the genome detection information
        (--alpha, --genome-detection-uncertainty)

            The coverage_values_per_nt is used to calculate the detection value (portion of nucleotides
            covered) for a sample. Then, a cutoff for detection values is used to determine the presence
            or absence of the genome in each sample.
        """

        # FIXME: some of the following variables are never used.
        MCG_samples_information_table_name      = 'MCG_classifier_samples_information'
        MCG_samples_information_table_structure = ['samples', 'presence', 'detection', 'number_of_taxon_specific_core_detected']
        MCG_samples_information_table_types     = ['str', 'bool', 'int', 'int']

        # create an empty dataframe
        samples_information = pd.DataFrame(index=self.samples, columns=MCG_samples_information_table_structure[1:])
        positive_samples = []
        negative_samples = []

        self.progress.new("Setting presence/absence in samples")
        progress.update('...')
        num_samples, counter = len(self.samples), 1
        detection = {}
        for sample in self.samples:
            # TODO: if coverage_values_per_nt is not available then use gene_detections to estimate the detection value of a sample
            if num_samples > 100 and counter % 100 == 0:
                self.progress.update('%d of %d samples...' % (counter, num_samples))
            print("total length for %s is %s" % (sample, self.total_length))
            print("the length of the vector: %s" % len(self.coverage_values_per_nt[sample])) # FIXME: after testing this module, delete this line. it is only here to make sure that anvio is not lying to us.
            print("number of nucleotide positions with non zero coverage in %s is %s " % (sample, np.count_nonzero(self.coverage_values_per_nt[sample])))
            detection[sample] = np.count_nonzero(self.coverage_values_per_nt[sample]) / self.total_length
            samples_information['presence'][sample] = get_presence_absence_information(detection[sample], self.alpha)
            if samples_information['presence'][sample]:
                positive_samples.append(sample)
            elif samples_information['presence'][sample] == False:
                negative_samples.append(sample)

            samples_information['detection'][sample] = detection[sample]
            counter += 1
        self.progress.end()

        self.positive_samples = positive_samples
        self.number_of_positive_samples = len(self.positive_samples)
        self.negative_samples = negative_samples
        self.samples_detection_information = samples_information
        self.run.warning('The number of positive samples is %s' % self.number_of_positive_samples)
        self.run.warning('The number of negative samples is %s' % len(self.negative_samples))
        self.sample_detection_information_was_initiated = True


    def init_samples_coverage_stats_dict(self):
        """ populate the samples_coverage_stats_dict.
        
            This dataframe is used to calculate the gene consistency information.
            It is also used for plotting purposes (both for the nucleotide-coverage-distribution plots and the gene-consistency plots).
        """
        # TODO: move this to wrapper class
        if not self.sample_detection_information_was_initiated:
            self.init_sample_detection_information()

        self.samples_coverage_stats_dicts = pd.DataFrame(index=self.samples, columns=columns_for_samples_coverage_stats_dict)

        num_samples, counter = len(self.samples), 1
        self.progress.new("Finding nucleotide positions in samples with outlier coverage values")
        progress.update('...')
        for sample in self.positive_samples:
            if num_samples > 100 and counter % 100 == 0:
                self.progress.update('%d of %d samples...' % (counter, num_samples))

            # loop through positive samples
            # get the non-outlier information
            non_outlier_indices, self.samples_coverage_stats_dicts.loc[sample,] = get_non_outliers_information(self.coverage_values_per_nt[sample], MAD_threshold=self.outliers_threshold)
            self.non_outlier_indices[sample] = non_outlier_indices
            # TODO: in manual mode this will either be supplied or it will be calculated from gene coverages

            self.run.info_single('The mean and std of non-outliers in sample %s are: %s, %s respectively' % (sample, self.samples_coverage_stats_dicts['non_outlier_mean_coverage'][sample], self.samples_coverage_stats_dicts['non_outlier_coverage_std'][sample]))
            number_of_non_outliers = len(self.non_outlier_indices[sample])
            self.run.info_single('The number of non-outliers is %s of %s (%.2f%%)' % (number_of_non_outliers, self.total_length, 100.0 * number_of_non_outliers / self.total_length))
        self.samples_coverage_stats_dicts_was_initiated = True
        self.progress.end()


    def store_samples_coverage_stats_dict(self):
        """ sotre samples_coverage_stats_dict into TAB-delimited file"""

        if not self.samples_coverage_stats_dicts_was_initiated:
            raise ConfigError("can't dtore samples_coverage_stats_dict, because it wasn't initiated")

        output_file_path = self.output_file_prefix + '-samples-coverage-stats.txt'
        self.samples_coverage_stats_dicts.to_csv(output_file_path, sep='\t', index_label='sample')


    def plot_nucleotide_coverage_distribution(self):
        """ Creates a pdf file with the following plots for each sample the sorted nucleotide coverages \
        (with the outliers in red and non-outliers in blue), and a histogram of coverages for the non-outliers"""
        # Creating a dircetory for the plots. If running on bins, each bin would be in a separate sub-directory

        if not self.samples_coverage_stats_dicts_was_initiated:
            self.init_samples_coverage_stats_dict()

        additional_description = ''
        if self.additional_description:
            additional_description = '-' + self.additional_description

        plot_dir = self.output_file_prefix + '-nucleotide-coverage-distribution-plots' + '/'
        os.makedirs(plot_dir, exist_ok=True)
        self.progress.new('Saving figures of taxon specific distributions to pdf')
        progress.update('...')
        number_of_fininshed = 0
        for sample in self.positive_samples:
            coverages_pdf_output = plot_dir + sample + additional_description + '-coverages.pdf'
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
            y2 = v[self.non_outlier_indices[sample]]
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
        if self.gene_detections is None:
            # making sure that a gene detections table exists
            raise ConfigError("gene presence/absence in samples cannot be determined without a gene detections table,\
                                but it seems that you don't have a gene detections table.")

        gene_callers_id = self.gene_detections.index
        self.gene_presence_absence_in_samples = pd.DataFrame(index=gene_callers_id, columns=self.samples)

        num_samples, counter = len(self.samples), 1
        self.progress.new('Computing gene presence/absence in samples')
        progress.update('...')
        for sample in self.samples:
            if num_samples > 100 and counter % 100 == 0:
                self.progress.update('%d of %d samples...' % (counter, num_samples))
            for gene_id in gene_callers_id:
                self.gene_presence_absence_in_samples.loc[gene_id, sample] = get_presence_absence_information(self.gene_detections.loc[gene_id, sample], self.alpha)
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
        num_genes, counter = len(self.gene_coverages.index), 1
        for gene_id in self.gene_coverages.index:
            if num_genes > 100 and counter % 100 == 0:
                self.progress.update('%d of %d genes...' % (counter, num_genes))

            # samples in which the gene is present
            _samples = self.gene_presence_absence_in_samples.loc[gene_id,self.gene_presence_absence_in_samples.loc[gene_id,]==True].index
            # mean and std of non-outlier nt in each sample
            x = self.samples_coverage_stats_dicts.loc[_samples,'non_outlier_mean_coverage']
            if "non_outlier_coverage_std" in self.samples_coverage_stats_dicts:
                # we only expect to have the sample coverage std in "full" mode
                std_x = self.samples_coverage_stats_dicts.loc[_samples,'non_outlier_coverage_std'].values
            else:
                std_x = None

            if len(_samples) > 1:
                # mean and std of non-outlier nt in the gene (in each sample)
                y = self.gene_non_outlier_coverages.loc[gene_id, _samples].values
                if self.gene_non_outlier_coverage_stds is None:
                    # checking if 
                    std_y = None
                else:
                    std_y = self.gene_non_outlier_coverage_stds.loc[gene_id, _samples].values

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
            return self.gene_coverage_consistency_dict[gene_id]['slope_precision'] > 0.5
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
            # return True if the the gene occurs in all positive samples.
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
        self.gene_class_df = pd.DataFrame(index=list(self.gene_coverages.index))
        for gene_id in self.gene_coverages.index:
            # determine the number of occurences in positive samples
            self.gene_class_df.loc[gene_id, 'occurence_in_positive_samples'] = np.sum(self.gene_presence_absence_in_samples.loc[gene_id, self.positive_samples])
            # determine the number of occurences in negative samples
            self.gene_class_df.loc[gene_id, 'occurence_in_negative_samples'] = np.sum(self.gene_presence_absence_in_samples.loc[gene_id, self.negative_samples])
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


    def get_gene_classes(self):
        """ The main process of this class - computes the class information for each gene"""
        # copmute coverage stats values
        self.init_samples_coverage_stats_dict()

        # find occurence of genes in the samples
        self.init_gene_presence_absence_in_samples()

        # compute gene consistency information
        self.init_gene_coverage_consistency_information()

        if self.gen_figures:
            # Create the plots for nucleotide-level coverage data per sample.
            self.plot_nucleotide_coverage_distribution()
            # generate plots for coverage consistency information for each gene.
            self.gen_gene_consistency_plots() 

        # create the gene_class_df
        self.init_gene_class_df()


    def gen_gene_consistency_plots(self):
        """ generate and save the gene consistency plots for each gene."""

        if not self.gene_coverage_consistency_dict_initiated:
            self.init_gene_coverage_consistency_information()

        num_genes, counter = len(self.gene_coverages.index), 1
        progress.new('Plotting gene consistency information')
        progress.update('...')
        for gene_id in self.gene_coverages.index:
            if num_genes > 100 and counter % 100 == 0:
                self.progress.update('%d of %d genes...' % (counter, num_genes))
            p = MCGPlots(self, gene_id, run=run, progress=progress)
            p.plot()

        progress.end()


    def get_coverage_and_detection_dict(self,bin_id):
        _bin = summarizer.Bin(self.summary, bin_id)
        self.coverage_values_per_nt = get_coverage_values_per_nucleotide(_bin.split_coverage_values_per_nt_dict, self.samples)


        # getting the total length of all contigs
        self.total_length = _bin.total_length

        self.init_coverage_and_detection_dataframes(_bin.gene_coverages, _bin.gene_detection, _bin.gene_non_outlier_coverages, _bin.gene_non_outlier_coverage_stds)

        # populate gene per-position information
        self.gene_coverage_values_per_nt = _bin.gene_coverage_values_per_nt
        self.gene_non_outlier_positions = _bin.gene_non_outlier_positions

    def save_gene_class_information_in_additional_layers(self, additional_description=''):
        additional_column_titles = []
        additional_layers_df = self.gene_class_df
        output_file_path = self.output_file_prefix + '-additional-layers.txt'
        additional_layers_df.to_csv(output_file_path, sep='\t', index_label='gene_callers_id')


    def save_samples_information(self, additional_description=''):
        # TODO: there used to be this here:
        #self.run.warning(self.samples_information)
        samples_information_df = self.samples_information

        if additional_description:
            additional_description = '-' + additional_description

        samples_information_file_name = self.output_file_prefix + additional_description + '-samples-information.txt'
        samples_information_df.to_csv(samples_information_file_name, sep='\t', index_label='samples')

    def classify(self):
        
        # FIXME:this is for debugging and should be removed
        print(self.gene_level_coverage_stats_dict)
        print(self.split_coverage_values_per_nt_dict)

        return
        self.get_gene_classes()
        #self.save_gene_class_information_in_additional_layers(bin_id)
        #self.save_samples_information(bin_id)


def get_coverage_values_per_nucleotide(split_coverage_values_per_nt_dict, samples=None):
    """ Helper function that accepts a split_coverage_values_per_nt_dict and returns a dictionary with
    samples as keys and the concatenated coverage values for all splits as one array
    """
    progress.new('Merging coverage values accross splits')
    progress.update('...')

    d = {}
    if samples is None:
        samples = split_coverage_values_per_nt_dict[next(iter(split_coverage_values_per_nt_dict.keys()))].keys()

    number_of_samples = len(samples)
    number_of_finished = 0

    for sample in samples:
        d[sample] = np.empty([1,0])
        for split in split_coverage_values_per_nt_dict:
            d[sample] = np.append(d[sample],split_coverage_values_per_nt_dict[split][sample])
        #d[sample] = np.array(d[sample])
        number_of_finished += 1
        progress.update("Finished sample %d of %d" % (number_of_finished,number_of_samples))

    progress.end()

    return d


def get_non_outliers_information(v, MAD_threshold=2.5):
    """ returns the non-outliers for the input pandas series using MAD"""

    d = pd.Series(index=columns_for_samples_coverage_stats_dict)
    outliers = get_list_of_outliers(v, threshold=MAD_threshold)
    non_outliers = np.logical_not(outliers)

    if not(len(non_outliers)):
        non_outlier_indices = np.array([])
        d['non_outlier_mean_coverage'] = 0.0
        d['non_outlier_coverage_std'] = 0.0

    else:
        non_outlier_indices = np.where(non_outliers)[0]
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


def get_presence_absence_information(detection, alpha):
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
    if detection >= 0.5 + alpha:
        return True
    elif detection <= 0.5 - alpha:
        return False
    else:
        return None

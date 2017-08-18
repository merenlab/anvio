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
import matplotlib.pyplot as plt
import anvio.terminal as terminal
import anvio.summarizer as summarizer
import anvio.filesnpaths as filesnpaths

from anvio.errors import ConfigError
from anvio.dbops import ProfileSuperclass
from matplotlib.backends.backend_pdf import PdfPages

matplotlib.use('pdf')

__author__ = "Alon Shaiber"
__copyright__ = "Copyright 2017, The anvio Project"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "Alon Shaiber"
__email__ = "alon.shaiber@gmail.com"


run = terminal.Run()
progress = terminal.Progress()
pp = terminal.pretty_print


class MetagenomeCentricGeneClassifier:
    def __init__(self, args, run=run, progress=progress):
        self.run = run
        self.progress = progress

        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        self.gene_coverages_data_file_path = A('data_file')
        self.gene_detections_data_file_path = A('gene_detection_data_file')
        self.profile_db_path = A('profile_db')
        self.output_file_prefix = A('output_file_prefix')
        self.alpha = A('alpha')
        self.beta = A('beta')     #
        self.gamma = A('gamma')   # FIXME: beta, gamma, eta, and zeta are not 
        self.eta = A('eta')       #        used anywhere in the code anyore :)
        self.zeta = A('zeta')     #
        self.additional_layers_to_append = A('additional_layers_to_append')
        self.samples_information_to_append = A('samples_information_to_append')
        self.collection_name = A('collection_name')
        self.bin_id = A('bin_id')
        self.bin_ids_file_path = A('bin_ids_file')
        self.exclude_samples = A('exclude_samples')
        self.include_samples = A('include_samples')
        self.profile_db = {}
        self.coverage_values_per_nt = {}
        self.gene_coverages = pd.DataFrame.empty
        self.gene_detections = pd.DataFrame.empty
        self.samples = {}
        self.positive_samples = []
        self.number_of_positive_samples = None
        self.negative_samples = pd.DataFrame.empty
        self.number_of_negative_samples = None
        self.gene_class_information = pd.DataFrame.empty
        self.samples_information = pd.DataFrame.empty
        self.gene_presence_absence_in_samples = pd.DataFrame.empty
        self.gene_coverages_filtered = pd.DataFrame.empty
        self.additional_description = ''
        self.total_length = None

        if self.exclude_samples:
            # check that there is a file like this
            filesnpaths.is_file_exists(self.exclude_samples)
            self.samples_to_exclude = set([l.split('\t')[0].strip() for l in open(self.exclude_samples, 'rU').readlines()])

            if not self.samples_to_exclude:
                raise ConfigError("You asked to exclude samples, but provided an empty list.")

            run.info('Excluding Samples', 'The following samples will be excluded: %s' % self.samples_to_exclude,)
        else:
            self.samples_to_exclude = set([])

        if self.include_samples:
            # check that there is a file like this
            filesnpaths.is_file_exists(self.include_samples)
            self.samples_to_include = set([l.split('\t')[0].strip() for l in open(self.include_samples, 'rU').readlines()])

            if not self.samples_to_include:
                raise ConfigError("You provided an empty list of samples to include.")

            run.info('Including Samples', 'The following samples will be included: %s' % self.samples_to_include,)
        else:
            self.samples_to_include = set([])

        # run sanity check on all input arguments
        self.sanity_check()

        if self.profile_db_path is None:
            # TODO: this will probably be removed because we don't save the coverage information in nucleotide level.
            pass
        else:
            # load sample list and gene_coverage_dict from the merged profile db
            args.init_gene_coverages = True
            if self.collection_name:
                self.summary = summarizer.ProfileSummarizer(args)
                self.summary.init()
                self.init_samples(self.summary.p_meta['samples'])
            else:
                self.profile_db = ProfileSuperclass(args)
                self.init_samples(self.profile_db.p_meta['samples'])
                self.profile_db.init_split_coverage_values_per_nt_dict()
                self.profile_db.init_gene_level_coverage_stats_dicts()
                self.coverage_values_per_nt = get_coverage_values_per_nucleotide(self.profile_db.split_coverage_values_per_nt_dict, self.samples)

                # comply with the new design and get gene_coverages and gene_detection dicsts from
                # gene_level_coverage_stats_dict.
                gene_coverages, gene_detection = self.get_gene_coverages_and_gene_detection_dicts()

                self.init_coverage_and_detection_dataframes(gene_coverages, gene_detection)

                # getting the total length of all contigs
                self.total_length = self.profile_db.p_meta['total_length']


    def get_gene_coverages_and_gene_detection_dicts(self):
        gene_coverages = {}
        gene_detection = {}

        A = lambda x: self.profile_db.gene_level_coverage_stats_dict[gene_callers_id][sample_name][x]

        gene_caller_ids = list(self.profile_db.gene_level_coverage_stats_dict.keys())

        # populate gene coverage and detection dictionaries
        if self.profile_db.gene_level_coverage_stats_dict:
            for gene_callers_id in gene_caller_ids:
                gene_coverages[gene_callers_id], gene_detection[gene_callers_id] = {}, {}

                for sample_name in self.profile_db.p_meta['samples']:
                    gene_coverages[gene_callers_id][sample_name] = A('mean_coverage')
                    gene_detection[gene_callers_id][sample_name] = A('detection')

        return gene_coverages, gene_detection


    def check_if_valid_portion_value(self, arg_name,arg_value):
        """ Helper function to verify that an argument has a valid value for a non-zero portion (i.e. greater than zero and a max of 1)"""
        if arg_value <= 0 or arg_value > 1:
            raise ConfigError("%s value must be greater than zero and a max of 1, the value you supplied %s" % (arg_name,arg_value))


    def sanity_check(self):
        """Basic sanity check for class inputs"""

        if self.profile_db_path is None and self.gene_coverages_data_file_path is None:
            raise ConfigError("You must provide either a profile.db or a gene coverage self.gene_coverages_filtered data file")

        if self.profile_db_path and self.gene_coverages_data_file_path:
            raise ConfigError("You provided both a profile database and a gene coverage self.gene_coverages_filtered data file, you \
            must provide only one or the other (hint: if you have a profile database, the use that")

        # checking output file
        filesnpaths.is_output_file_writable(self.output_file_prefix + '-additional-layers.txt', ok_if_exists=False)

        # checking alpha
        if not isinstance(self.alpha, float):
            raise ConfigError("alpha value must be a type float.")
        # alpha must be a min of 0 and smaller than 0.5
        if self.alpha < 0 or self.alpha >= 0.5:
            raise ConfigError("alpha must be a minimum of 0 and smaller than 0.5")

        # Checking beta
        if not isinstance(self.beta, float):
            raise ConfigError("beta value must be a type float.")
        self.check_if_valid_portion_value("beta", self.beta)

        # Checking gamma
        if not isinstance(self.gamma, float):
            raise ConfigError("Gamma value must be a type float.")
        self.check_if_valid_portion_value("gamma", self.gamma)

        # Checking eta
        self.check_if_valid_portion_value("eta", self.eta)

        if self.collection_name:
            if not self.profile_db_path:
                raise ConfigError("You specified a collection name %s, but you provided a gene coverage self.gene_coverages_filtered data file \
                 collections are only available when working with a profile database." % self.collection_name)

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


    def init_coverage_and_detection_dataframes(self, gene_coverages_dict, gene_detection_dict):
        """ Populate the following: self.gene_coverages, self.Ng, self.gene_detections.

            Notice that this function could get as input either an object of ProfileSuperclass or of summarizer.Bin
        """
        self.gene_coverages = pd.DataFrame.from_dict(gene_coverages_dict, orient='index', dtype=float)
        self.Ng = len(self.gene_coverages.index)
        self.gene_detections = pd.DataFrame.from_dict(gene_detection_dict, orient='index', dtype=float)

        if self.include_samples or self.exclude_samples:
            # Only include samples that the user want
            self.gene_coverages = self.gene_coverages[list(self.samples)]
            self.gene_detections = self.gene_detections[list(self.samples)]


    def init_sample_detection_information(self):
        """ Determine  positive, negative, and ambiguous samples with the genome detection information
        (--alpha, --genome-detection-uncertainty)
        """

        # FIXME: some of the following variables are never used.
        MDG_samples_information_table_name      = 'MDG_classifier_samples_information'
        MDG_samples_information_table_structure = ['samples', 'presence', 'detection', 'number_of_taxon_specific_core_detected']
        MDG_samples_information_table_types     = ['str', 'bool', 'int', 'int']

        # create an empty dataframe
        samples_information = pd.DataFrame(index=self.samples, columns=MDG_samples_information_table_structure[1:])
        positive_samples = []
        negative_samples = []

        self.progress.new("Setting presence/absence in samples")
        num_samples, counter = len(self.samples), 1
        detection = {}
        for sample in self.samples:
            if num_samples > 100 and counter % 100 == 0:
                self.progress.update('%d of %d samples...' % (counter, num_samples))
            print("total length for %s is %s" % (sample, self.total_length))
            print("the length of the vector: %s" % len(self.coverage_values_per_nt[sample]))
            print("number of non zero in %s is %s " % (sample, np.count_nonzero(self.coverage_values_per_nt[sample])))
            detection[sample] = np.count_nonzero(self.coverage_values_per_nt[sample]) / self.total_length
            if detection[sample] >= 0.5 + self.alpha:
                positive_samples.append(sample)
                samples_information['presence'][sample] = True
            elif detection[sample] <= 0.5 - self.alpha:
                negative_samples.append(sample)
                samples_information['presence'][sample] = False
            else:
                samples_information['presence'][sample] = None
            samples_information['detection'][sample] = detection[sample]
            counter += 1
        self.progress.end()

        self.positive_samples = positive_samples
        self.number_of_positive_samples = len(self.positive_samples)
        self.negative_samples = negative_samples
        self.samples_information = samples_information
        self.run.warning('The number of positive samples is %s' % self.number_of_positive_samples)
        self.run.warning('The number of negative samples is %s' % len(self.negative_samples))


    def plot_TS(self, non_outliers_indices, mean_TS, std_TS):
        """ Creates a pdf file with the following plots for each sample the sorted nucleotide coverages \
        (with a the outliers in red and non-outliers in blue), and a histogram of coverages for the non-outliers"""
        # Creating a dircetory for the plots. If running on bins, each bin would be in a separate sub-directory
        additional_description = ''
        if self.additional_description:
            additional_description = '-' + self.additional_description
        plot_dir = self.output_file_prefix + '-TS-plots' + '/'
        os.makedirs(plot_dir, exist_ok=True)
        self.progress.new('Saving figures of taxon specific distributions to pdf')
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
            x2 = reverse_sorted_indices[non_outliers_indices[sample]]
            y2 = v[non_outliers_indices[sample]]
            # plot all in red
            ax.semilogy(x1,v[sorting_indices],'r.', rasterized=True)
            # plot on top the non-outliers in blue
            ax.semilogy(x2,v[non_outliers_indices[sample]],'b.', rasterized=True)
            fig.suptitle("%s - sorted coverage values with outliers" % sample)
            plt.savefig(pdf_output_file, format='pdf')
            plt.close()

            # plotting a histogram of the non-outliers
            # This would allow to see if they resemble a normal distribution
            hist_range = (min(v[non_outliers_indices[sample]]),max(v[non_outliers_indices[sample]]))
            # computing the number of bins so that the width of a bin is ~1/4 of the standard deviation
            # FIXME: need to make it so the bins are only of integers (so the smallest bin is of width 1
            # and that bins are integers)
            number_of_hist_bins = np.ceil((hist_range[1] - hist_range[0]) / (std_TS[sample]/4)).astype(int) # setting the histogram bins to be of the width of a quarter of std
            fig = plt.figure()
            ax = fig.add_subplot(111, rasterized=True)
            ax.set_xlabel = 'Coverage'
            ax.hist(v[non_outliers_indices[sample]], number_of_hist_bins,hist_range, rasterized=True)
            fig.suptitle("%s - histogram of non-outliers" % sample)
            # adding the mean and std of the non-outliers as text to the plot
            text_for_hist = u'$\mu = %d$\n $\sigma = %d$' % (mean_TS[sample], std_TS[sample])
            ax.text(0.8, 0.9, text_for_hist, ha='center', va='center', transform=ax.transAxes)
            plt.savefig(pdf_output_file, format='pdf')
            plt.close()
            # close the pdf file
            pdf_output_file.close()
            number_of_fininshed += 1
            self.progress.update("Finished %d of %d" % (number_of_fininshed, self.number_of_positive_samples))
        self.progress.end()


    def get_taxon_specific_genes_in_samples(self):
        """ Use only positive samples to identify the single copy taxon specific genes in each sample:
        """
        non_outliers_indices = {}
        mean_TS = {}
        std_TS = {}
        num_samples, counter = len(self.samples), 1
        self.progress.new("Finding taxon specific genes in samples")
        for sample in self.positive_samples:
            if num_samples > 100 and counter % 100 == 0:
                self.progress.update('%d of %d samples...' % (counter, num_samples))

            # loop through positive samples
            # get the indexes of the non outliers and a pdf for the coverage of the single copy core genes
            non_outliers_indices[sample], mean_TS[sample], std_TS[sample] = get_non_outliers(self.coverage_values_per_nt[sample])

            self.run.info_single('The mean and std in sample %s are: %s, %s respectively' % (sample, mean_TS[sample], std_TS[sample]))
            self.run.info_single('The number of non_outliers is %s of %s' % (len(non_outliers_indices[sample]), self.total_length))
        self.progress.end()
        self.plot_TS(non_outliers_indices,mean_TS,std_TS)


    def get_gene_classes(self):
        """ The main process of this class - computes the class information for each gene"""
        # need to start a new gene_class_information dict
        # this is due to the fact that if the algorithm is ran on a list of bins then this necessary
        self.gene_class_information = pd.DataFrame(index=self.gene_coverages.index,columns=['gene_class'])

        # set the presence/absence values for samples
        self.init_sample_detection_information()

        # find the taxon-specific genes for each sample
        self.get_taxon_specific_genes_in_samples()


    def get_coverage_and_detection_dict(self,bin_id):
        _bin = summarizer.Bin(self.summary, bin_id)
        self.coverage_values_per_nt = get_coverage_values_per_nucleotide(_bin.split_coverage_values_per_nt_dict, self.samples)

        # getting the total length of all contigs
        self.total_length = _bin.total_length

        self.init_coverage_and_detection_dataframes(_bin.gene_coverages, _bin.gene_detection)


    def classify(self):
        if self.collection_name:
            bin_names_in_collection = self.summary.bin_ids
            if self.bin_ids_file_path:
                filesnpaths.is_file_exists(self.bin_ids_file_path)
                bin_names_of_interest = [line.strip() for line in open(self.bin_ids_file_path).readlines()]

                missing_bins = [b for b in bin_names_of_interest if b not in bin_names_in_collection]
                if len(missing_bins):
                    raise ConfigError("Some bin names you declared do not appear to be in the collection %s. \
                                        These are the bins that are missing: %s, these are the bins that are \
                                        actually in your collection: %s" % (self.collection_name,missing_bins,bin_names_in_collection))
            elif self.bin_id:
                if self.bin_id not in bin_names_in_collection:
                    raise ConfigError("The bin you declared, %s, does not appear to be in the collection %s." \
                                      % (self.bin_id, self.collection_name))
                bin_names_of_interest = [self.bin_id]
            else:
                bin_names_of_interest = bin_names_in_collection

            for bin_id in bin_names_of_interest:
                self.run.info_single('Classifying genes in bin: %s' % bin_id)
                self.get_coverage_and_detection_dict(bin_id)
                self.additional_description = bin_id
                self.get_gene_classes()
                #self.save_gene_class_information_in_additional_layers(bin_id)
                #self.save_samples_information(bin_id)
                #if self.store_gene_detections_and_gene_coverages_tables:
                #    self.save_gene_detection_and_coverage(bin_id)

        else:
            # No collection provided so running on the entire detection table
            self.get_gene_classes()
            #self.save_gene_class_information_in_additional_layers()
            #self.save_samples_information()


def get_coverage_values_per_nucleotide(split_coverage_values_per_nt_dict, samples=None):
    """ Helper function that accepts a split_coverage_values_per_nt_dict and returns a dictionary with
    samples as keys and the concatenated coverage values for all splits as one array
    """
    progress.new('Merging coverage values accross splits')

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


def get_non_outliers(v):
    """ returns the non-outliers according to the interqurtile range (IQR) for the input pandas series"""
    q1 = np.percentile(v, 25)
    q3 = np.percentile(v, 75)
    IQR = q3 - q1

    # The non-outliers are non-zero values that are in the IQR (positions that are zero are considered outliers
    # even if the IQR includes zero)
    non_outliers_indices = np.where((v >= q1 - 1.5 * IQR) & (v <= q3 + 1.5 * IQR) & (v > 0))[0]
    mean = np.mean(v[non_outliers_indices])
    std = np.std(v[non_outliers_indices])

    return non_outliers_indices, mean, std


def get_new_mean(_mean, x, N):
    """ Helper function to calculate a new mean after removing one data point."""
    new_mean = N/(N-1)*_mean - 1/(N-1)*x

    return new_mean

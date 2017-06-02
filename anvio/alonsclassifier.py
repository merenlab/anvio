# -*- coding: utf-8
# pylint: disable=line-too-long
"""
    Classes to classify genes based on coverages across metagenomes.

    anvi-alons-classifier is the default client using this module
"""


import anvio
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import anvio.terminal as terminal
import anvio.summarizer as summarizer
import anvio.filesnpaths as filesnpaths

from math import ceil
from math import floor
from anvio.errors import ConfigError
from anvio.dbops import ProfileSuperclass
from matplotlib.backends.backend_pdf import PdfPages

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

# Defining the structure of the sample information dictionary

def get_non_outliers(v):
    """ returns the interqurtile range (IQR) for the input pandas series"""
    q1 = np.percentile(v, 25)
    q3 = np.percentile(v, 75)
    IQR = q3 - q1
    non_outliers = (v >= q1 - 1.5 * IQR) & (v <= q3 + 1.5 * IQR) & v > 0
    return non_outliers


def plot_outliers(pdf_output_file, v, non_outliers, sample_id):
    """ Takes a PdfFile object (matplotlib.backends.backend_pdf), a vector of coverage, and a boolean vector of 
    the non-outlier indexes, and plots the non-outliers in blue and the outliers in red"""
    mu = np.mean(v[non_outliers])
    sigma  = np.std(v[non_outliers])

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_xlabel = 'Gene Number (ordered)'
    ax.set_ylabel = r'$Gene Coverage^2$'
    v_sqr = np.sqrt(v)
    min_y = min(v_sqr)
    range_y = max(v_sqr) - min_y
    print(sample_id)
    print(mu)
    print(sigma)
    ax.text(0.25 * max(v.index), min_y + 10**0.75*range_y, u'Mean = %d\n Standard deviation = %d' % (mu, sigma))
    # sorting the coverage values
    x_values = v.index
    v.sort_values(inplace=True)
    ax.semilogy(v.values,'r.')
    ax.semilogy(v.reset_index().index[non_outliers[v.index]],v[non_outliers].values,'b.')
    fig.suptitle("%s - sorted coverage values with outliers" % sample_id)
    plt.savefig(pdf_output_file, format='pdf')
    plt.close()

    # plotting a histogram of the non-outliers
    # This would allow to see if they resemble a normal distribution
    number_of_hist_bins = 100
    hist_range = (max(0,mu - 3 * sigma), mu + 3*sigma)
    plt.hist(v[non_outliers.values], number_of_hist_bins,hist_range)
    plt.title("%s - histogram of non-outliers" % sample_id)
    plt.savefig(pdf_output_file, format='pdf')
    plt.close()


def get_sliding_window(N, window_portion, overlap=None):
    """ Accepts an integer (number of items), a window size, and window overlap portion, and returns a list of tuples \
    where each item is a pair of idexes, a start and end index"""
    window = ceil(N * window_portion)
    if overlap is None:
        overlap = window * 0.5
    overlap_before = ceil(overlap)
    overlap_after = floor(overlap)
    window_indexes = [(0,window - 1)]
    max_index = window - 1
    while max_index + overlap_after < N - 1:
        window_indexes.append((max_index - overlap_before, max_index + overlap_after))
        max_index += overlap_after
        print(max_index)
    window_indexes.append((N - window - 1, N-1))
    return window_indexes


class AlonsClassifier:
    def __init__(self, args, run=run, progress=progress):
        self.run = run
        self.progress = progress

        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        self.gene_coverages_data_file_path = A('data_file')
        self.gene_detections_data_file_path = A('gene_detection_data_file')
        self.profile_db_path = A('profile_db')
        self.output_file_prefix = A('output_file_prefix')
        self.alpha = A('alpha')
        self.beta = A('beta')
        self.gamma = A('gamma')
        self.eta = A('eta')
        self.zeta = A('zeta')
        self.additional_layers_to_append = A('additional_layers_to_append')
        self.samples_information_to_append = A('samples_information_to_append')
        self.collection_name = A('collection_name')
        self.bin_id = A('bin_id')
        self.bin_ids_file_path = A('bin_ids_file')
        self.store_gene_detections_and_gene_coverages_tables = A('store_gene_detections_and_gene_coverages_tables')
        self.exclude_samples = A('exclude_samples')
        self.gene_coverages = pd.DataFrame.empty
        self.gene_detections = pd.DataFrame.empty
        self.samples = {}
        self.positive_samples = pd.DataFrame.empty
        self.number_of_positive_samples = None
        self.negative_samples = pd.DataFrame.empty
        self.number_of_negative_samples = None
        self.gene_class_information = pd.DataFrame.empty
        self.samples_information = pd.DataFrame.empty
        self.profile_db = {}
        self.gene_presence_absence_in_samples = pd.DataFrame.empty
        self.gene_coverages_filtered = pd.DataFrame.empty

        # check that there is a file like this
        if self.exclude_samples:
            filesnpaths.is_file_exists(self.exclude_samples)
            self.samples_to_exclude = set([l.split('\t')[0].strip() for l in open(args.exclude_samples, 'rU').readlines()])
            run.info('Excluding Samples', 'The following samples will be excluded: %s' % self.samples_to_exclude,)
        else:
            self.samples_to_exclude = set([])

        self.sanity_check()
        if self.profile_db_path is None:
            self.get_data_from_txt_file()
        else:
            # load sample list and gene_coverage_dict from the merged profile db
            args.init_gene_coverages = True
            if self.collection_name:
                self.summary = summarizer.ProfileSummarizer(args)
                self.summary.init()
            else:
                self.profile_db = ProfileSuperclass(args)
                self.profile_db.init_gene_coverages_and_detection_dicts()
                self.gene_coverages = pd.DataFrame.from_dict(self.profile_db.gene_coverages_dict, orient='index', dtype=float)
                self.gene_coverages.drop(self.samples_to_exclude, axis=1, inplace=True)
                self.Ng = len(self.gene_coverages.index)
                self.gene_detections = pd.DataFrame.from_dict(self.profile_db.gene_detection_dict, orient='index', dtype=float)
                self.gene_detections.drop(self.samples_to_exclude, axis=1, inplace=True)
                self.samples = set(self.gene_coverages.columns)


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
        self.check_if_valid_portion_value("alpha", self.alpha)

        # Checking beta
        if not isinstance(self.beta, float):
            raise ConfigError("beta value must be a type float.")
        self.check_if_valid_portion_value("beta", self.beta)
        if self.beta > self.alpha:
            raise ConfigError("beta value must be smaller than alpha value. The beta value you specified is %s while the alpha value\
            is %s" % (self.beta, self.alpha))

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


    def get_data_from_txt_file(self):
        """ Reads the coverage data from TAB delimited file """
        self.gene_coverages = pd.read_table(self.gene_coverages_data_file_path, sep='\t', header=0, index_col=0)
        self.gene_coverages.drop(self.samples_to_exclude, axis=1, inplace=True)
        self.Ng = len(self.gene_coverages.index)
        self.samples = set(self.gene_coverages.columns.values)
        # checking if a gene_detection file was also supplied
        if self.gene_detections_data_file_path:
            self.gene_detections = pd.read_table(self.gene_coverages_data_file_path, sep='\t', header=0, index_col=0)
            self.gene_detections.drop(self.samples_to_exclude, axis=1, inplace=True)
            # making sure that the tables are compatible, notice we're only checking if gene_detection contains everything that's in gene_coverages (and not vise versa)
            for gene_id in self.gene_coverages.index:
                if gene_id not in self.gene_detections.index:
                    raise ConfigError("Your tables are not compatible. For example gene_id %s is in %s, but not in %s" % (gene_id, self.gene_coverages_data_file_path,
                                                                                                                         self.gene_detections_data_file_path))
            for sample_id in self.samples:
                if sample_id not in self.gene_detections.columns.values:
                    raise ConfigError("Your tables are not compatible. For example sample_id %s is in %s, but not in %s" % (sample_id, self.gene_coverages_data_file_path,
                                                                                                                         self.gene_detections_data_file_path))


    def set_gene_presence_absence_in_samples(self):
        """ Determines the presence/absense of genes according to the gene detection threshold.
            The following arguments are populated:
                self.gene_presence_absence_in_samples
                self.gene_coverages_filtered
        """
        self.gene_presence_absence_in_samples = self.gene_detections > self.zeta
        self.gene_coverages_filtered = self.gene_coverages.copy()
        self.gene_coverages_filtered[~self.gene_presence_absence_in_samples] = 0


    def set_sample_detection_information(self):
        """ Using the --genome-presence-threshold and the --genome-absence-threhold the samples are devided to three groups:
                positive samples: samples in which the number of genes that are present (according to the --min-gene-detection threshold) is
                    greater than --genome-presence-threshold
                negative samples: samples in which the number of genes that are present (according to the --min-gene-detection threshold) is
                    smaller than --genome-absence-threhold
                ambiguous samples: all other samples (i.e. samples in which the number of genes that are present is between the two thresholds)
            This function populates the following arguments of self:
                self.positive_samples - a set of the positive sample ids self.negative_samples - a set of the negative sample ids
                self.samples_information - dictionary with True, False, None for positive, negative and ambiguous samples respectively
        """
        MDG_samples_information_table_name      = 'MDG_classifier_samples_information'
        MDG_samples_information_table_structure = ['samples', 'presence', 'number_of_detected_genes', 'number_of_taxon_specific_core_detected']
        MDG_samples_information_table_types     = ['str', 'bool', 'int', 'int']
        # create an empty dataframe
        samples_information = pd.DataFrame(index=self.samples, columns=MDG_samples_information_table_structure[1:])
        
        # Compute the number of detected genes per samples
        number_of_genes_in_sampels = self.gene_presence_absence_in_samples.sum(axis=0)
        for sample in self.samples:
            samples_information['number_of_detected_genes'][sample]= number_of_genes_in_sampels[sample]

        self.positive_samples = set(number_of_genes_in_sampels[number_of_genes_in_sampels/ self.Ng > self.alpha].index)
        self.negative_samples = set(number_of_genes_in_sampels[number_of_genes_in_sampels / self.Ng > self.beta].index) - self.positive_samples
        for sample in self.positive_samples:
            samples_information['presence'][sample] = True
        for sample in self.negative_samples:
            samples_information['presence'][sample] = False
        self.samples_information = samples_information

        self.number_of_positive_samples = len(self.positive_samples)
        self.number_of_negative_samples = len(self.negative_samples)
        self.run.warning('The number of positive samples is: %d ' % self.number_of_positive_samples)
        self.run.warning('The number of negative samples is: %d ' % self.number_of_negative_samples)
        self.run.warning('The number of ambiuous samples is: %d ' % (len(self.samples - self.negative_samples - self.positive_samples)))


    def get_taxon_specific_genes_in_samples(self, additional_description=''):
        """ Use only positive samples to identify the single copy core genes:
            Assumption: At least 25% of the genes are single copy and taxon specific (so that it is included in the Q1Q3 range).
                a. Sort coverage vaules
                b. Sliding window (with 50% overlapping windows) that has least outliers with the smallest range. The reason for 50% \
                overlap is that this guarentees that if the assumption that at least 25% are single copy taxon specific then, it is \
                guarentees that at least one window will have Q1Q3 composed only from taxon specific genes.
                    i.Run sliding window and for each step calculate the number of outliers and the range
                    ii.Sort according the number of outliers and then according to the 
        """
        # get the indixes that sort each column separately TODO: delete this part if you don't need it
        sorting_indexes = self.gene_coverages_filtered.apply(lambda x: x.sort_values(ascending=True).index)
        # TODO: I'm creating a copy but probably I don't need to
        coverages_pdf_output = self.output_file_prefix + additional_description + '-coverages.pdf'
        coverages_pdf_output_file = PdfPages(coverages_pdf_output)
        # creating dataframe for non-outliers (see usage below)
        non_outliers = pd.DataFrame().reindex_like(self.gene_coverages)
        sorting_indexes = self.gene_coverages_filtered.apply(lambda x: x.sort_values(ascending=True).index)
        print(self.samples_information)
        print(self.positive_samples)
        print(self.negative_samples)
        print(self.samples)
        for sample in self.positive_samples:
            # a vector of the coverages for the sample
            v = self.gene_coverages_filtered[sample].copy()
            # set coverage to zero if not above gene detection threshold
            v[self.gene_detections[sample] < self.zeta] = 0
            # get non-outliers according to the Interquartile range
            non_outliers[sample] = get_non_outliers(v)
            plot_outliers(coverages_pdf_output_file, v, non_outliers[sample], sample)
        coverages_pdf_output_file.close()

        # finding the genes that are non-outliers in the majority of samples
        # the required majority is defined by the user defined argument self.eta (--core-min-detection)
        non_outliers_all = non_outliers[non_outliers.sum(axis=1) >= self.number_of_positive_samples * self.eta].index
        return non_outliers_all


    def get_gene_specificity(self, detection_of_genes):
        """ Find all genes that are detected in """
        gene_specificity = {}
        for gene_id in detection_of_genes:
            if detection_of_genes[gene_id]['detected_in_non_positive_samples']:
            # if the gene is detected in at least one negative sample then gene_specificity is False
                gene_specificity[gene_id] = False
            elif detection_of_genes[gene_id]['number_of_detections'] == 0:
            # if the gene is not detected in any sample then the gene_specificity is None
                gene_specificity[gene_id] = None
            else:
            # if the gene is only detected in positive samples the gene_specificity is True
                gene_specificity[gene_id] = True
        return gene_specificity


    def get_coverage_consistency(self, adjusted_stds, detection_of_genes, beta):
        """For each gene if the adjusted standard deviation (to understand what this is refer to Alon Shaiber) is smaller
        than beta then coverage_consistency is True, otherwise, coverage_consistency is False. If the gene is not
        detected in positive samples then coverage consistency is not defined (get's a value 'None')"""
        coverage_consistency = {}

        for gene_id in adjusted_stds:
            # if the gene is not detected in any sample then return None
            if detection_of_genes[gene_id]['number_of_detections'] <= 1:
                coverage_consistency[gene_id] = None
            elif adjusted_stds[gene_id] == -1:
                coverage_consistency[gene_id] = None
            else:
                if adjusted_stds[gene_id] < beta:
                    coverage_consistency[gene_id] = True
                else:
                    coverage_consistency[gene_id] = False
        return coverage_consistency


    def get_taxon_specificity(self, coverage_consistency, gene_specificity):
        taxon_specificity = {}
        for gene_id in coverage_consistency:
            if coverage_consistency[gene_id] is None or gene_specificity[gene_id] is None:
                taxon_specificity[gene_id] = None
            elif coverage_consistency[gene_id] and gene_specificity[gene_id]:
                taxon_specificity[gene_id] = 'TS'
            else:
                taxon_specificity[gene_id] = 'TNS'
        return taxon_specificity


    def get_number_of_detections_for_gene(self, detection_of_genes, gene_id, samples):
        detections = 0
        for sample_id in samples:
            detections += detection_of_genes[gene_id][sample_id]
        return detections


    def get_core_accessory_info(self, detection_of_genes, gene_id, eta):
        """ Returns 'core'/'accessory' classification for each gene. This is done using only the samples in which the
        genome is detected """
        if detection_of_genes[gene_id]['number_of_detections'] == 0:
            return 'None'
        elif self.get_number_of_detections_for_gene(detection_of_genes, gene_id, self.positive_samples) < eta * len(
            self.positive_samples):
            return 'accessory'
        else:
            return 'core'


    def get_gene_class(self, taxon_specificity, core_or_accessory):
        if taxon_specificity == None or core_or_accessory == 'None':
            return 'None'
        elif taxon_specificity == 'TS':
            if core_or_accessory == 'core':
                return 'TSC'
            elif core_or_accessory == 'accessory':
                return 'TSA'
            else:
                print('%s is not valid. Value should be \'core\' or \'accessory\'' % core_or_accessory)
                exit(1)
        elif taxon_specificity == 'TNS':
            if core_or_accessory == 'core':
                return 'TNC'
            elif core_or_accessory == 'accessory':
                return 'TNA'
            else:
                print('%s is not valid. Value should be \'core\' or \'accessory\'' % core_or_accessory)
                exit(1)
        else:
            print('%s is not valid. Value should be \'TS\' or \'TNS\'' % taxon_specificity)
            exit(1)


    def report_gene_class_information(self):
        # TODO: change to pandas
        C = lambda dictionary, field, value : len([dict_id for dict_id in dictionary if dictionary[dict_id][field]==value])

        for gene_class in ['TSC', 'TSA', 'TNC', 'TNA', 'None']:
            self.run.info('Num class %s' % gene_class, C(self.gene_class_information, 'gene_class', gene_class))

        # TODO: report the number of negative samples and the number of NA samples
        self.run.info('Num samples in which the genome is detected', C(self.samples_information, 'detection', True), mc='green')


    def get_gene_classes(self):
        """ returning the classification per gene along with detection in samples (i.e. for each sample, whether the
        genome has been detected in the sample or not """
        # need to start a new gene_class_information dict
        # this is due to the fact that if the algorithm is ran on a list of bins then this necessary
        self.gene_class_information = pd.DataFrame(index=self.gene_coverages.index,columns=['gene_class'])
        # use a gene detection threshold to determine gene presence/absence in samples
        self.set_gene_presence_absence_in_samples()
        # use a gene presence threshold classify samples as positive, negative or ambiguous
        self.set_sample_detection_information()
        # 
        non_outliers_all = self.get_taxon_specific_genes_in_samples()

        # set the gene classes
        for gene_id in self.gene_coverages_filtered.index:
            if gene_id in non_outliers_all:
                self.gene_class_information['gene_class'][gene_id] = self.get_gene_class('TS', 'core')
            else:
                self.gene_class_information['gene_class'][gene_id] = self.get_gene_class(None,None)


    def get_specificity_from_class_id(self, class_id):
        try:
            class_id = int(class_id)
        except:
            raise ConfigError("Classes must be of type integer. You sent this: ", class_id)

        classes = {0: 'None',
                   1: 'TS',
                   2: 'TS',
                   3: 'TS',
                   4: 'TNS',
                   5: 'TNS'}

        try:
            return classes(class_id)
        except:
            raise ConfigError("The class id '%d' is not a valid one. Try one of these: '%s'" % (class_id, ', '.join(list(classes.keys()))))


    def save_gene_class_information_in_additional_layers(self, additional_description=''):
        if not self.additional_layers_to_append:
            additional_column_titles = []
            additional_layers_df = self.gene_class_information
        else:
            additional_layers_df = pd.read_table(self.additional_layers_to_append)
            try:
                # concatinating the gene_class_information with the user provided additional layers
                additional_layers_df.join(self.gene_class_information, how='outer')
            except ValueError as e:
                raise ConfigError("Something went wrong. This is what we know: %s. This could be happening because \
                you have columns in your --additional-layers file with the following title: %s" % (e, self.gene_class_information.columns.tolist()))

        if additional_description:
            additional_description = '-' + additional_description
        additional_layers_file_name = self.output_file_prefix + additional_description + '-additional-layers.txt'
        additional_layers_df.to_csv(additional_layers_file_name, sep='\t', index_label='gene_callers_id')


    def save_samples_information(self, additional_description=''):
        # TODO: there used to be this here:
        #self.run.warning(self.samples_information)
        if not self.samples_information_to_append:
            samples_information_df = self.samples_information
        else:
            samples_information_df = pd.read_table(self.samples_information_to_append)
            try:
                # concatinating the samples_information with the user provided samples_information file 
                samples_information_df.join(self.samples_information, how='outer')
            except ValueError as e:
                raise ConfigError("Something went wrong. This is what we know: %s. This could be happening because \
                you have columns in your --additional-layers file with the following title: %s" % (e, self.samples_information.columns.tolist()))

        if additional_description:
            additional_description = '-' + additional_description

        samples_information_file_name = self.output_file_prefix + additional_description + '-samples-information.txt'
        samples_information_df.to_csv(samples_information_file_name, sep='\t', index_label='samples')

    def save_gene_detection_and_coverage(self, additional_description=''):
        if additional_description:
            prefix = self.output_file_prefix + '-' + additional_description
        else:
            prefix = self.output_file_prefix
        gene_coverages_file_name = prefix + '-gene-coverages.txt'
        gene_detections_file_name = prefix + '-gene-detections.txt'
        self.gene_coverages.to_csv(gene_coverages_file_name, sep='\t', index_label='gene_callers_id')
        self.gene_detections.to_csv(gene_detections_file_name, sep='\t', index_label='gene_callers_id')


    def get_coverage_and_detection_dict(self,bin_id):
        _bin = summarizer.Bin(self.summary, bin_id)
        self.gene_coverages = pd.DataFrame.from_dict(_bin.gene_coverages, orient='index', dtype=float)
        print(self.gene_coverages)
        self.gene_coverages.drop(self.samples_to_exclude, axis=1, inplace=True)
        self.Ng = len(self.gene_coverages.index)
        self.gene_detections = pd.DataFrame.from_dict(_bin.gene_detection, orient='index', dtype=float)
        self.gene_detections.drop(self.samples_to_exclude, axis=1, inplace=True)
        self.samples = set(self.gene_coverages.columns.values)


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
                self.get_gene_classes()
                self.save_gene_class_information_in_additional_layers(bin_id)
                self.save_samples_information(bin_id)
                if self.store_gene_detections_and_gene_coverages_tables:
                    self.save_gene_detection_and_coverage(bin_id)

        else:
            # No collection provided so running on the entire detection table
            self.get_gene_classes()
            self.save_gene_class_information_in_additional_layers()
            self.save_samples_information()
            if self.store_gene_detections_and_gene_coverages_tables:
                self.save_gene_detection_and_coverage()


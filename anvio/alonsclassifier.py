# -*- coding: utf-8
# pylint: disable=line-too-long
"""
    Classes to classify genes based on coverages across metagenomes.

    anvi-alons-classifier is the default client using this module
"""


import anvio
import numpy as np
import anvio.utils as utils
import anvio.terminal as terminal
import anvio.summarizer as summarizer
import anvio.filesnpaths as filesnpaths

from anvio.errors import ConfigError
from anvio.dbops import ProfileSuperclass


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


class AlonsClassifier:
    def __init__(self, args, run=run, progress=progress):
        self.run = run
        self.progress = progress

        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        self.gene_coverages_data_file_path = A('data_file')
        self.gene_detection_data_file_path = A('gene_detection_data_file')
        self.profile_db_path = A('profile_db')
        self.output_file_prefix = A('output_file_prefix')
        self.alpha = A('alpha')
        self.beta = A('beta')
        self.gamma = A('gamma')
        self.eta = A('eta')
        self.zeta = A('zeta')
        self.additional_layers_to_append = A('additional_layers_to_append')
        self.samples_information_to_append = A('samples_information_to_append')
        self.number_of_positive_samples = None
        self.collection_name = A('collection_name')
        self.bin_id = A('bin_id')
        self.bin_ids_file_path = A('bin_ids_file')
        self.store_gene_detections_and_gene_coverages_tables = A('store_gene_detections_and_gene_coverages_tables')
        self.gene_coverages = {}
        self.gene_detection = {}
        self.samples = {}
        self.positive_samples = {}
        self.negative_samples = {}
        self.gene_class_information = {}
        self.samples_information = {}
        self.profile_db = {}

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
                self.gene_coverages = self.profile_db.gene_coverages_dict
                self.gene_detection = self.profile_db.gene_detection_dict
                self.samples = set(next(iter(self.gene_coverages.values())).keys())


    def sanity_check(self):
        """Basic sanity check for class inputs"""

        if self.profile_db_path is None and self.gene_coverages_data_file_path is None:
            raise ConfigError("You must provide either a profile.db or a gene coverage matrix data file")

        if self.profile_db_path and self.gene_coverages_data_file_path:
            raise ConfigError("You provided both a profile database and a gene coverage matrix data file, you \
            must provide only one or the other (hint: if you have a profile database, the use that")

        # checking output file
        filesnpaths.is_output_file_writable(self.output_file_prefix + '-additional-layers.txt', ok_if_exists=False)
        # checking alpha
        if not isinstance(self.alpha, float):
            raise ConfigError("alpha value must be a type float.")
        if self.alpha <= 0 or self.alpha > 1:
            raise ConfigError("alpha value must be greater than zero and a max of 1, the value you supplied %s" % self.alpha)

        # Checking beta
        if not isinstance(self.beta, float):
            raise ConfigError("beta value must be a type float.")
        if self.beta <= 0:
            raise ConfigError("beta value must be greater than zero, the value you supplied %s" % self.beta)

        # Checking gamma
        if not isinstance(self.gamma, float):
            raise ConfigError("Gamma value must be a type float.")
        if self.gamma <= 0:
            raise ConfigError("Gamma value must be greater than zero, the value you supplied %s" % self.gamma)

        # Checking eta
        if self.eta <= 0 or self.eta > 1:
            raise ConfigError("eta value must be greater than zero and a max of 1, the value you supplied %s" % self.eta)

        if self.collection_name:
            if not self.profile_db_path:
                raise ConfigError("You specified a collection name %s, but you provided a gene coverage matrix data file \
                 collections are only available when working with a profile database." % self.collection_name)


    def get_data_from_txt_file(self):
        """ Reads the coverage data from TAB delimited file """
        self.samples = set(utils.get_columns_of_TAB_delim_file(self.gene_coverages_data_file_path))
        self.gene_coverages = utils.get_TAB_delimited_file_as_dictionary(self.gene_coverages_data_file_path, column_mapping=[int] + [float] * len(self.samples))
        # checking if a gene_detection file was also supplied
        if self.gene_detection_data_file_path:
            self.gene_detection = utils.get_TAB_delimited_file_as_dictionary(self.gene_detection_data_file_path, column_mapping=[int] + [float] * len(self.samples))
            # making sure that the tables are compatible, notice we're only checking if gene_detection contains everything that's in gene_coverages (and not vise versa)
            for gene_id in self.gene_coverages:
                if gene_id not in self.gene_detection:
                    raise ConfigError("Your tables are not compatible. For example gene_id %s is in %s, but not in %s" % (gene_id, self.gene_coverages_data_file_path,
                                                                                                                         self.gene_detection_data_file_path))
            gene_detection_sample_list = next(iter(self.gene_detection.values())).keys()
            for sample_id in next(iter(self.gene_coverages.values())).keys():
                if sample_id not in gene_detection_sample_list:
                    raise ConfigError("Your tables are not compatible. For example sample_id %s is in %s, but not in %s" % (sample_id, self.gene_coverages_data_file_path,
                                                                                                                         self.gene_detection_data_file_path))


    def apply_func_to_genes_in_sample(self, func, list_of_genes=None):
        """ Apply the give function on the list of genes in each sample. The function is expected to accept a list """
        if not list_of_genes:
            list_of_genes = self.gene_coverages.keys()
        d = dict(zip(self.samples, [next(map(func, [[self.gene_coverages[gene_id][sample_id] for gene_id in list_of_genes]])) for
                           sample_id in self.samples]))
        return d


    def get_mean_coverage_in_samples(self, list_of_genes=None):
        """ Returns a dictionary with of the average coverage value of the list of genes per sample. if no list of genes is
        supplied then the average is calculated over all genes """
        if not self.samples:
            # if all samples don't contain the genome then return 0 for mean value
            return 0
        else:
            mean_coverage_in_samples = self.apply_func_to_genes_in_sample(np.mean, list_of_genes)
            return mean_coverage_in_samples


    def get_std_in_samples(self, list_of_genes=None):
        """ Returns a dictionary with of the standard deviation of the coverage values of the list of genes per sample.
        if no list of genes is supplied then the average is calculated over all genes """
        std_in_samples = self.apply_func_to_genes_in_sample(np.std, list_of_genes)
        return std_in_samples


    def get_detection_of_genes(self, mean_coverage_in_samples, std_in_samples):
        """ Returns a dictionary (of dictionaries), where for each gene_id, and each sample_id the detection of the gene
        is determined. The criteria for detection is having coverage that is greater than 0 and also that is not more
        than gamma (default is gamma=3) standard deviations below the mean coverage in the sample.
        Notice that the mean coverage isn't the mean of all genes in the sample necesarilly. In fact it would be the mean of
        only the taxon-specific genes."""
        detection_of_genes = {}
        non_zero_non_detections = False
        for gene_id in self.gene_coverages:
            detection_of_genes[gene_id] = {}
            detection_of_genes[gene_id]['number_of_detections'] = 0
            detection_of_genes[gene_id]['detected_in_non_positive_samples'] = False
            for sample in self.samples:
                if sample in self.positive_samples or sample not in self.negative_samples:
                # getting gene detection according to coverage criteria
                # samples that are ambiguous (neither negative nor positive)
                    detection_of_genes[gene_id][sample] = self.gene_coverages[gene_id][sample] > max(0,mean_coverage_in_samples[sample] -
                                                                                 self.gamma*std_in_samples[sample])
                else:
                # if the sample is a negative sample or the presence of the genome is unknown then the threshold for detections is
                # more stringent and requires coverage significantly higher than the coverage of the TSC genes
                    detection_of_genes[gene_id][sample] = self.gene_coverages[gene_id][sample] > max(0,mean_coverage_in_samples[sample] +
                                                                                 self.gamma*std_in_samples[sample])
                    detection_of_genes[gene_id]['detected_in_non_positive_samples'] += detection_of_genes[gene_id][sample]
                if self.gene_detection:
                    # if we have the gene detection (previously known as "percent covered") information then we will also use it to determine detection in samples:
                    gene_detection_above_threshold = self.gene_detection[gene_id][sample] > self.zeta
                    detection_of_genes[gene_id][sample] = detection_of_genes[gene_id][sample] * gene_detection_above_threshold
                detection_of_genes[gene_id]['number_of_detections'] += detection_of_genes[gene_id][sample]
                if self.gene_coverages[gene_id][sample] > 0 and self.gene_coverages[gene_id][sample] < mean_coverage_in_samples[sample] - \
                        self.gamma*std_in_samples[sample]:
                    non_zero_non_detections = True
        if non_zero_non_detections:
            # print('gene %s, in some sample has non-zero coverage %s, and it has been marked as not detected due '
                  # 'to the detection criteria' % (gene_id, data[gene_id][sample]))
                  print('some genes in some samples were marked as not detected due to the detection criteria')

        return detection_of_genes


    def get_samples_information(self, detection_of_genes, alpha, genes_to_consider=None):
        '''
        Setting the values for the samples_information dictionary with the logical variable detection and the number of TSC genes 
        for each sample. In addition the total number of positive samples (i.e. samples that contain the reference genome) saved in 
        number_of_positive_samples.
        '''
        if not genes_to_consider:
            # if no list of genes is supplied then considering all genes
            genes_to_consider = detection_of_genes.keys()
        samples_information = {}
        self.number_of_positive_samples = 0
        for sample_id in self.samples:
            samples_information[sample_id] = {}
            number_of_detected_genes_in_sample = len([gene_id for gene_id in genes_to_consider if detection_of_genes[gene_id][sample_id]])
            samples_information[sample_id]['detection'] = number_of_detected_genes_in_sample > alpha * len(
                genes_to_consider)
            # There should be a better way to do this, but for now, if there is an intermediate value of genes that is
            # detected then the sample would be considered neither positive nor negative
            if number_of_detected_genes_in_sample > alpha * len(genes_to_consider):
                samples_information[sample_id]['detection'] = True
            elif number_of_detected_genes_in_sample > 0.5 * alpha * len(genes_to_consider):
                samples_information[sample_id]['detection'] = None
            else:
                samples_information[sample_id]['detection'] = False
            self.number_of_positive_samples += 1 if samples_information[sample_id]['detection'] else 0
            samples_information[sample_id]['number_of_detected_genes'] = number_of_detected_genes_in_sample
        self.samples_information = samples_information
        self.run.warning('The number of positive samples is: %s ' % self.number_of_positive_samples)


    def get_adjusted_std_for_gene_id(self, gene_id, mean_coverage_in_samples, detection_of_genes):
        """Returns the adjusted standard deviation for a gene_id. The adjusted standard deviation is calculated 
        only according to the positive samples (ones that contain the genome), and only in positive samples that 
        contain the gene (according to the gene detection criteria)."""
        # Note: originally I thought I would only consider samples in which the genome was detected, but in fact,
        # if a gene is detected in a sample in which the genome is not detected then that is a good sign that this is
        #  a TNS gene. But I still kept here the original definition of adjusted_std
        # adjusted_std = np.std([d[gene_id, sample_id] / mean_coverage_in_samples[sample_id] for sample_id in samples if (
        #         detection_of_genes[gene_id][sample_id] and samples_information[sample_id])])
        
        # First find all positive samples that also contain the gene
        positive_samples_with_gene = []
        for sample_id in self.positive_samples:
            if detection_of_genes[gene_id][sample_id]:
                positive_samples_with_gene.append(sample_id)
        if not positive_samples_with_gene:
            # If no positive samples contain the gene then we want the coverage consistency to be None
            return -1, 0
        else:
            adjusted_gene_coverages = [self.gene_coverages[gene_id][sample_id]/mean_coverage_in_samples[sample_id] for sample_id in
                                   positive_samples_with_gene]
            adjusted_mean = np.mean(adjusted_gene_coverages)
            adjusted_std = np.std(adjusted_gene_coverages)
            return adjusted_std, adjusted_mean


    def get_adjusted_stds(self, mean_coverage_in_samples, detection_of_genes):
        """ Returns a dictionary with the gene ids as keys and the adjusted standard deviation per gene"""
        adjusted_std = {}
        adjusted_mean = {}
        # Make sure that there is at least one positive sample
        if not self.number_of_positive_samples:
            raise ConfigError("It seems that the reference genome is not detected in any of your samples. \
            You might consider changing some of the parameters (for example --min-gene-detection, --min-portion-of-TSC \
            or --core-min-detection).")
            
        for gene_id in self.gene_coverages:
            adjusted_std[gene_id], adjusted_mean[gene_id] = self.get_adjusted_std_for_gene_id(gene_id, mean_coverage_in_samples, detection_of_genes)
        return adjusted_std, adjusted_mean


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

    def get_loss_function_value(self, taxon_specificity, adjusted_stds, beta):
        loss = 0
        for gene_id in taxon_specificity:
            if taxon_specificity[gene_id] == 'TS':
                # Notice: here adjusted std includes the samples that don't have the genome detected in them (it kind of
                # makes sense, because if the gene is detected even though the genome is not, then maybe it's not
                # taxon-specific
                loss += adjusted_stds[gene_id]
            else:
                loss += beta
        return loss


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
        C = lambda dictionary, field, value : len([dict_id for dict_id in dictionary if dictionary[dict_id][field]==value])

        for gene_class in ['TSC', 'TSA', 'TNC', 'TNA', 'None']:
            self.run.info('Num class %s' % gene_class, C(self.gene_class_information, 'gene_class', gene_class))

        # TODO: report the number of negative samples and the number of NA samples
        self.run.info('Num samples in which the genome is detected', C(self.samples_information, 'detection', True), mc='green')


    def get_gene_classes(self):
        """ returning the classification per gene along with detection in samples (i.e. for each sample, whether the
        genome has been detected in the sample or not """
        TSC_genes = set(self.gene_coverages.keys())
        converged = False
        loss = None
        self.gene_class_information = {}
        # Initializing all the samples to be positive
        self.positive_samples = self.samples
        self.negative_samples = {}
        self.adjusted_stds = 0
        self.adjusted_mean = dict.fromkeys(TSC_genes,1)

        while not converged:
            # mean of coverage of all TS genes in each sample
            mean_coverage_of_TS_in_samples = self.get_mean_coverage_in_samples(TSC_genes)

            # Get the standard deviation of the taxon-specific genes in a sample
            # TODO: right now, single copy, and multi-copy genes would be treated identically. Hence, multi-copy genes
            # would skew both the mean and the std of the taxon-specific genes.
            std_of_TS_in_samples = self.get_std_in_samples(TSC_genes)
            detection_of_genes = self.get_detection_of_genes(mean_coverage_of_TS_in_samples, std_of_TS_in_samples)
            self.get_samples_information(detection_of_genes, self.alpha, TSC_genes)
            self.positive_samples = {sample_id for sample_id in self.samples if self.samples_information[sample_id]['detection']}
            self.negative_samples = {sample_id for sample_id in self.samples if self.samples_information[sample_id]['detection'] == False}
            self.adjusted_stds, self.adjusted_mean = self.get_adjusted_stds(mean_coverage_of_TS_in_samples,detection_of_genes)
            coverage_consistency = self.get_coverage_consistency(self.adjusted_stds, detection_of_genes, self.beta)
            gene_specificity = self.get_gene_specificity(detection_of_genes)
            taxon_specificity = self.get_taxon_specificity(coverage_consistency, gene_specificity)
            new_loss = self.get_loss_function_value(taxon_specificity, self.adjusted_stds, self.beta)
            epsilon = 2 * self.beta

            if loss is not None:
                converged = True
            loss = new_loss

            self.run.warning('current value of loss function: %s ' % loss)

            for gene_id in self.gene_coverages:
                # setup a dict for gene id:
                g = {}

                g['gene_specificity'] = gene_specificity[gene_id]
                g['gene_coverage_consistency'] = coverage_consistency[gene_id]
                g['number_of_detections'] = detection_of_genes[gene_id]['number_of_detections']
                g['core_or_accessory'] = self.get_core_accessory_info(detection_of_genes, gene_id, self.eta)
                g['gene_class'] = self.get_gene_class(taxon_specificity[gene_id], g['core_or_accessory'])
                g['adjusted_stds'] = self.adjusted_stds[gene_id]
                g['adjusted_mean'] = self.adjusted_mean[gene_id]

                # counting the number of positive samples that contain the gene
                g['detection_in_positive_samples'] = len([sample_id for sample_id in self.positive_samples if detection_of_genes[gene_id][sample_id]])

                # Getting the portion of positive samples that contain the gene
                g['portion_detected'] = g['detection_in_positive_samples'] / len(self.positive_samples) if g['detection_in_positive_samples'] else 0


                self.gene_class_information[gene_id] = g

            TSC_genes = {gene_id for gene_id in self.gene_class_information if self.gene_class_information[gene_id]['gene_class']=='TSC'}

            self.report_gene_class_information()

        self.get_samples_information(detection_of_genes, self.alpha, genes_to_consider=TSC_genes)


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
            additional_layers_dict = self.gene_class_information
        else:
            additional_column_titles = utils.get_columns_of_TAB_delim_file(self.additional_layers_to_append)
            additional_layers_dict = utils.get_TAB_delimited_file_as_dictionary(self.additional_layers_to_append,
                                                                                dict_to_append=self.gene_class_information,
                                                                                assign_none_for_missing=True,
                                                                                column_mapping=[int] + [str] * len(additional_column_titles))

        if additional_description:
            additional_description = '-' + additional_description

        additional_layers_file_name = self.output_file_prefix + additional_description + '-additional-layers.txt'
        headers = headers=['gene_callers_id', 'gene_class', 'number_of_detections', 'portion_detected','gene_specificity','gene_coverage_consistency','core_or_accessory', 'adjusted_mean', 'adjusted_stds'] + additional_column_titles

        utils.store_dict_as_TAB_delimited_file(additional_layers_dict, additional_layers_file_name, headers=headers)


    def save_samples_information(self, additional_description=''):
        if not self.samples_information_to_append:
            samples_information_column_titles = list(self.samples_information[next(iter(self.samples_information))])
            samples_information_dict = self.samples_information
        else:
            samples_information_column_titles = utils.get_columns_of_TAB_delim_file(self.samples_information_to_append)
            column_mapping = [str] * (len(samples_information_column_titles) + 2)
            self.run.warning(self.samples_information)
            samples_information_dict = utils.get_TAB_delimited_file_as_dictionary(self.samples_information_to_append,
                                                                                  dict_to_append=self.samples_information,
                                                                                  assign_none_for_missing=True,
                                                                                  column_mapping=column_mapping)

        if additional_description:
            additional_description = '-' + additional_description

        samples_information_file_name = self.output_file_prefix + additional_description + '-samples-information.txt'
        utils.store_dict_as_TAB_delimited_file(samples_information_dict, samples_information_file_name,
                                                   headers=['samples'] + samples_information_column_titles)


    def save_gene_detection_and_coverage(self, additional_description=''):
        if additional_description:
            prefix = self.output_file_prefix + '-' + additional_description
        else:
            prefix = self.output_file_prefix
        gene_coverages_file_name = prefix + '-gene-coverages.txt'
        gene_detections_file_name = prefix + '-gene-detections.txt'
        utils.store_dict_as_TAB_delimited_file(self.gene_coverages, gene_coverages_file_name)
        utils.store_dict_as_TAB_delimited_file(self.gene_detection, gene_detections_file_name)


    def get_coverage_and_detection_dict(self,bin_id):
        _bin = summarizer.Bin(self.summary, bin_id)
        self.gene_coverages = _bin.gene_coverages
        self.gene_detection = _bin.gene_detection
        self.samples = set(next(iter(self.gene_coverages.values())).keys())


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


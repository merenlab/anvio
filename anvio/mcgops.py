# -*- coding: utf-8
# pylint: disable=line-too-long
"""
    Classes to classify genes based on coverages across metagenomes.

    anvi-mcg-classifier is the default client using this module
"""

import os
import anvio
import numpy as np
import matplotlib.pyplot as plt
import anvio.terminal as terminal


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


class MCGPlots:
    ''' Class to handle the plots generated for the MCG-classifier'''
    def __init__(self, mcg, gene_id, args=object, run=run, progress=progress):
        self.run = run
        self.progress = progress

        A = lambda x,y: args.__dict__[x] if x in args.__dict__ else y
        self.plot_err = A('plot_error', True)
        self.plot_fit = A('plot_fit', False)
        self.plot_confidence_intervals = A('plot_confidence_intervals', True)
        self.plot_gene_nuc_coverage = A('plot_gene_nuc_coverage', False)
        self.plot_sample_nuc_coverage = A('plot_sample_nuc_coverage', False)
        self.save = A('save', True)
        self.fit = {}
        self.gene_nuc_coverage_dict = {}
        self.samples_nuc_coverage_dict = {}
        self.gene_id = gene_id

        self.samples = set(mcg.gene_presence_absence_in_samples.loc[gene_id, mcg.gene_presence_absence_in_samples.loc[gene_id,]==True].index)
        # taking only the positive samples in which the gene is present
        self.samples = self.samples.intersection(mcg.positive_samples)
        if not self.samples:
            # if the gene is not detected in any positive samples then dont save the plot
            self.save = False

        self.x = mcg.samples_coverage_stats_dicts.loc[self.samples,'non_outlier_mean_coverage']
        self.y = mcg.gene_level_coverage_stats_dict_of_dataframes['non_outlier_mean_coverage'].loc[gene_id, self.samples]
        self.std_x = mcg.samples_coverage_stats_dicts.loc[self.samples,'non_outlier_coverage_std']
        self.std_y = mcg.gene_level_coverage_stats_dict_of_dataframes['non_outlier_coverage_std'].loc[gene_id, self.samples]

        if self.plot_gene_nuc_coverage:
            for sample in self.samples:
            # getting outlier and non outlier nucleotide coverage for the gene
                non_outliers = mcg.gene_level_coverage_stats_dict_of_dataframes['gene_coverage_values_per_nt'][gene_id][sample][mcg.gene_level_coverage_stats_dict_of_dataframes['non_outlier_positions'][gene_id][sample]]
                outlier_positions = np.logical_not(mcg.gene_level_coverage_stats_dict_of_dataframes['non_outlier_positions'][gene_id][sample])
                outliers = mcg.gene_level_coverage_stats_dict_of_dataframes['gene_coverage_values_per_nt'][gene_id][sample][outlier_positions]

                self.gene_nuc_coverage_dict[sample] = {}
                self.gene_nuc_coverage_dict[sample]['non_outliers'] = non_outliers
                self.gene_nuc_coverage_dict[sample]['x_non_outliers'] = np.array([mcg.samples_coverage_stats_dicts.loc[sample,'non_outlier_mean_coverage']] * len(non_outliers))
                self.gene_nuc_coverage_dict[sample]['outliers'] = outliers
                self.gene_nuc_coverage_dict[sample]['x_outliers'] = np.array([mcg.samples_coverage_stats_dicts.loc[sample,'non_outlier_mean_coverage']] * len(outliers))

        if self.plot_sample_nuc_coverage:
           for sample in self.samples:
                sample_non_outliers = mcg.coverage_values_per_nt[sample][mcg.non_outlier_indices[sample]]
                self.samples_nuc_coverage_dict[sample] = {}
                self.samples_nuc_coverage_dict[sample]['non_outliers'] = sample_non_outliers
                self.samples_nuc_coverage_dict[sample]['y'] = np.array([self.y[sample]] * len(sample_non_outliers))

        if gene_id in mcg.gene_coverage_consistency_dict:
            self.fit = mcg.gene_coverage_consistency_dict[gene_id]
            self.plot_fit = self.fit['converged']
        self.output_file_prefix = mcg.output_file_prefix


    def plot(self):
        plt.scatter(self.x,self.y,c='black')
        axes = plt.gca()
        if self.plot_err and len(self.samples) > 1:
            # plot error box - a box with width and hight corresponding to the
            # standard deciation in sample nucleotide coverage and gene nucleotide
            # coverage correspondingly.
            #x_err_width = (axes.get_xlim()[1] - axes.get_xlim()[0]) / 20
            #y_err_width = (axes.get_ylim()[1] - axes.get_ylim()[0]) / 20

            for i in range(len(self.x)):
                plot_err(self.x[i], self.std_x[i], self.y[i], self.std_y[i])
                plot_err(self.y[i], self.std_y[i], self.x[i], self.std_x[i], axis=1)

        if self.plot_fit:
            # plot the curve of the linear fit
            fit_slope = self.fit['slope']
            fit_slope_std = self.fit['slope_std']
            fit_R_squered = self.fit['R_squered']
            fit_slope_precision = self.fit['slope_precision']

            # add slope
            f = lambda b: lambda x: b*x
            plt.plot(self.x, np.apply_along_axis(f(fit_slope),0,self.x),'c')

            if self.plot_confidence_intervals:
                # add confidence intervals to the plot
                plt.plot(self.x, np.apply_along_axis(f(fit_slope + 2*fit_slope_std),0,self.x),'g')
                plt.plot(self.x, np.apply_along_axis(f(fit_slope - 2*fit_slope_std),0,self.x),'g')

            # adding the text to the plot
            text_for_fit = u'$R^2 = %.2f$\n $slope = %.2f$\n $slope\ precision = %.2f$' % (fit_R_squered, fit_slope, fit_slope_precision)
            axes.text(0.2, 0.9, text_for_fit, ha='center', va='center', transform=axes.transAxes)

        if self.plot_gene_nuc_coverage:
            for sample in self.samples:
                plt.scatter(self.gene_nuc_coverage_dict[sample]['x_non_outliers'], self.gene_nuc_coverage_dict[sample]['non_outliers'], c='b',s=1)
                plt.scatter(self.gene_nuc_coverage_dict[sample]['x_outliers'], self.gene_nuc_coverage_dict[sample]['outliers'], c='r',s=1)

        if self.plot_sample_nuc_coverage:
            for sample in self.samples:
                plt.scatter(self.samples_nuc_coverage_dict[sample]['non_outliers'], self.samples_nuc_coverage_dict[sample]['y'],c='b',s=1)

        if self.save:
            output_dir = self.output_file_prefix + '_gene_consistency_plots'
            os.makedirs(output_dir, exist_ok=True)
            output_file = output_dir + '/' + str(self.gene_id) + '.png'
            plt.savefig(output_file, format='png')
        plt.close()


def plot_err(x1,err_x1,x2,err_width,color='m',axis=0):
    a = [x1 - err_x1, x1 + err_x1]
    b1 = [x2 - err_width, x2 - err_width]
    b2 = [x2 + err_width, x2 + err_width]
    if axis == 0:
        plt.plot(a,b1,color)
        plt.plot(a,b2,color)
    elif axis == 1:
        plt.plot(b1,a,color)
        plt.plot(b2,a,color)

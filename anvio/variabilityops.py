# -*- coding: utf-8
# pylint: disable=line-too-long

"""Classes to make sense of sequence variation"""


import os
import sys
import copy
import random
import inspect
import argparse
import numpy as np
import pandas as pd
import operator as op

from scipy.stats import entropy

import anvio
import anvio.tables as t
import anvio.dbops as dbops
import anvio.utils as utils
import anvio.terminal as terminal
import anvio.constants as constants
import anvio.filesnpaths as filesnpaths
import anvio.ccollections as ccollections
import anvio.structureops as structureops
import anvio.auxiliarydataops as auxiliarydataops

from anvio.errors import ConfigError


__author__ = "Developers of anvi'o (see AUTHORS.txt)"
__copyright__ = "Copyleft 2015-2018, the Meren Lab (http://merenlab.org/)"
__credits__ = ['Alon Shaiber']
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "A. Murat Eren"
__email__ = "a.murat.eren@gmail.com"


pd.options.display.max_columns=100
pd.options.display.max_rows=100
pp = terminal.pretty_print
progress = terminal.Progress()
run = terminal.Run(width=62)


class VariabilityFilter:
    def __init__(self, args={}, p=progress, r=run):
        np.seterr(invalid='ignore')
        self.stealth_filtering = False

        self.known_kwargs = [
            "min_filter", "min_condition",
            "max_filter", "max_condition",
            "subset_filter", "subset_condition",
        ]


    def __setattr__(self, key, value):
        self.__dict__[key] = value
        if key == "df":
            self.__dict__[self.name] = value


    def filter_by_occurrence(self, min_occurrence=None):
        if not min_occurrence:
            min_occurrence = self.min_occurrence

        if min_occurrence <= 1:
            return

        occurrences = self.df["unique_pos_identifier_str"].value_counts()
        postions_with_greater_than_min_occurrence = occurrences[occurrences >= min_occurrence].index

        func = self.subset_filter
        filter_args = {"criterion": "unique_pos_identifier_str",
                       "subset_values": postions_with_greater_than_min_occurrence}

        self.filter_wrapper(func, descriptor="Minimum occurrence of positions", value_for_run_info=min_occurrence, **filter_args)


    def filter_by_minimum_coverage_in_each_sample(self, min_coverage_in_each_sample=None):
        """To remove any unique entry from the variable positions table that describes a variable position
           and yet is not helpful to distinguish samples from each other."""
        if not min_coverage_in_each_sample:
            min_coverage_in_each_sample = self.min_coverage_in_each_sample

        if min_coverage_in_each_sample < 1:
            return

        min_cov_each_pos = self.df.groupby("unique_pos_identifier")["coverage"].min()
        pos_to_keep = list(min_cov_each_pos[min_cov_each_pos >= min_coverage_in_each_sample].index)

        func = self.subset_filter
        filter_args = {"criterion": "unique_pos_identifier",
                       "subset_values": pos_to_keep}

        self.filter_wrapper(func, descriptor="Minimum coverage across all samples", value_for_run_info=min_coverage_in_each_sample, **filter_args)


    def filter_by_num_positions_from_each_split(self, num_positions_from_each_split=None):
        if not num_positions_from_each_split:
            num_positions_from_each_split = self.num_positions_from_each_split

        if num_positions_from_each_split <= 0:
            return

        subsample_func = lambda x: pd.Series(x.unique()) if len(x.unique()) <= num_positions_from_each_split else\
                                   pd.Series(np.random.choice(x.unique(), size=num_positions_from_each_split, replace=False))
        unique_positions_to_keep = self.df.groupby('split_name')['unique_pos_identifier'].apply(subsample_func)

        func = self.subset_filter
        filter_args = {"criterion": "unique_pos_identifier",
                       "subset_values": unique_positions_to_keep}

        self.filter_wrapper(func, descriptor="Num positions from each split", value_for_run_info=num_positions_from_each_split, **filter_args)


    def filter_by_scattering_factor(self):
        """ this is what filter by scattering factor is supposed to do (it has never done this)
        Using the --min-scatter parameter you can eliminate some SNV positions based on how they
        partition samples. This one is a bit tricky, but Meren wants to keep it in the code base. If
        you skip this you will not lose anything, but for the nerd kind, this is how it goes: If you
        have N samples in your dataset, a given variable position x in one of your splits can split
        your N samples into t groups based on the identity of the variation they harbor at position
        x. For instance, t would have been 1, if all samples had the same type of variation at
        position x (which would not be very interesting, because in this case position x would have
        zero contribution to a deeper understanding of how these samples differ based on
        variability. When t > 1, it would mean that identities at position x across samples do
        differ. But how much scattering occurs based on position x when t > 1? If t=2, how many
        samples would end up in each group? Obviously, the even distribution of samples across
        groups may tell us something different than uneven distribution of samples across groups.
        So, this parameter filters out any x if ‘the number of samples in the second largest group’
        (=scatter) is less than the value you choose. Here is an example: lets assume you have 7
        samples. While 5 of those have AG, 2 of them have TC at position x. This would mean the
        scatter of x is 2. If you set -m to 2, this position would not be reported in your output
        matrix. The default value for -m is 0, which means every x found in the database and
        survived previous filtering criteria will be reported. Naturally, -m can not be more than
        half of the number of samples.)"""
        raise ConfigError("Woops. The function that handles --min-scatter doesn't do \
                           what we thought it did. This will be fixed maybe. Sorry for \
                           the inconvenience.")


    def are_passed_arguments_valid(self, kwargs):
        if not self.criterion and not self.passed_function:
            raise ConfigError("VariabilityFilter :: pass a filter criterion (criterion=) if you want to apply a standard filter,\
                               or if you want to apply a special filter, pass a filter function (function=).\
                               You passed neither.")

        if self.criterion and self.passed_function:
            raise ConfigError("VariabilityFilter :: pass a filter criterion (criterion=) if you want to apply a standard filter,\
                               or if you want to apply a special filter, pass a filter function (function=).\
                               You passed both.")

        self.is_df_exists()
        self.is_df_a_dataframe()
        self.is_filter_criterion_valid()
        self.is_passed_function_valid()
        self.is_passed_kwargs_valid(kwargs)


    def filter_data(self, name='data', criterion=None, function=None, verbose=False, **kwargs):
        """ FIXME add example here """
        self.name = name
        self.filter_info_log = {}
        self.criterion = criterion
        self.passed_function = function
        self.append_info_log("criterion", self.criterion)
        self.append_info_log("passed_function", self.passed_function)

        try:
            self.are_passed_arguments_valid(kwargs)
            filter_routine = self.get_filter_function()
            params = self.get_filter_params(filter_routine, kwargs)
            filter_routine(**params)

        except ConfigError as e:
            header = "Filtering failed. Here is the log:"
            self.report_filter_info_log(header)
            print(e)
            sys.exit()

        else:
            if verbose:
                header = "Filter passed without raising error. Here is the log:"
                self.report_filter_info_log(header, success=True)


    def remove_condition_params(self, params):
        params_to_remove = [x for x in params.keys() if x.endswith("_condition")]
        return {k: v for k, v in params.items() if k not in params_to_remove}


    def get_filter_params(self, function, kwargs):
        if self.passed_function:
            self.append_info_log("final filter params", kwargs)
            return kwargs

        else:
            filter_and_condition_params = kwargs
            if filter_and_condition_params:
                filter_and_condition_params = self.remove_filter_params_if_condition_params_are_not_true(filter_and_condition_params)
                filter_params = self.remove_condition_params(filter_and_condition_params)
            else:
                filter_params = self.infer_filter_params_from_criterion(function)

            self.append_info_log("final filter params", filter_params)
            return filter_params


    def remove_filter_params_if_condition_params_are_not_true(self, params):
        filter_condition_keys = [x for x in params.keys() if x.endswith("_condition")]

        params_to_remove = []
        for filter_condition in filter_condition_keys:
            if not params[filter_condition]:
                filt = filter_condition.replace("_condition", "_filter")
                params_to_remove.append(filt)
                self.append_info_log("`%s` was not true. no filter for" % (filter_condition), filt)

        return {k: v for k, v in params.items() if k not in params_to_remove}


    def infer_filter_params_from_criterion(self, function):
        """This function should never be called because the programmer should always provide explicit
           filtering parameters to filter_data. But if none are passed because they are not known,
           this function tries to infer the filtering parameters"""
        self.append_info_log("attempting filter param inference", True)

        candidate_args = {
            "min_filter": ["min_" + self.criterion], # e.g. ["min_departure_from_consensus"]
            "max_filter": ["max_" + self.criterion],
            "subset_filter": [self.criterion + "s_of_interest", self.criterion + "_of_interest"],
        }

        candidate_attribute_names = [x for y in candidate_args for x in candidate_args[y]]
        self.append_info_log("looking for attribute names", ", ".join(candidate_attribute_names))

        filter_params = {}
        arg_val_names_found = []
        for arg_key in candidate_args:
            for arg_val_name in candidate_args[arg_key]:
                if hasattr(self, arg_val_name):
                    arg_val_names_found.append(arg_val_name)
                    filter_params[arg_key] = getattr(self, arg_val_name)
                    break

        self.append_info_log("found attribute names", ", ".join(arg_val_names_found))
        if not filter_params:
            raise ConfigError("VariabilityFilter :: Explicit parameters were not passed to\
                               filter_data, so anvi'o tried its best to find attributes that made\
                               sense given your filter criterion `%s`, but failed.")
        return filter_params


    def get_filter_function(self):

        if self.passed_function:
            function = self.passed_function
        else:
            function = self.general_serial_filtering

        self.append_info_log("filter function", function)
        return function


    def subset_filter(self, criterion, subset_values):
        subset_values = self.get_subset_values_as_set(subset_values)
        return self.df.loc[self.df[criterion].isin(subset_values)]


    def get_subset_values_as_set(self, subset_values):
        try:
            return set(subset_values)
        except:
            try:
                return set([subset_values])
            except:
                raise TypeError


    def extremum_filter(self, criterion, extremum_value, extremum_type, inclusive = True):
        """extremum_value: float or int
               The value you are using as a cut off.
           extremum_type: str
               Either 'max' or 'min'.
           inclusive: bool (default True)
               If True, values equal to extremum_value will not be removed"""
        extremum_value = float(extremum_value)
        op = self.get_bound_filter_operator(extremum_type, inclusive)
        return self.df[op(self.df[criterion], extremum_value)]


    def get_bound_filter_operator(self, extremum_type, inclusive):
        if extremum_type == 'min' and inclusive:
            return op.ge
        if extremum_type == 'min' and not inclusive:
            return op.gt
        if extremum_type == 'max' and inclusive:
            return op.le
        if extremum_type == 'max' and not inclusive:
            return op.lt


    def filter_wrapper(self, filter_func, descriptor=None, value_for_run_info=None, **filter_args):
        if self.stealth_filtering or not (descriptor and str(value_for_run_info)):
            self.df = filter_func(**filter_args)

        else:
            key, val = self.generate_message_for_run_info(descriptor, value_for_run_info)
            self.run.info(key, val)
            self.progress.new('Filtering based on %s' % (descriptor))
            entries_before = len(self.df.index)
            self.df = filter_func(**filter_args)
            entries_after = len(self.df.index)
            self.progress.end()
            self.report_change_in_entry_number(entries_before, entries_after, reason=descriptor.lower())

        self.check_if_data_is_empty()


    def gen_filter_wrapper_args_for_serial_filtering(self, keyword_name, keyword_value):
        D = lambda x, y: "{} {}".format(x, y).replace("_", " ").capitalize() if not self.stealth_filtering else lambda x, y: None
        d = {}

        if keyword_name == "min_filter":
            d["filter_function"]    = self.extremum_filter
            d["descriptor"]         = D("minimum", self.criterion)
            d["value_for_run_info"] = keyword_value
            d["filter_args"]        = {"extremum_value" : keyword_value,
                                       "criterion"      : self.criterion,
                                       "extremum_type"  : "min",
                                       "inclusive"      : True}

        elif keyword_name == "max_filter":
            d["filter_function"]    = self.extremum_filter
            d["descriptor"]         = D("maximum", self.criterion)
            d["value_for_run_info"] = keyword_value
            d["filter_args"]        = {"extremum_value" : keyword_value,
                                       "criterion"      : self.criterion,
                                       "extremum_type"  : "max",
                                       "inclusive"      : True}

        elif keyword_name == "subset_filter":
            d["filter_function"]    = self.subset_filter
            d["descriptor"]         = D(self.criterion+"(s)", "of interest")
            d["value_for_run_info"] = keyword_value
            d["filter_args"]        = {"subset_values" : keyword_value,
                                       "criterion"     : self.criterion}

        else:
            raise ConfigError("VariabilitySuper :: `%s` is not a keyword recognized by the built in\
                               function methods of this class." % (keyword_name))

        return d


    def general_serial_filtering(self, **params):
        for param in params:
            d = self.gen_filter_wrapper_args_for_serial_filtering(param, params[param])
            self.filter_wrapper(d["filter_function"],
                                d["descriptor"],
                                d["value_for_run_info"],
                                **d["filter_args"])

    def is_passed_kwargs_compatible_with_passed_function(self, kwargs):
        params_inspection = inspect.signature(self.passed_function).parameters

        for param in params_inspection:
            has_default = True if params_inspection[param].default is not inspect._empty else False
            if not has_default and param not in kwargs.keys():
                raise ConfigError("`%s` was passed to filter_data. All its arguments without defaults must\
                                   also be passed, but `%s` was not. Do so with \"%s = ...\"" \
                                       % (self.passed_function, param, param))
            else:
                self.append_info_log("`%s` of `%s` was passed or has default" % (param, self.passed_function), True)

        bad_args = [arg for arg in kwargs.keys() if arg not in params_inspection]
        if bad_args:
            raise ConfigError("VariabilityFilter :: The args [%s] were passed to filter_data,\
                               but are not args of `%s`. Available args for this function are [%s]" \
                               % (", ".join(bad_args), self.passed_function, ", ".join(params_inspection)))
        for kwarg in kwargs:
            self.append_info_log("filter argument: `%s`" % (kwarg), "valid argument for `%s`" % (self.passed_function))


    def is_passed_kwargs_compatible_with_built_in_function(self, kwargs):
        for kwarg in kwargs:
            if kwarg in self.known_kwargs:
                if kwarg.endswith("_condition") and kwarg.replace("_condition","_filter") not in kwargs:
                    self.append_info_log("condition argument: `%s`" % (kwarg), "missing filter partner: `%s`" % (kwarg.replace("_condition","_filter")))
                    raise ConfigError("VariabilityFilter :: The argument `%s` takes something that\
                                       evaluates to True or False that determines whether or not the\
                                       data will be filtered by the argument `%s`, which must\
                                       also be passed." % (kwarg, kwarg.replace("_condition","_filter")))
                else:
                    self.append_info_log("filter argument: `%s`" % (kwarg), "known kwarg")

            else:
                self.append_info_log("filter argument: `%s`" % (kwarg), "unknown kwarg")
                raise ConfigError("VariabilityFilter :: Argument `%s` was passed to filter_data,\
                                   but is unknown. Unknown arguments are only allowed if\
                                   a function is passed to filter_data. Otherwise generalized function\
                                   methods will be employed based on known arguments you supply.\
                                   Here are all known arguments: [%s]" % (kwarg, ", ".join(self.known_kwargs)))


    def is_passed_kwargs_valid(self, kwargs):
        if self.passed_function:
            self.is_passed_kwargs_compatible_with_passed_function(kwargs)
        else:
            self.is_passed_kwargs_compatible_with_built_in_function(kwargs)


    def is_passed_function_valid(self):
        if not self.passed_function:
            return

        if not callable(self.passed_function):
            self.append_info_log("function is callable", False)
            raise ConfigError("Function %s, which was passed to filter_data, is not callable" \
                                   % (self.passed_function))
        else:
            self.append_info_log("function is callable", True)


    def is_df_exists(self):
        try:
            self.df = getattr(self, self.name)
        except AttributeError as e:
            self.append_info_log("%s is an object" % self.name, False)
            raise ConfigError("VariabilityFilter :: You tried to filter the object `%s` which is\
                               not an attribute of VariabilitySuper or its parent classes." \
                                   % (self.name))
        self.append_info_log("%s is an object" % self.name, True)


    def is_df_a_dataframe(self):
        if type(self.df) != pd.core.frame.DataFrame:
            self.append_info_log("name points to a valid dataframe", False)
            raise ConfigError("VariabilityFilter :: You tried to filter the object `%s` which is\
                               of type `%s`. You can only filter pandas dataframes." \
                                   % (self.name, type(self.df)))
        self.append_info_log("name points to a valid dataframe", True)


    def is_filter_criterion_valid(self):
        if not self.criterion:
            self.append_info_log("criterion was passed", False)
            return

        if not self.is_filter_criterion_a_column_in_dataframe():
            raise ConfigError("VariabilityFilter :: The filter criterion `%s` does not exist as\
                               a column in self.df. It must." % (self.criterion))


    def is_filter_criterion_a_column_in_dataframe(self):
        is_column = True if self.criterion in self.df.columns else False
        self.append_info_log("%s is a column in self.df" % (self.criterion), is_column)
        return is_column


    def append_info_log(self, key, val):
        self.filter_info_log[key] = val


    def report_filter_info_log(self, header="Filter progress log:", success=False):
        run.info_single(header, nl_before=1, nl_after=0, mc="yellow" if success else "red")
        for key, val in self.filter_info_log.items():
            run.info(key, val)
        run.info("end of log", True)


    def generate_message_for_run_info(self, descriptor, value):
        if isinstance(value, str):
            return descriptor, value

        if isinstance(value, float) or isinstance(value, int):
            return descriptor, str(value)

        try:
            if not len(value):
                return descriptor, "None"
            elif len(value) < 100:
                return descriptor, ", ".join([str(x) for x in value])
            else:
                return "Number of {}".format(descriptor), len(value)
        except TypeError:
            return descriptor, value


class VariabilitySuper(VariabilityFilter, object):
    def __init__(self, args={}, p=progress, r=run):
        np.seterr(invalid='ignore')
        self.args = args

        if args.engine not in variability_engines:
            raise ConfigError("You are doing something wrong :/ Focus '%s' does not correspond to an available engine." % args.engine)

        A = lambda x, t: t(args.__dict__[x]) if x in args.__dict__ else None
        null = lambda x: x
        # variability
        self.data = args.data if 'data' in args.__dict__ else pd.DataFrame({})
        # splits
        self.bin_id = A('bin_id', null)
        self.collection_name = A('collection_name', null)
        self.splits_of_interest = A('splits_of_interest_set', set)
        self.splits_of_interest_path = A('splits_of_interest', null)
        # database
        self.profile_db_path = A('profile_db', null)
        self.contigs_db_path = A('contigs_db', null)
        self.structure_db_path = A('structure_db', null)
        # genes
        self.gene_caller_ids = A('gene_caller_ids', null)
        self.genes_of_interest = A('genes_of_interest_set', set)
        self.genes_of_interest_path = A('genes_of_interest', null)
        # samples
        self.sample_ids_of_interest = A('samples_of_interest_set', set)
        self.sample_ids_of_interest_path = A('samples_of_interest', null)
        # filtering
        self.min_scatter = A('min_scatter', int) or 0
        self.min_occurrence = A('min_occurrence', int) or 1
        self.min_coverage = A('min_coverage', int) or 0
        self.min_coverage_in_each_sample = A('min_coverage_in_each_sample', int) or 0
        self.min_departure_from_reference = A('min_departure_from_reference', float) or 0
        self.max_departure_from_reference = A('max_departure_from_reference', float) or 1
        self.min_departure_from_consensus = A('min_departure_from_consensus', float) or 0
        self.max_departure_from_consensus = A('max_departure_from_consensus', float) or 1
        self.num_positions_from_each_split = A('num_positions_from_each_split', int) or 0
        # output
        self.quince_mode = A('quince_mode', bool)
        self.compute_gene_coverage_stats = A('compute_gene_coverage_stats', bool)
        self.output_file_path = A('output_file', null)
        self.only_if_structure = A('only_if_structure', bool)
        self.skip_sanity_check = A('skip_sanity_check', bool) or False
        self.include_split_names_in_output = A('include_split_names', null)
        self.include_contig_names_in_output = A('include_contig_names', null)
        self.skip_comprehensive_variability_scores = A('skip_comprehensive_variability_scores', bool) or False

        self.append_structure_residue_info = True if self.structure_db_path else False
        self.table_provided = False if self.data.empty else True
        self.load_all_genes = True
        self.load_all_samples = True
        self.min_departure_from_reference_filtered = False
        self.max_departure_from_reference_filtered = False
        self.substitution_scoring_matrices = None
        self.merged_split_coverage_values = None
        self.unique_pos_identifier = 0
        self.split_name_position_dict = {}
        self.unique_pos_id_to_entry_id = {}
        self.contig_sequences = None
        self.input_file_path = None

        self.comprehensive_stats_headers = []
        self.comprehensive_variability_scores_computed = False

        VariabilityFilter.__init__(self, self.args, r=self.run, p=self.progress)

        # f = function called in self.process
        # **kwargs = parameters passed to function
        F = lambda f, **kwargs: (f, kwargs)
        self.process_functions = [F(self.init_commons),
                                  F(self.load_variability_data),
                                  F(self.load_structure_data),
                                  F(self.apply_preliminary_filters),
                                  F(self.set_unique_pos_identification_numbers),
                                  F(self.filter_data, function=self.filter_by_num_positions_from_each_split),
                                  F(self.compute_additional_fields),
                                  F(self.filter_data, criterion="departure_from_consensus",
                                                      min_filter=self.min_departure_from_consensus,
                                                      min_condition=self.min_departure_from_consensus > 0,
                                                      max_filter=self.max_departure_from_consensus,
                                                      max_condition=self.max_departure_from_consensus < 1),
                                  F(self.recover_base_frequencies_for_all_samples),
                                  F(self.filter_data, function=self.filter_by_minimum_coverage_in_each_sample),
                                  F(self.compute_comprehensive_variability_scores),
                                  F(self.compute_gene_coverage_fields),
                                  F(self.get_residue_structure_information,)]

        if not self.skip_sanity_check:
            self.sanity_check()

        # Initialize the contigs super
        if self.contigs_db_path:
            filesnpaths.is_file_exists(self.contigs_db_path)
            dbops.ContigsSuperclass.__init__(self, self.args, r=self.run, p=self.progress)
            self.init_contig_sequences()

        # these lists are dynamically extended
        self.columns_to_report = {
            'position_identifiers': [
                ('entry_id', int),
                ('unique_pos_identifier', int),
                ('pos', int),
                ('pos_in_contig', int),
                ('contig_name', str),
                ('split_name', str),
            ],
            'sample_info': [
                ('sample_id', str),
            ],
            'gene_info': [
                ('corresponding_gene_call', int),
                ('in_partial_gene_call', int),
                ('in_complete_gene_call', int),
                ('base_pos_in_codon', int),
                ('codon_order_in_gene', int),
                ('codon_number', int),
                ('gene_length', int),
                ('gene_coverage', float),
                ('non_outlier_gene_coverage', float),
                ('non_outlier_gene_coverage_std', float),
            ],
            'coverage_info': [
                ('coverage', int),
                ('mean_normalized_coverage', float),
                ('cov_outlier_in_split', int),
                ('cov_outlier_in_contig', int),
            ],
            'sequence_identifiers': [
                ('reference', str),
                ('consensus', str),
                ('competing_nts', str),
                ('competing_aas', str),
                ('competing_codons', str),
            ],
            'statistical': [
                ('departure_from_reference', float),
                ('departure_from_consensus', float),
                ('n2n1ratio', float),
                ('entropy', float),
                ('kullback_leibler_divergence_raw', float),
                ('kullback_leibler_divergence_normalized', float),
                ('synonymity', float),
            ],
            'SSMs': [
            ],
            'structural': [
            ],
        }
        self.columns_to_report_order = ['position_identifiers', 'sample_info', 'gene_info', 'coverage_info',
                                        'sequence_identifiers', 'statistical', 'SSMs', 'structural']



    def sanity_check(self):
        if self.engine not in variability_engines:
            raise ConfigError("The VariabilitySuper class is inherited with an unknown engine.\
                               WTF is '%s'? Anvi'o needs an adult :(" % self.engine)

        if not self.table_provided:
            if not self.contigs_db_path:
                raise ConfigError("You need to provide a contigs database.")

            if not self.profile_db_path:
                raise ConfigError("You need to provide a profile database.")

        if self.table_provided and (self.contigs_db_path or self.profile_db_path):
            raise ConfigError("You provided a variability table (--variability-profile), but you\
                               also provided a contigs database and/or a profile database. You need\
                               to supply either a variability table, or, both a profile database and\
                               a contigs database combo.")

        if self.table_provided and (self.collection_name or self.bin_id):
            raise ConfigError("You provided a variability table (--variability-profile), but you\
                               also provided a collection name and/or a bin id, which are parameters\
                               only used in conjunction with a profile database and a contigs\
                               database. No big deal, just remove these parameters from your\
                               command.")

        if self.table_provided and (self.splits_of_interest_path or self.splits_of_interest):
            if "split_name" not in self.data.columns:
                raise ConfigError("Your variability profile does not have a split_name column, and\
                                   therefore you cannot provide splits of interest\
                                   (--splits-of-interest).")

        if self.genes_of_interest and (self.genes_of_interest_path or self.gene_caller_ids):
            raise ConfigError("VariabilitySuper: you initialized me with self.genes_of_interest\
                               because you're a programmer and you know what you're doing. But you\
                               also initialized me with self.genes_of_interest_path and/or\
                               self.gene_caller_ids. Because you didn't skip sanity_check(), I am\
                               complaining. If you skip sanity_check, self.genes_of_interest_path\
                               and self.gene_caller_ids will be ignored.")

        if self.sample_ids_of_interest and self.sample_ids_of_interest_path:
            raise ConfigError("VariabilitySuper: you initialized me with self.sampes_of_interest\
                               because you're a programmer and you know what you're doing.  But you\
                               also initialized me with self.sampes_of_interest_path. Because you\
                               didn't skip sanity_check(), I am complaining. If you skip\
                               sanity_check, self.sampes_of_interest_path.")

        if self.splits_of_interest and (self.splits_of_interest_path or self.bin_id or self.collection_name):
            raise ConfigError("VariabilitySuper: you initialized me with self.splits_of_interest\
                               because you're a programmer and you know what you're doing. But you\
                               also initialized me with one/all/some of\
                               self.splits_of_interest_path, self.bin_id, self.collection_name.\
                               Because you didn't skip sanity_check(), I am complaining. If you skip\
                               sanity_check, self.splits_of_interest_path, self.bin_id,\
                               self.collection_name will all be ignored.")


    def convert_counts_to_frequencies(self, retain_counts = False):
        if retain_counts:
            freq_columns = [x + '_freq' for x in self.items]
            self.data[freq_columns] = self.data[self.items].divide(self.data['coverage'], axis = 0)
        else:
            self.data[self.items] = self.data[self.items].divide(self.data['coverage'], axis = 0)


    def get_sample_ids_of_interest(self, sample_ids_of_interest_path=""):
        """It is essential to note that sample_ids_of_interest_path is "" by default, not None. The
           programmer can pass None to avoid the argument defaulting to a class-wide attribute
           (first code block of this method), which may not exist for classes inheriting this
           method.
        """
        # use class-wide attribute if no parameters is passed
        if sample_ids_of_interest_path is "":
            sample_ids_of_interest_path = self.sample_ids_of_interest_path

        if sample_ids_of_interest_path:
            filesnpaths.is_file_tab_delimited(sample_ids_of_interest_path, expected_number_of_fields=1)
            sample_ids_of_interest = set([s.strip() for s in open(sample_ids_of_interest_path).readlines()])

        else:
            sample_ids_of_interest = set([])

        return sample_ids_of_interest


    def get_genes_of_interest(self, genes_of_interest_path="", gene_caller_ids=""):
        """It is essential to note that genes_of_interest_path and gene_caller_ids are "" by
           default, not None. The programmer can pass None to avoid the argument defaulting to a
           class-wide attribute (first code block of this method), which may not exist for classes
           inheriting this method.
        """
        # use class-wide attributes if no parameters are passed
        if genes_of_interest_path is "" and gene_caller_ids is "":
            genes_of_interest_path = self.genes_of_interest_path
            gene_caller_ids = self.gene_caller_ids

        if genes_of_interest_path and gene_caller_ids:
            self.progress.end()
            raise ConfigError("You can't provide gene caller ids from the command line, and a list\
                               of gene caller ids as a file at the same time, obviously.")

        if gene_caller_ids:
            if "," in gene_caller_ids:
                gene_caller_ids = [g.strip() for g in gene_caller_ids.split(",")]
            else:
                gene_caller_ids = [gene_caller_ids]
            for index, gene_caller_id in enumerate(gene_caller_ids):
                try:
                    gene_caller_ids[index] = int(gene_caller_id)
                except:
                    raise ConfigError("Anvi'o does not like your gene caller id '%s'..." % gene_caller_id)
            genes_of_interest = set(gene_caller_ids)

        elif genes_of_interest_path:
            filesnpaths.is_file_tab_delimited(genes_of_interest_path, expected_number_of_fields=1)

            try:
                genes_of_interest = set([int(s.strip()) for s in open(genes_of_interest_path).readlines()])
            except ValueError:
                self.progress.end()
                raise ConfigError("Well. Anvi'o was working on your genes of interest ... and ... \
                                   those gene IDs did not look like anvi'o gene caller ids :/ Anvi'o\
                                   is now sad.")

        else:
            # looks like no genes were specified
            genes_of_interest = set([])

        return genes_of_interest


    def get_splits_of_interest(self, splits_of_interest_path="", split_source=""):
        """It is essential to note that splits_of_interest_path and split_source are "" by default,
           not None. The programmer can pass None to avoid the argument defaulting to a class-wide
           attribute (first code block of this method), which may not exist for classes inheriting
           this method.
        """
        # use class-wide attributes if no parameters are passed
        if split_source is "" and splits_of_interest_path is "":
            split_source = self.split_source
            splits_of_interest_path = self.splits_of_interest_path

        if not split_source:
            splits_of_interest = set([])

        elif split_source == "gene_caller_ids":
            splits_of_interest = list(set([self.gene_callers_id_to_split_name_dict[g] for g in self.genes_of_interest]))

        elif split_source == "bin_id":
            if self.collection_name and not self.bin_id:
                raise ConfigError('When you declare a collection name, you must also declare a bin id\
                                   (from which the split names of interest will be acquired).')
            if self.bin_id and not self.collection_name:
                raise ConfigError("You declared a bin id but anvi'o doesn't know which collection\
                                   it comes from. Please provide a collection name.")
            splits_of_interest = ccollections.GetSplitNamesInBins(self.args).get_split_names_only()

        elif split_source == "split_names":
            filesnpaths.is_file_tab_delimited(splits_of_interest_path, expected_number_of_fields=1)
            splits_of_interest = set([c.strip().replace('\r', '') for c in open(splits_of_interest_path).readlines()])

        else:
            raise ConfigError("Invalid split source '%s'. Expected 'split_names', 'bin_id', or\
                               'gene_caller_ids'." % split_source)

        return splits_of_interest


    def get_items(self):
        # Ensure self.items is sorted alphabetically.  This is required for resolving ties in
        # coverage alphabetically, which is described in the docstring of
        # self.compute_additional_fields.
        if self.engine == 'NT':
            self.items = sorted(constants.nucleotides)
        elif self.engine == 'CDN':
            self.items = sorted(constants.codons)
        elif self.engine == 'AA':
            self.items = sorted(constants.amino_acids)
        else:
            raise ConfigError("https://goo.gl/sx6JHg :(")
        self.columns_to_report['coverage_info'].extend([(x, str) for x in self.items])


    def get_substitution_scoring_matrices(self):
        import anvio.data.SSMs as SSMs
        self.substitution_scoring_matrices = SSMs.get(self.engine, self.run)
        for m in self.substitution_scoring_matrices:
            self.columns_to_report['SSMs'].extend([(m, int), (m + '_weighted', float)])


    def init_commons(self):
        self.progress.new('Init')

        if self.only_if_structure and not self.structure_db_path:
            self.progress.end()
            raise ConfigError("You can't ask to only include genes with structures \
                               (--only-if-structure) without providing a structure database.")

        self.progress.update('Checking the output file path ...')
        if self.output_file_path:
            filesnpaths.is_output_file_writable(self.output_file_path)

        if not self.sample_ids_of_interest:
            self.progress.update('Setting up samples of interest ...')
            self.sample_ids_of_interest = self.get_sample_ids_of_interest()

        if not self.genes_of_interest:
            self.progress.update('Setting up genes of interest ...');
            # self.genes_of_interest can be injected into this class programatically; this
            # conditional method call prevents overwriting
            self.genes_of_interest = self.get_genes_of_interest()

        # ways to get splits of interest: 1) genes of interest, 2) bin id, 3) directly
        self.progress.update('Setting up splits of interest ...')
        self.check_how_splits_are_found()

        if not self.splits_of_interest:
            # self.splits_of_interest can be injected into this class programatically; this
            # conditional method call prevents overwriting
            self.splits_of_interest = self.get_splits_of_interest()

        if self.genes_of_interest:
            genes_available = self.gene_callers_id_to_split_name_dict if not self.table_provided else self.data["corresponding_gene_call"].unique()
            # check for genes that do not appear in the contigs database
            bad_gene_caller_ids = [g for g in self.genes_of_interest if g not in genes_available]
            if bad_gene_caller_ids:
                self.progress.end()
                some_to_report = bad_gene_caller_ids[:5] if len(bad_gene_caller_ids) <= 5 else bad_gene_caller_ids
                raise ConfigError("{} of the gene caller ids you provided {} not {}. {}: {}. You only have 2 lives left.\
                                   2 more mistakes, and anvi'o will automatically uninstall itself. Yes, seriously :(".\
                                   format(len(bad_gene_caller_ids),
                                          "is" if len(bad_gene_caller_ids) == 1 else "are",
                                          "in this variability table" if self.table_provided else "known to this contigs database",
                                          "Here are a few of those ids" if len(some_to_report) > 1 else "Its id is",
                                          ", ".join([str(x) for x in some_to_report])))

        self.progress.update('Making sure you are not playing games ..')
        if self.engine not in ['NT', 'CDN', 'AA']:
            raise ConfigError("Anvi'o doesn't know what to do with a engine on '%s' yet :/" % self.engine)
        self.table_structure = t.variable_nts_table_structure if self.engine ==  'NT' else t.variable_codons_table_structure

        # set items of interest
        self.get_items()
        # populate substitution scoring matrices
        self.progress.end()
        self.get_substitution_scoring_matrices()
        self.progress.new('Init')

        # if we already have variability data we are almost done
        if self.table_provided:
            self.check_if_data_is_empty()
            self.available_sample_ids = self.data["sample_id"].unique()
            self.is_available_samples_compatible_with_sample_ids_of_interest()
            self.progress.end()
            return

        # otherwise we have more work to do
        self.progress.update('Making sure our databases are here ..')
        if not self.profile_db_path:
            raise ConfigError('You need to provide a profile database.')

        if not self.contigs_db_path:
            raise ConfigError('You need to provide a contigs database.')

        if self.append_structure_residue_info and self.engine not in ["AA", "CDN"]:
            raise ConfigError('You provided a structure database, which is only compatible with --engine AA and --engine CDN')

        self.progress.update('Making sure our databases are compatible ..')
        utils.is_profile_db_and_contigs_db_compatible(self.profile_db_path, self.contigs_db_path)
        if self.append_structure_residue_info:
            utils.is_structure_db_and_contigs_db_compatible(self.structure_db_path, self.contigs_db_path)

        if self.min_coverage_in_each_sample and not self.quince_mode:
            self.progress.end()
            raise ConfigError("When you specify a coverage value through --min-coverage-in-each-sample, you must also\
                                use --quince-mode flag, since the former parameter needs to know the coverage values in all\
                                samples even if variation is reported for only one sample among others. This is the only way\
                                to figure out whether variation is not reported for other samples due to low or zero coverage,\
                                or there was no variation to report despite the high coverage. Anvi'o could turn --quince-mode\
                                flat automatically for you, but it is much better if you have full control and understanding\
                                of what is going on.")

        if self.quince_mode:
            self.progress.update('Accessing auxiliary data file ...')
            auxiliary_data_file_path = dbops.get_auxiliary_data_path_for_profile_db(self.profile_db_path)
            if not os.path.exists(auxiliary_data_file_path):
                raise ConfigError("Anvi'o needs the auxiliary data file to run this program with '--quince-mode' flag.\
                                    However it wasn't found at '%s' :/" % auxiliary_data_file_path)
            self.merged_split_coverage_values = auxiliarydataops.AuxiliaryDataForSplitCoverages(auxiliary_data_file_path, None, ignore_hash=True)

        self.input_file_path = '/' + '/'.join(os.path.abspath(self.profile_db_path).split('/')[:-1])

        self.progress.update('Reading the profile database ...')
        profile_db = dbops.ProfileDatabase(self.profile_db_path)
        self.available_sample_ids = sorted(list(profile_db.samples))
        self.is_available_samples_compatible_with_sample_ids_of_interest()

        if not profile_db.meta['SNVs_profiled']:
            self.progress.end()
            raise ConfigError("Well well well. It seems SNVs were not characterized for this profile database.\
                                Sorry, there is nothing to report here!")

        profile_db.disconnect()
        self.progress.end()


    def gen_sqlite_where_clause_for_variability_table(self):
        """It is impractical to load the entire variability table and then filter it according to our
           splits_of_interest, sample_ids_of_interest, and genes_of_interest. For example, what if
           genes_of_interest = set([0]) in a profile database with 50,000 genes? Why is splits of interest
           not included here? Because split_name is not a column in the variable codon table."""
        R = lambda x, y: self.run.info("%s that variability data will be fetched for" % \
                         (x.capitalize() if len(y)<200 else "Num of "+x), ", ".join([str(z) for z in y]) if len(y)<200 else len(y))

        conditions = {}
        if self.sample_ids_of_interest:
            conditions['sample_id'] = ' IN (%s)' % ','.join(['"{}"'.format(x) for x in self.sample_ids_of_interest])
            R("sample ids", self.sample_ids_of_interest)
            self.load_all_samples = False

        if self.genes_of_interest:
            conditions['corresponding_gene_call'] = ' IN (%s)' % ','.join(['{}'.format(x) for x in self.genes_of_interest])
            R("gene ids", self.genes_of_interest)
            self.load_all_genes = False

        if self.min_departure_from_reference > 0:
            self.run.info("Only fetch entries with min departure from reference", self.min_departure_from_reference)
            conditions['departure_from_reference'] = ' >= {}'.format(self.min_departure_from_reference)
            self.min_departure_from_reference_filtered = True

        if self.max_departure_from_reference < 1:
            self.run.info("variability data fetched only for entries with max departure from reference", self.min_departure_from_reference)
            conditions['departure_from_reference'] = ' <= {}'.format(self.max_departure_from_reference)
            self.max_departure_from_reference_filtered = True

        if not conditions:
            return ""

        return " AND ".join([col_name + col_condition for col_name, col_condition in conditions.items()])


    def load_variability_data(self):
        """Populates self.data (type pandas.DataFrame) from profile database tables."""
        if self.table_provided:
            return

        sqlite_where_clause = self.gen_sqlite_where_clause_for_variability_table()

        self.progress.new('Generating variability')
        self.progress.update('Reading the profile database ...')
        profile_db = dbops.ProfileDatabase(self.profile_db_path)

        if self.engine == 'NT':
            self.data = profile_db.db.get_table_as_dataframe(t.variable_nts_table_name,
                                                             table_structure=self.table_structure,
                                                             where_clause=sqlite_where_clause)

        elif self.engine == 'CDN' or self.engine == 'AA':
            if not profile_db.meta['SCVs_profiled']:
                raise ConfigError("It seems codon frequencies were not computed for this profile database,\
                                   therefore there is nothing to report here for codon or amino acid variability\
                                   profiles :(")

            self.data = profile_db.db.get_table_as_dataframe(t.variable_codons_table_name,
                                                             where_clause=sqlite_where_clause)
            self.check_if_data_is_empty()

            # this is where magic happens for the AA engine. we just read the data from the variable codons table, and it
            # needs to be turned into AAs if the engine is AA.
            if self.engine == 'AA':
                self.convert_item_coverages()
                self.convert_reference_info()

            # FIXME should this be in the variability table to begin with?
            self.data["split_name"] = self.data["corresponding_gene_call"].apply(lambda x: self.gene_callers_id_to_split_name_dict[x])

        self.data["codon_number"] = utils.convert_sequence_indexing(self.data["codon_order_in_gene"], source="M0", destination="M1")

        # we're done here. bye.
        profile_db.disconnect()
        self.progress.end()


    def load_structure_data(self):
        if not self.append_structure_residue_info:
            return

        # open up residue_info table from structure db
        self.progress.new('Loading structure information')
        self.progress.update('Reading the structure database ...')
        structure_db = structureops.StructureDatabase(self.structure_db_path)
        self.structure_residue_info = structure_db.db.get_table_as_dataframe(t.structure_residue_info_table_name)

        self.genes_with_structure = set(self.structure_residue_info["corresponding_gene_call"].unique())
        # genes_included = genes_of_interest, unless genes_of_interest weren't specified. then
        # it equals all genes in self.data. I don't overwrite self.genes_of_interest because a
        # needless gene filtering step will be carried out if self.genes_of_interest is not an
        # empty set
        genes_included = self.genes_of_interest if self.genes_of_interest else set(self.data["corresponding_gene_call"].unique())
        genes_with_var_and_struct = set([g for g in self.genes_with_structure if g in genes_included])

        if not genes_with_var_and_struct:
            self.progress.end()
            raise ConfigError("There is no overlap between genes that have structures and genes that have variability.\
                               Consider changing things upstream in your workflow or do not provide the structure db.\
                               Here are the genes in your structure database: {}".\
                               format(", ".join([str(x) for x in self.genes_with_structure])))

        if self.only_if_structure:
            # subset self.genes_of_interest to those that have structure
            self.genes_of_interest = genes_with_var_and_struct

        # we're done here. bye.
        structure_db.disconnect()
        self.progress.end()


    def check_if_data_is_empty(self):
        if self.data.empty:
            raise self.EndProcess


    def check_how_splits_are_found(self):
        """splits of interest are specified either by providing the splits of interest directly, or
           by providing a collection and bin. Alternatively, splits can be inferred from genes of
           interest. These three routes for determining splits of interest are mutually exclusive and
           we make sure the user/programmer provides parameters for one route only.
        """
        if self.table_provided:
            # splits of interest can still be specified if table was provided
            self.split_source = "split_names" if self.splits_of_interest_path else None
            return

        requested_split_source = {
            "gene_caller_ids": True if self.genes_of_interest_path or self.gene_caller_ids or self.genes_of_interest else False,
            "split_names":     True if self.splits_of_interest_path or self.splits_of_interest else False,
            "bin_id":          True if self.bin_id or self.collection_name else False
           }

        if not any(list(requested_split_source.values())):
            raise ConfigError("You must specify a list of genes (with --gene-caller-ids or\
                               --genes-of-interest), OR a list of splits (--splits-of-interest), OR\
                               a collection and bin combo (--collection-name and bin-id). You\
                               supplied none of these parameters and so anvi'o doesn't know what you\
                               want. If you are truly interested in everything, you\
                               should run the script anvi-script-add-default-collection, and then\
                               supply the collection name 'DEFAULT' and the bin id 'EVERYTHING'.")

        if sum(list(requested_split_source.values())) > 1:
            raise ConfigError("You must specify a list of genes (with --gene-caller-ids or\
                               --genes-of-interest), OR a list of splits (--splits-of-interest), OR a\
                               collection and bin combo (--collection-name and bin-id). You\
                               supplied too many of these parameters, and now anvi'o doesn't\
                               know what you want.")

        self.split_source = None
        for source in requested_split_source:
            if requested_split_source[source]:
                self.split_source = source
                break


    def is_available_samples_compatible_with_sample_ids_of_interest(self):
        self.run.info("Samples available", ", ".join(sorted(self.available_sample_ids)), progress=self.progress)
        if self.sample_ids_of_interest:
            samples_missing = [sample_id for sample_id in self.sample_ids_of_interest if sample_id not in self.available_sample_ids]
            if len(samples_missing):
                self.progress.end()
                raise ConfigError('One or more samples you are interested in seem to be missing from\
                                   the %s: %s' % ('variability table' if self.table_provided else 'profile database',
                                                  ', '.join(samples_missing)))

            self.available_sample_ids = sorted(list(self.sample_ids_of_interest))


    def apply_preliminary_filters(self):
        self.run.info('Variability data', '%s entries in %s splits across %s samples'\
                % (pp(len(self.data)), pp(len(self.splits_basic_info)), pp(len(self.available_sample_ids))))

        self.filter_data(criterion = "sample_id",
                         subset_filter = self.sample_ids_of_interest,
                         subset_condition = self.sample_ids_of_interest and self.load_all_samples)

        self.filter_data(criterion = "corresponding_gene_call",
                         subset_filter = self.genes_of_interest,
                         subset_condition = (self.genes_of_interest and self.load_all_genes) or self.only_if_structure)

        self.filter_data(criterion = "split_name",
                         subset_filter = self.splits_of_interest,
                         subset_condition = self.splits_of_interest)

        # let's report the number of positions reported in each sample before filtering any further:
        num_positions_each_sample = dict(self.data.sample_id.value_counts())
        self.run.info('Total number of variable positions in samples', '; '.join(['%s: %s' % (s, num_positions_each_sample.get(s, 0)) for s in sorted(self.available_sample_ids)]))

        self.filter_data(criterion = "departure_from_reference",
                         min_filter = self.min_departure_from_reference,
                         min_condition = self.min_departure_from_reference > 0 and not self.min_departure_from_reference_filtered,
                         max_filter = self.max_departure_from_reference,
                         max_condition = self.max_departure_from_reference < 1 and not self.max_departure_from_reference_filtered)

        if self.engine == 'NT':
            self.data['unique_pos_identifier_str'] = self.data['split_name'] + "_" + self.data['pos'].astype(str)
        elif self.engine == 'CDN' or self.engine == 'AA':
            self.data['unique_pos_identifier_str'] = self.data['split_name'] + "_" + self.data['corresponding_gene_call'].astype(str) + "_" + self.data['codon_order_in_gene'].astype(str)
        else:
            pass

        self.filter_data(function = self.filter_by_occurrence,
                         min_occurrence = self.min_occurrence)

        # this guy has no home
        self.data['gene_length'] = self.data['corresponding_gene_call'].apply(self.get_gene_length)


    def set_unique_pos_identification_numbers(self):
        self.progress.new('Further processing')
        self.progress.update('re-setting unique identifiers to track split/position pairs across samples')

        self.data['unique_pos_identifier'] = self.data['unique_pos_identifier_str'].apply(self.get_unique_pos_identification_number)
        self.data['contig_name'] = self.data['split_name'].apply(lambda split: self.splits_basic_info[split]['parent'])

        self.progress.end()


    def get_unique_pos_identification_number(self, unique_pos_identifier_str):
        if unique_pos_identifier_str in self.split_name_position_dict:
            return self.split_name_position_dict[unique_pos_identifier_str]
        else:
            self.split_name_position_dict[unique_pos_identifier_str] = self.unique_pos_identifier
            self.unique_pos_identifier += 1
            return self.split_name_position_dict[unique_pos_identifier_str]


    def gen_unique_pos_identifier_to_entry_id_dict(self):
        self.progress.new('Generating the `unique pos` -> `entry id` dict')
        self.progress.update('...')

        self.unique_pos_id_to_entry_id = {}

        for entry_id in self.data:
            v = self.data[entry_id]
            u = v['unique_pos_identifier_str']
            if u in self.unique_pos_id_to_entry_id:
                self.unique_pos_id_to_entry_id[u].add(entry_id)
            else:
                self.unique_pos_id_to_entry_id[u] = set([entry_id])

        self.progress.end()


    def compute_additional_fields(self, entry_ids=[]):
        """
        This calculates the following columns: consensus, n2n1 ratio, competing_aas,
        departure_from_consensus, and all substitution scoring matrices.

        NOTE For defining the "consensus" column and the "competing_aas" column (or related column
        for different --engine values), it is important to make it explict how we resolve ties. If
        you finish reading this comment and still do not understand, than I have failed you. It's
        an issue because suppose Ala and Trp are tied for sharing the most reads. Is the consensus
        Ala or Trp? Is competing_aas AlaTrp or TrpAla? In a separate example, if Ser is most
        common and Gly and Glu are tied for second, should Gly or Glu be a part of competing_aas
        and should it be GlxSer or SerGlx? There are three rules that define our conventions:

            1. Competing_aas ALWAYS appear in alphabetical order. Even if Cys is most common, and
               Ala is second most commond, competing_aas = AlaCys.
            2. Ties are always resolved alphabetically. If there is a 3-way tie for second between
               His, Met, and Thr, the item included in competing_aas will be His.
            3. If the coverage of the second-most common item is 0, the most common is paired with
               itself.

        NOTE if the coverage is 0 (which is stupid, but it can exist if both
        --min-coverage-in-each-sample is 0 and --quince-mode is active), departure_from_consensus
        and n2n1ratio are set to 0.

        NOTE To get a first glance at variation, we display competing_nts during inspect mode in the
        interactive interface, which means they are stored. I (Evan) don't know how the
        competing_nts are computed during this process, but after looking at 3 datasets (Infant Gut,
        E. Faecalis; TARA Oceans, HIMB083; and Mushroom Spring, Synechococcus), there were no "N"
        counts (--engine NT), which makes me think that consensus and competing_nts are calculated
        without consideration of "N" values. In contrast, the variability table is not prejudiced to
        N, and treats it like any other item in self.items. This means it can be the consensus
        value, or be one of the items in competing_nts, etc.
        """
        self.progress.new("Computing additional fields")
        self.progress.update("...")

        # First, we just make sure that whatever operations we have performed on self.data, we have
        # not altered the types of coverage values from int. I am learning that with pandas
        # DataFrames, it is WORTH being explicit with types, even at the cost of redundancy.
        self.data.loc[:, self.items] = self.data.loc[:, self.items].astype(int)

        if not len(entry_ids):
            entry_ids = list(self.data.index)

        # index entries with and without non-zero coverage (introduced by --quince-mode)
        coverage_zero = self.data.index.isin(entry_ids) & (self.data["coverage"] == 0)
        coverage_nonzero = self.data.index.isin(entry_ids) & (self.data["coverage"] > 0)

        # rank the items of each entry from 1 - 21 (for --engine AA) based on item coverage.
        # method="first" ensures alphabetic ordering in the case of ties. Convert the rank DataFrame
        # into a numpy array, and find the order of indices that sort each entry's items based on
        # their ranks. type(ranks) = pd.DataFrame, type(item_index_order) = numpy array
        ranks = self.data.loc[entry_ids, self.items].rank(ascending=False, axis=1, method="first").astype(int)
        item_index_order = np.argsort(ranks.values, axis=1)

        # the first and second most common items, according to the 2nd convention in the docstring,
        # are now defined for each entry in two pandas Series.
        items_first_and_second = self.data.loc[entry_ids, self.items].columns[item_index_order[:,:2]]

        # we also calculate the coverage values for the first and second most common items
        sorted_coverage = np.sort(self.data.loc[entry_ids, self.items].values, axis=1)
        coverages_first =  pd.Series(sorted_coverage[:,-1], index=self.data.loc[entry_ids].index)
        coverages_second = pd.Series(sorted_coverage[:,-2], index=self.data.loc[entry_ids].index)

        # define consensus as the first most common item (hence `[:,0]`)
        self.data.loc[entry_ids, "consensus"] = items_first_and_second[:,0]

        # if the coverage is zero, departure_from_consensus = 0
        self.data.loc[coverage_zero, "departure_from_consensus"] = 0
        self.data.loc[coverage_nonzero, "departure_from_consensus"] = \
                    (self.data.loc[coverage_nonzero, "coverage"] - coverages_first) / self.data.loc[coverage_nonzero, "coverage"]

        # if the coverage is zero, n2n1ratio = 0
        self.data.loc[coverage_zero, "n2n1ratio"] = 0
        self.data.loc[coverage_nonzero, "n2n1ratio"] = coverages_second[coverage_nonzero] / coverages_first[coverage_nonzero]

        # we name the competing_items column differently depending on the engine used
        if self.engine == 'NT':
            self.competing_items = 'competing_nts'
        elif self.engine == 'CDN':
            self.competing_items = 'competing_codons'
        elif self.engine == 'AA':
            self.competing_items = 'competing_aas'
        else:
            pass

        # This step computes the "competing_items" column. First, convert items_first_and_second
        # [which has shape (len(entry_ids), 2)] into a numpy array, then sort them alphabetically
        # (in order to comply with convention 1 in docstring). Then, sum along the sorted axis.
        # Since the elements are of type str, the sum operator concatenates the strings together.
        # Finally, if the coverage of the most common item is equal to the total coverage, we pair
        # the most common item with itself.
        self.data.loc[entry_ids, self.competing_items] = np.sum(np.sort(items_first_and_second.values, axis=1), axis=1) # V/\
        self.data.loc[self.data.index.isin(entry_ids) & (self.data["coverage"] == self.data[self.items].max(axis=1)), self.competing_items] = self.data["consensus"]*2

        # Loop through each SSM, filling each corresponding column entry by entry using the `apply`
        # operator. Instead of using self.substitution_scoring_matrices[m], we speed things up by
        # converting to a dictionary we call `substitution_scoring_matrix` that takes as its input
        # an entry from competing_items instead of having to pull from the first AND second most
        # common item, and then indexing a triple-nested dictionary.
        for m in self.substitution_scoring_matrices:
            substitution_scoring_matrix = utils.convert_SSM_to_single_accession(self.substitution_scoring_matrices[m])
            self.data.loc[entry_ids, m] = self.data.loc[entry_ids, self.competing_items].apply(lambda x: substitution_scoring_matrix.get(x, None))

        self.progress.end()


    def report_change_in_entry_number(self, num_before, num_after, reason="unknown reason"):
        """Reports how many entries were removed (or added) during a filtering step."""
        changed = "removed" if num_after <= num_before else "added"

        genes_remaining = self.data["corresponding_gene_call"].unique()
        if self.append_structure_residue_info:
            structures_remaining = [gene for gene in self.genes_with_structure if gene in genes_remaining]

        extra_msg = ", %d with structure" % (len(structures_remaining)) if self.append_structure_residue_info else ""
        self.run.info('Entries after %s filter' % (reason),
                      '%s (%s were %s) [%s genes remaining%s]' % (pp(num_after),
                                                                  pp(abs(num_before - num_after)),
                                                                  changed,
                                                                  len(genes_remaining),
                                                                  extra_msg),
                      mc='green')

        self.check_if_data_is_empty()


    def compute_comprehensive_variability_scores(self):
        """
            Comprehensive stats are defined as scores that take into consideration the entire vector of variability and
            not only the two most competing items (and thus it is comprehensive).  Currently the scores that are
            included are: site-entropy, site-Kullback-Leibler divergence (both a raw score and a normalized score (see
            below)), and weighted substitution scores (i.e. BLOSUM).

            site-entropy
            ============
            The entropy of the items at a single position in a single sample. This means that each entry of the
            variability table receives its own site-entropy score (in the table we just called it "entropy"). If a site
            has no coverage, which can happen if --quince-mode is enabled, the value of entropy is -np.inf, something we
            maintain until exporting the table as a tab-delimited file, at which point we recast them to something
            reasonable.

            Kullback-Leibler divergence raw
            ===============================
            the Kullback-Leibler divergence of the frequencies in a sample compared to the raw frequencies of the sum of
            occurrences in the same site accross samples.

            Kullback-Leibler divergence normalized
            ======================================
            The Kullback-Leibler divergence of the frequencies in a sample compared to the frequencies of the sum of
            normalized occurances in the same site accross samples. Where the normalization is such that the occernce of
            items is normalized to sum to one in each sample. This method eliminates the effect of coverage on the
            score. The disadvantage of this method is that if there is a sample with low coverage then any noise (like a
            single sequencing error) could have a major effect. It is recommended to use this score in combination with
            the --min-coverage-in-each-sample.

            Weighted substitution scores
            ============================
            The weights per substitution score is weighted by the joint frequency of the items i.e. sum(S_{i,j}*pi*pj)
            where i does not equal j (i.e. the substitution of an item with itself is not considered)
        """

        if self.skip_comprehensive_variability_scores:
            self.run.warning("Anvi'o will skip comprehensive variability score computations.")
            return

        # suppress division by zero runtime warning
        np.seterr(invalid='ignore')

        if not self.quince_mode:
            self.run.warning("Only some comprehensive variability metrics can be computed without `--quince-mode`")

        self.progress.new("Comprehensive stats")
        self.progress.update("Those that don't require --quince-mode")

        self.comprehensive_stats_headers = [m + '_weighted' for m in self.substitution_scoring_matrices] + ['entropy']

        # Pandas is fun, but numpy is fast. Here we convert the coverage table information from the DataFrame to a
        # numpy array. The transpose is required because scipy.stats entropy function calculates along an
        # unspecifiable axis that we must conform to. But before any of this is done we order the entries according
        # to unique_pos_identifier (and for a given unique_pos_identifier, entries are ordered alphabetically by
        # sample_id). The reason for this is aesthetic but also required for vectorized operations that occur after
        # self.progress.update("Those that do require --quince-mode")
        self.data = self.data.sort_values(by=["unique_pos_identifier", "sample_id"])
        coverage_table = self.data[self.items].T.astype(int).values

        # Now we compute the entropy, which is defined at a per position, per sample basis. There is a reason we
        # pass coverage_table instead of a normalized table. If we pass a normalized table scipy.stat.entropy complains
        # that there is division by zero for entries introduced by --quince-mode that have 0 coverage.  By passing
        # coverage_table instead, scipy.stats.entropy does the normalization itself ensures that such entries return
        # -inf instead of raising an error.
        self.data["entropy"] = entropy(coverage_table)

        # Next we worry ourselves with weighted substitution scores. We convert each SSM to a numpy array and store
        # each of them in a dictionary with keys equal to the SSM names, e.g. BLOSUM90. The index of numpy arrays
        # are organized alphabetically. For example, for an amino acid substitution matrix, the substitution
        # Ala->Ala is indexed by [0,0], whereas Val->Val is indexed by [-1,-1]. We decided it makes sense to set the
        # substitution score of an item with itself to zero. This way we only consider substitutions to other items
        # (and don't consider substitution of an item with itself)
        numpy_substitution_scoring_matrices = {SSM: np.array([[matrix[i][j] if j != i else 0 for j in sorted(matrix[i])] for i in sorted(matrix)]) for SSM, matrix in self.substitution_scoring_matrices.items()}

        # We loop through each of the substitution scoring matrices available. It's possible the items in the
        # substitution matrices aren't a complete set of the items we have coverage data for. For example, BLOSUM90
        # doesn't have substitution scores for the STP codon (yet we have coverage data for STP). Available items
        # can vary from SSM to SSM. To deal with this, we subset the coverage_table to only include available items
        # for a given SSM. We then transform this subset of the coverage data to a frequency table that's
        # potentially unique to the SSM. To avoid calculating frequency tables unnecessarily, the previous SSM's
        # items are stored in `previous_indices` so the table is only recalculated if the current SSM items vary
        # from the previous SSM items.
        previous_indices = None
        for SSM, matrix in numpy_substitution_scoring_matrices.items():

            # initialize the entropy column in self.data with np.nans
            self.data[SSM + "_weighted"] = np.nan

            # We find which items in self.items are in the SSM and record the index at which they appear in
            # self.items. By definition this is also the index ordering of coverage_table. For example, the index of
            # Ala is 0 so the coverage data for Ala is coverage_table[0,:]. Hence, the subset of coverage_table for
            # the given SSM is coverage_table[indices, :].
            indices = sorted([self.items.index(item) for item in self.substitution_scoring_matrices[SSM]])
            if indices != previous_indices:

                # To get the frequency table we divide each column by the sum of the rows (the array called
                # total_coverage). But we must be careful, because it's possible for some entries of total_coverage
                # to be 0, which can occur if both a) --quince-mode is enabled and b) --min-coverage-in-each-sample
                # is 0 (the default). If this is the case we set the total_coverage entry to -1 so that the
                # frequency for each item in the entry becomes 0/-1 = 0, instead of producing a NaN due to division
                # by zero.
                total_coverage = np.sum(coverage_table[indices,:], axis=0)
                freq_table = np.divide(coverage_table[indices,:], np.where(total_coverage==0, -1, total_coverage))

                # While in this loop, we define a Boolean array with length equal to the number of entries in
                # freq_table. It is True if the substitution score should be calculated and False otherwise. What
                # can cause an entry not to be calculated? First of all, if there is no coverage data then there is
                # no substitution to report. Secondly, since we don't consider the substitution of an item with
                # itself (Remember? We set all diagonals in each SSM to 0), we don't report scores for entries that
                # only have 1 item with non-zero coverage. Both of these conditions are encapsulated nicely with the
                # np.count_nonzero function. We then subset the freq_table to only include these entries but we hang
                # onto the Boolean array for when we put the scores into self.data
                keep = np.count_nonzero(freq_table, axis=0) > 1
                freq_table = freq_table[:, keep]

                # The last thing we do in here is calculate a normalization factor for the entries, since it is
                # unique to each freq_table. A normalization score is needed since we don't consider a substitution
                # of an item with itself. Hence, the sum of frequencies doesn't sum to 1, and so to make sure they
                # sum to 1 we multiply by this normalization factor. We want these to sum to 1 because otherwise
                # these are not valid weights.
                normalization = 1 / (1 - np.sum(np.square(freq_table), axis=0))

            # This is legitimately legendary status array broadcasting. I can't explain how it works but the
            # quantity sum(S_{i,j}*pi*pj) is being calculated for each entry, where pi is the frequency of the ith
            # item, pj is the frequency of the jth item, and S_{i,j} is the substitution matrix score for the i->j
            # substitution. If you remember, we already set i=j entries to zero so they contribute zero weight to
            # the score.
            self.data.loc[keep, SSM + "_weighted"] = normalization * np.sum(freq_table[np.newaxis, :].T * freq_table[:, np.newaxis].T * matrix, axis=(1,2)) # V/\

        if not self.quince_mode:
            self.progress.end()
            self.comprehensive_variability_scores_computed = True
            return

        self.progress.update("Those that do require --quince-mode")
        self.comprehensive_stats_headers.extend(['kullback_leibler_divergence_raw', 'kullback_leibler_divergence_normalized'])

        # Due to --quince-mode every unique position identifier has the same number of entries (this screams
        # vectorization). We abuse this to make a 3 dimensional numpy array, coverage_by_pos. The zeroth axis
        # indexes the items, the first indexes the unique_pos_identifier, and the second indexes the sample number.
        # This reshaping operation works only because at the start of this function we ordered entries by
        # unique_pos_identifier.
        numpy_dimension = (len(self.items), self.data["unique_pos_identifier"].nunique(), len(self.available_sample_ids))
        coverage_by_pos = coverage_table.reshape(numpy_dimension)

        # We also define a normalized version of coverage_by_pos so that each entries defines the frequency
        # (probability) of a certain item occuring, rather than the raw counts. This is used for the normalized
        # Kullback-Leibler divergence. If the entry is all zeros, the frequencies for the entry are all defined to
        # be 0. NOTE If this is not done, any position where one or more of the samples has 0 coverage yields a
        # normalized kullback-leibler value of inf. I think it makes more sense this way.
        counts_per_pos_per_sample = np.sum(coverage_by_pos, axis=0)
        freq_by_pos = np.divide(coverage_by_pos, counts_per_pos_per_sample)

        # The entropy from scipy.stats operates on axis = 0, so the returned array dimension is flattened in the
        # zeroth dimension compared to the input array. For example, if the input is shape (X,Y,Z), the output is
        # shape (Y,Z). That means when we pass coverage_by_pos, an entropy array is returned with the zeroth axis
        # indexed by unique_pos_identifiers and the first axis indexed by sample_ids. The entropy function insists
        # that our reference distributions (the normalized and unnormalized mean frequency counts across samples,
        # see docstring for more details) must ALSO have the shape (X,Y,Z). The issue with this is that by
        # calculating the mean over samples we collapse the Z dimension so the shape of our reference distribution
        # array is (X,Y). We solve this issue by stacking Z identical reference distribution arrays to create a
        # pseudo-second axis so the final shape is (X,Y,Z). coverage_summed_over_samples is the reference
        # distribution array for the raw Kullback-Leibler divergence and freq_summed_over_samples is for the
        # normalized Kullback-Leibler divergence.
        coverage_summed_over_samples = np.repeat(np.sum(coverage_by_pos, axis=2)[:,:,np.newaxis], len(self.available_sample_ids), axis=2)
        freq_summed_over_samples     = np.repeat(np.sum(freq_by_pos,     axis=2)[:,:,np.newaxis], len(self.available_sample_ids), axis=2)

        # As mentioned in last comment, a 2D array is returned by entropy, which we flatten into 1D.  The flattened
        # array is equal to the length of entries in self.data. Furthermore, the order of entropy values is the same
        # order as the entries they correspond to, so all we do is assign the arrays to new columns in self.data and
        # we're done.
        self.data['kullback_leibler_divergence_raw'] = entropy(coverage_by_pos, coverage_summed_over_samples).flatten()
        self.data['kullback_leibler_divergence_normalized'] = entropy(freq_by_pos, freq_summed_over_samples).flatten()

        self.progress.end()
        self.comprehensive_variability_scores_computed = True


    def get_unique_positions_and_frequencies_dict(self):
        """From the self.data object, creates a dict that contains item frequencies for
           each sample for each unique position identifier."""

        unique_positions_and_frequencies_dict = {}

        template = dict.fromkeys(self.items, 0)

        if not self.unique_pos_id_to_entry_id:
            self.gen_unique_pos_identifier_to_entry_id_dict()

        self.progress.new('The unique positions and frequencies dict')
        self.progress.update('generating ..')

        for entry_ids in self.unique_pos_id_to_entry_id.values():
            unique_pos_identifier = self.data[list(entry_ids)[0]]['unique_pos_identifier_str']
            unique_positions_and_frequencies_dict[unique_pos_identifier] = {}

            for entry_id in entry_ids:
                v = self.data[entry_id]
                unique_positions_and_frequencies_dict[unique_pos_identifier][v['sample_id']] = copy.deepcopy(template)

                for item in self.items:
                    if v[item]:
                        unique_positions_and_frequencies_dict[unique_pos_identifier][v['sample_id']][item] = v[item]

        self.progress.end()

        return unique_positions_and_frequencies_dict


    def process(self, process_functions=None, exit_if_data_empty=True):
        """self.data is checked if empty after each function call. if exit_if_data_empty, exists,
           otherwise returns prematurely."""
        if not process_functions:
            process_functions = self.process_functions

        try:
            for func, kwargs in process_functions:
                func(**kwargs)

        except self.EndProcess as e:
            msg = 'Nothing left in the variability data to work with. Quitting :/' if exit_if_data_empty else ''
            e.end(exit_if_data_empty, msg)


    def get_histogram(self, column, fix_offset=False, **kwargs):
        """fix_offset can be provided if you're interested in returning the centre point of each bin
           rather than the edges of each bin.

           **kwargs are the optional arguments of np.histogram
           (https://docs.scipy.org/doc/numpy-1.14.0/reference/generated/numpy.histogram.html)
        """
        if not pd.api.types.is_numeric_dtype(self.data[column]):
            raise ConfigError("get_histogram :: %s is not of numeric type" % (column))

        if fix_offset:
            range_offset = (kwargs["range"][1] - kwargs["range"][0]) / (kwargs["bins"] - 1) / 2
            kwargs["range"] = (kwargs["range"][0] - range_offset, kwargs["range"][1] + range_offset)

        # define numpy array; filter infinities and nans
        column_data = self.data[column].values
        column_data = column_data[np.isfinite(column_data)]

        # histogram
        values, bins = np.histogram(self.data[column], **kwargs)

        if fix_offset:
            bins = bins[:-1] + range_offset
            # now bins have the same length as values and represent the midpoint of each bin (e.g.
            # the first bin value is the original minimum value passed to this function

        return values, bins


    def get_residue_structure_information(self):
        """
        If by the end of all filtering there is no overlap between genes with variability and genes
        with structure, this function raises a warning and the structure columns are not added to
        the table. Otherwise this function appends the columns from residue_info to self.data
        """
        if not self.append_structure_residue_info:
            return

        genes_with_var = list(self.data["corresponding_gene_call"].unique())
        genes_with_var_and_struct = [g for g in self.genes_with_structure if g in genes_with_var]
        if not genes_with_var_and_struct:
            run.warning("Before filtering entries, there was an overlap between genes with\
                         variability and genes with structuresr, however that is no longer the case. As a result,\
                         no structure information columns will be added. Above you can see the number of genes\
                         remaining with structures after each filtering step.")
            self.append_structure_residue_info = False
            return

        self.progress.new("Adding structure db info")
        self.progress.update("Appending residue_info columns")

        skip_columns = ["entry_id"]
        include_columns = [x for x in self.structure_residue_info.columns if x not in skip_columns]
        self.data = pd.merge(self.data,
                             self.structure_residue_info[include_columns],
                             on = ["corresponding_gene_call", "codon_order_in_gene", "codon_number"],
                             how = "left")

        # add all known residue info sources to columns_to_report
        C = {'text': str, 'real': float, 'integer': int}
        redundant_columns = ['entry_id', 'corresponding_gene_call', 'codon_order_in_gene', 'aa', 'amino_acid', 'codon', 'codon_number']

        # append mandatory columns
        self.columns_to_report['structural'].extend([(x, C[y]) for x, y in zip(t.structure_residue_info_table_structure, t.structure_residue_info_table_types) if x not in redundant_columns])

        # append non-mandatory columns
        for source in t.residue_info_sources:
            self.columns_to_report['structural'].extend([(x, C[y]) for x, y in zip(t.residue_info_sources[source]["structure"], t.residue_info_sources[source]["types"]) if x not in redundant_columns])

        self.progress.end()


    def compute_gene_coverage_fields(self):
        if not self.compute_gene_coverage_stats:
            return

        # Initialize the profile super FIXME This bastard spits out
        #       Auxiliary Data ...............................: Found: SAR11/AUXILIARY-DATA.db (v. 2)
        #       Profile Super ................................: Initialized with all 1393 splits: SAR11/PROFILE.db (v. 27)
        # and it isn't silenced even if self.Run(verbose=False) is passed to the VariabilitySuper class
        profile_super = dbops.ProfileSuperclass(argparse.Namespace(profile_db = self.profile_db_path))

        self.progress.new('Computing gene coverage stats')
        self.progress.update('... {consider --skip-gene-coverage-stats if taking too long}')

        # obtain gene coverage info per gene/sample combo
        gene_cov_dict = {}
        for split_name in self.splits_of_interest:
            entry_ids = self.split_name_to_genes_in_splits_entry_ids[split_name]
            for entry_id in entry_ids:
                gene_cov_dict.update(profile_super.get_gene_level_coverage_stats(
                    self.genes_in_splits[entry_id]['split'],
                    self,
                    gene_caller_ids_of_interest = set([self.genes_in_splits[entry_id]['gene_callers_id']])
                ))

        gene_coverage_columns = ['gene_coverage',
                                 'non_outlier_gene_coverage',
                                 'non_outlier_gene_coverage_std']

        J = lambda row, g: (g[row.iloc[0]][row.iloc[1]]['mean_coverage'],
                            g[row.iloc[0]][row.iloc[1]]['non_outlier_mean_coverage'],
                            g[row.iloc[0]][row.iloc[1]]['non_outlier_coverage_std']) \
                            if row.iloc[0] in self.gene_lengths else (-1, -1, -1)

        self.data = utils.apply_and_concat(df = self.data,
                                           fields = ['corresponding_gene_call', 'sample_id'],
                                           func = J,
                                           column_names = gene_coverage_columns,
                                           func_args = (gene_cov_dict,))

        # this guy piggybacks in this method
        self.data['mean_normalized_coverage'] = self.data['coverage'] / self.data['gene_coverage']
        self.data.loc[self.data['gene_coverage'] == -1, 'mean_normalized_coverage'] = -1

        self.progress.end()


    def get_gene_length(self, gene_callers_id):
        if gene_callers_id in self.gene_lengths:
            return self.gene_lengths[gene_callers_id]
        else:
            return -1


    def get_unique_pos_identifier_to_corresponding_gene_id(self):
        self.progress.update('populating a dict to track corresponding gene ids for each unique position')

        # key = unique_pos_identifier, val = corresponding_gene_call
        return self.data[["unique_pos_identifier","corresponding_gene_call"]].\
               drop_duplicates().set_index("unique_pos_identifier").to_dict()["corresponding_gene_call"]


    def get_unique_pos_identifier_to_codon_order_in_gene(self):
        self.progress.update('populating a dict to track codon order in genes for each unique position')
        # key = unique_pos_identifier, val = codon_order_in_gene
        return self.data[["unique_pos_identifier","codon_order_in_gene"]].\
               drop_duplicates().set_index("unique_pos_identifier").to_dict()["codon_order_in_gene"]


    def report(self, data=None):
        if data is None:
            data = self.data

        self.progress.new('Reporting variability data')

        new_structure, _ = self.get_data_column_structure()

        if not self.include_contig_names_in_output:
            new_structure.remove('contig_name')

        if not self.include_split_names_in_output:
            new_structure.remove('split_name')

        # Update entry_id with sequential numbers based on the final ordering of the data:
        data.reset_index(drop=True, inplace=True)
        data["entry_id"] = data.index

        # order by [corresponding_gene_call, codon_order_in_gene]
        data = data.sort_values(by = ["corresponding_gene_call", "codon_order_in_gene"])

        self.progress.update('exporting variable positions table as a TAB-delimited file ...')
        utils.store_dataframe_as_TAB_delimited_file(data, self.args.output_file, columns=new_structure)
        self.progress.end()

        self.run.info('Num entries reported', pp(len(data.index)))
        self.run.info('Output File', self.output_file_path)
        self.run.info('Num %s positions reported' % self.engine, data["unique_pos_identifier"].nunique())


    def get_data_column_structure(self, data=None):
        if data is None:
            data = self.data

        structure = []
        data_types = []
        for column_group, columns in self.columns_to_report.items():
            for column, data_type in columns:
                if column in data.columns:
                    structure.append(column)
                    data_types.append(data_type)
        return structure, data_types


    class EndProcess(Exception):
        def end(self, exit, msg=None):
            """exit: bool
                   if True, sys.exit() is called.
               msg: str (default=None)
                   message to user
            """
            if msg:
                run.info_single(msg, 'red', 1, 1)
            if exit:
                sys.exit()


class NucleotidesEngine(dbops.ContigsSuperclass, VariabilitySuper):
    """This is the main class to make sense and report variability for a given set of splits,
       or a bin in a collection, across multiple or all samples. The user can scrutinize the
       nature of the variable positions to be reported dramatically given the ecology and/or
       other biologically-relevant considerations, or purely technical limitations such as
       minimum coverage of a given nucleotide position or the ratio of the competing nts at a
       given position. The default entry to this class is the `anvi-gen-variability-profile`
       program."""

    def __init__(self, args={}, p=progress, r=run):
        self.run = r
        self.progress = p

        self.engine = 'NT'

        # Init Meta
        VariabilitySuper.__init__(self, args=args, r=self.run, p=self.progress)


    def recover_base_frequencies_for_all_samples(self):
        """this function populates variable_nts_table dict with entries from samples that have no
           variation at nucleotide positions reported in the table"""
        if not self.quince_mode:
            return

        self.progress.new('Recovering NT data')

        samples_wanted = self.sample_ids_of_interest if self.sample_ids_of_interest else self.available_sample_ids
        splits_wanted = self.splits_of_interest if self.splits_of_interest else set(self.splits_basic_info.keys())
        next_available_entry_id = self.data["entry_id"].max() + 1

        unique_pos_identifier_to_corresponding_gene_id = self.get_unique_pos_identifier_to_corresponding_gene_id()

        unique_pos_identifier_to_codon_order_in_gene = self.get_unique_pos_identifier_to_codon_order_in_gene()
        self.progress.update('creating a dict to track missing base frequencies for each sample / split / pos')
        split_pos_to_unique_pos_identifier = {}
        splits_to_consider_dict = {}
        for split_name in splits_wanted:
            splits_to_consider_dict[split_name] = {}
            split_pos_to_unique_pos_identifier[split_name] = {}

        self.progress.update('populating the dict to track missing base frequencies for each sample / split / pos')
        for entry_id, v in self.data.iterrows():
            p = v['pos']
            d = splits_to_consider_dict[v['split_name']]
            u = split_pos_to_unique_pos_identifier[v['split_name']]

            if p in d:
                d[p].remove(v['sample_id'])
            else:
                d[p] = copy.deepcopy(samples_wanted)
                d[p].remove(v['sample_id'])

            if p not in u:
                u[p] = v['unique_pos_identifier']

        split_names_to_consider = list(splits_to_consider_dict.keys())
        num_splits = len(split_names_to_consider)
        new_entries = {}
        for split_index in range(num_splits):
            split = split_names_to_consider[split_index]
            self.progress.update('Accessing split covs, updating variable pos dict (%s of %s)' % (pp(split_index + 1), pp(num_splits)))

            split_coverage_across_samples = self.merged_split_coverage_values.get(split)

            split_info = self.splits_basic_info[split]
            contig_name_name = split_info['parent']

            for pos in splits_to_consider_dict[split]:
                unique_pos_identifier = split_pos_to_unique_pos_identifier[split][pos]
                contig_name_seq = self.contig_sequences[contig_name_name]['sequence']
                pos_in_contig = split_info['start'] + pos
                base_at_pos = contig_name_seq[pos_in_contig]
                corresponding_gene_call = unique_pos_identifier_to_corresponding_gene_id[unique_pos_identifier]
                gene_length = self.get_gene_length(corresponding_gene_call)
                codon_order_in_gene = unique_pos_identifier_to_codon_order_in_gene[unique_pos_identifier]

                in_partial_gene_call, in_complete_gene_call, base_pos_in_codon = self.get_nt_position_info(contig_name_name, pos_in_contig)

                for sample_id in splits_to_consider_dict[split][pos]:
                    new_entries[next_available_entry_id] = {'entry_id': next_available_entry_id,
                                                            'contig_name': contig_name_name,
                                                            'departure_from_reference': 0,
                                                            'reference': base_at_pos,
                                                            'A': 0, 'T': 0, 'C': 0, 'G': 0, 'N': 0,
                                                            'pos': pos,
                                                            'pos_in_contig': pos_in_contig,
                                                            'in_partial_gene_call': in_partial_gene_call,
                                                            'in_complete_gene_call': in_complete_gene_call,
                                                            'base_pos_in_codon': base_pos_in_codon,
                                                            'coverage': split_coverage_across_samples[sample_id][pos],
                                                            'sample_id': sample_id,
                                                            'cov_outlier_in_split': 0,
                                                            'cov_outlier_in_contig': 0,
                                                            'competing_nts': base_at_pos + base_at_pos,
                                                            'unique_pos_identifier': unique_pos_identifier,
                                                            'unique_pos_identifier_str': '%s_%d' % (split, pos),
                                                            'corresponding_gene_call': corresponding_gene_call,
                                                            'gene_length': gene_length,
                                                            'codon_order_in_gene': codon_order_in_gene,
                                                            'split_name': split}
                    new_entries[next_available_entry_id][base_at_pos] = split_coverage_across_samples[sample_id][pos]
                    next_available_entry_id += 1

        # convert to pandas DataFrame (its faster to build and convert a dictionary than to build
        # DataFrame row by row).
        new_entries = pd.DataFrame(new_entries).T

        # before concatenating the new entries, store the self.data column order. Also, check that
        # no columns exist in new_entries but not in self.data. This is unacceptable, and could have
        # happened if code for new_entries was changed or if the workflow in process() is
        # significantly reworked.
        column_order = self.data.columns.tolist()
        if len([x for x in new_entries.columns.tolist() if x not in self.data.columns.tolist()]):
            raise ValueError("Columns found in new_entries exist that aren't in self.data.")

        # concatenate new columns to self.data
        entries_before = len(self.data.index)
        self.data = pd.concat([self.data, new_entries], sort=True)
        new_entries.set_index("entry_id", drop=False, inplace=True)
        self.data = self.data[column_order]
        entries_after = len(self.data.index)

        # fill in additional fields for new entries. compute_additional_fields takes a list of
        # entry_ids to consider for self.data, which here is provided from new_entries (what I'm
        # saying is new_entries is not passed, only the entry_id's in new_entries
        self.progress.end()
        self.compute_additional_fields(list(new_entries["entry_id"]))

        self.report_change_in_entry_number(entries_before, entries_after, reason="quince mode")


class QuinceModeWrapperForFancyEngines(object):
    """A base class to recover quince mode data for both CDN and AA engines.

       This wrapper exists outside of the actual classes for these engines since
       the way they recover these frequencies is pretty much identical except one
       place where the engine needs to be specifically.
    """
    def __init__(self):
        if self.engine not in ['CDN', 'AA']:
            raise ConfigError("This fancy class is only relevant to be inherited from within CDN or AA engines :(")


    def recover_base_frequencies_for_all_samples(self):
        if not self.quince_mode:
            return

        self.progress.new('[%s] Recovering item variability data' % self.engine)

        samples_wanted = self.sample_ids_of_interest if self.sample_ids_of_interest else self.available_sample_ids
        splits_wanted = self.splits_of_interest if self.splits_of_interest else set(self.splits_basic_info.keys())
        next_available_entry_id = self.data["entry_id"].max() + 1

        unique_pos_identifier_str_to_consenus_item = {}
        unique_pos_identifier_str_to_unique_pos_identifier = {}
        for _, e in self.data.iterrows():
            upi = e['unique_pos_identifier_str']
            unique_pos_identifier_str_to_consenus_item[upi] = e['reference']
            unique_pos_identifier_str_to_unique_pos_identifier[upi] = e['unique_pos_identifier']

        self.progress.update('creating a dict to track missing item frequencies for each sample / split / pos')

        splits_to_consider_dict = {}
        for split_name in splits_wanted:
            splits_to_consider_dict[split_name] = {}

        self.progress.update('populating the dict to track missing item frequencies for each sample / split / pos')
        for entry_id, v in self.data.iterrows():
            gene_codon_key = '%d_%d' % (v['corresponding_gene_call'], v['codon_order_in_gene'])
            d = splits_to_consider_dict[v['split_name']]

            if gene_codon_key in d:
                d[gene_codon_key].remove(v['sample_id'])
            else:
                d[gene_codon_key] = copy.deepcopy(samples_wanted)
                d[gene_codon_key].remove(v['sample_id'])

        split_names_to_consider = list(splits_to_consider_dict.keys())
        num_splits = len(split_names_to_consider)
        new_entries = {}
        for split_index in range(num_splits):
            split_name = split_names_to_consider[split_index]
            self.progress.update('Accessing split covs, updating variable pos dict (%s of %s)' % (pp(split_index + 1), pp(num_splits)))

            split_coverage_across_samples = self.merged_split_coverage_values.get(split_name)

            split_info = self.splits_basic_info[split_name]
            contig_name = split_info['parent']

            for gene_codon_key in splits_to_consider_dict[split_name]:
                corresponding_gene_call, codon_order_in_gene = [int(k) for k in gene_codon_key.split('_')]
                gene_length = self.get_gene_length(corresponding_gene_call)

                for sample_name in splits_to_consider_dict[split_name][gene_codon_key]:
                    unique_pos_identifier_str = '_'.join([split_name, str(corresponding_gene_call), str(codon_order_in_gene)])
                    reference_item = unique_pos_identifier_str_to_consenus_item[unique_pos_identifier_str]

                    new_entries[next_available_entry_id] = {'entry_id': next_available_entry_id,
                                                            'unique_pos_identifier_str': unique_pos_identifier_str,
                                                            'unique_pos_identifier': unique_pos_identifier_str_to_unique_pos_identifier[unique_pos_identifier_str],
                                                            'sample_id': sample_name,
                                                            'split_name': split_name,
                                                            'contig_name': contig_name,
                                                            'corresponding_gene_call': corresponding_gene_call,
                                                            'gene_length': gene_length,
                                                            'codon_order_in_gene': codon_order_in_gene,
                                                            'departure_from_reference': 0,
                                                            'coverage': None,
                                                            'reference': reference_item}

                    # DEALING WITH COVERAGE ##################################################################
                    # some very cool but expensive shit is going on here, let me break it down for poor souls of the future.
                    # what we want to do is to learn the coverage of this codon or amino acid in the sample. all we have is
                    # the corresponding gene call id, and the order of this codon or amino acid in the gene. so here how it goes:
                    #
                    # learn the gene call
                    gene_call = self.genes_in_contigs_dict[corresponding_gene_call]

                    # the following dict converts codon  orders into nt positions in contig for a geven gene call
                    codon_order_to_nt_positions_in_contig = utils.get_codon_order_to_nt_positions_dict(gene_call)

                    # so the nucleotide positions for this codon or amino acid in the contig is the following:
                    nt_positions_for_codon_in_contig = codon_order_to_nt_positions_in_contig[codon_order_in_gene]

                    # but we need to convert those positions to the context of this split. so here is the start pos:
                    split_start = self.splits_basic_info[split_name]['start']

                    # here we map nt positions from the contig context to split context using the start position
                    nt_positions_for_codon_in_split = [p - split_start for p in nt_positions_for_codon_in_contig]

                    # we acquire coverages that match to these positions
                    coverages = split_coverage_across_samples[sample_name][nt_positions_for_codon_in_split]
                    coverage = int(round(sum(coverages) / 3))

                    # and finally update the data table
                    new_entries[next_available_entry_id]['coverage'] = coverage

                    # DEALING WITH IT ##################################################################
                    # here we need to put all the codons or amino acids into the data table for this sample
                    for item in set(constants.codons if self.engine == 'CDN' else constants.amino_acids):
                        new_entries[next_available_entry_id][item] = 0

                    # and finally update the frequency of the reference codon or amino acid with the coverage
                    # (WHICH IS VERY BAD WE HAVE NO CLUE WHAT IS THE ACTUAL COVERAGE OF TRIPLICATE LINKMERS):
                    new_entries[next_available_entry_id][reference_item] = coverage

                    next_available_entry_id += 1

        # convert to pandas DataFrame (its faster to build and convert a dictionary than to build
        # DataFrame row by row).
        new_entries = pd.DataFrame(new_entries).T

        # before concatenating the new entries, store the self.data column order. Also, check that
        # no columns exist in new_entries but not in self.data. This is unacceptable, and could have
        # happened if code for new_entries was changed or if the workflow in process() is
        # significantly reworked.
        column_order = self.data.columns.tolist()
        if len([x for x in new_entries.columns.tolist() if x not in self.data.columns.tolist()]):
            raise ValueError("Columns found in new_entries exist that aren't in self.data.")

        # concatenate new columns to self.data
        entries_before = len(self.data.index)
        self.data = pd.concat([self.data, new_entries], sort=True)
        new_entries.set_index("entry_id", drop=False, inplace=True)
        self.data = self.data[column_order]
        entries_after = len(self.data.index)

        self.progress.end()

        # fill in additional fields for new entries. compute_additional_fields takes a list of
        # entry_ids to consider for self.data, which here is provided from new_entries (what I'm
        # saying is new_entries is not passed, only the entry_id's in new_entries
        self.compute_additional_fields(list(new_entries["entry_id"]))

        self.report_change_in_entry_number(entries_before, entries_after, reason="quince mode")


class AminoAcidsEngine(dbops.ContigsSuperclass, VariabilitySuper, QuinceModeWrapperForFancyEngines):
    def __init__(self, args={}, p=progress, r=run):
        self.run = r
        self.progress = p

        self.engine = 'AA'

        # Init Meta
        VariabilitySuper.__init__(self, args=args, r=self.run, p=self.progress)

        # Init the quince mode recoverer
        QuinceModeWrapperForFancyEngines.__init__(self)


    def convert_item_coverages(self):
        for amino_acid in sorted(constants.amino_acids):
            codons = constants.AA_to_codons[amino_acid]
            self.data[amino_acid] = self.data.loc[:, self.data.columns.isin(codons)].sum(axis = 1)
            self.data = self.data.drop(codons, axis=1)


    def convert_reference_info(self):
        self.data['reference'] = self.data['reference'].map(constants.codon_to_AA)
        self.data['departure_from_reference'] = self.data.apply(lambda entry: 1.0 - entry[entry["reference"]] / entry["coverage"], axis=1)


class CodonsEngine(dbops.ContigsSuperclass, VariabilitySuper, QuinceModeWrapperForFancyEngines):
    def __init__(self, args={}, p=progress, r=run):
        self.run = r
        self.progress = p

        self.engine = 'CDN'
        A = lambda x, t: t(args.__dict__[x]) if x in args.__dict__ else None
        null = lambda x: x
        self.skip_synonymity = A('skip_synonymity', null)

        # Init Meta
        VariabilitySuper.__init__(self, args=args, r=self.run, p=self.progress)

        # Init the quince mode recoverer
        QuinceModeWrapperForFancyEngines.__init__(self)

        # add codon specific functions to self.process
        F = lambda f, **kwargs: (f, kwargs)
        self.process_functions.append(F(self.compute_synonymity))


    def compute_synonymity(self):
        """This method is currently prohibitively slow for large datasets."""
        if self.skip_synonymity:
            return

        coding_codons = constants.coding_codons

        number_of_pairs = len(coding_codons)*(len(coding_codons)+1)//2
        array = np.zeros((self.data.shape[0], number_of_pairs))

        array_index = 0
        s_or_ns = []
        for i in coding_codons:
            for j in coding_codons:
                if j > i:
                    break
                array[:, array_index] = self.data.loc[:, i] * self.data.loc[:, j]
                array_index += 1
                s_or_ns.append(constants.is_synonymous[i][j])

        # normalize
        array = array / np.sum(array, axis=1)[:,np.newaxis]

        # each row sums to 1. Synonymity is the sum of those that are synonymous
        synonymity = np.sum(array[:, s_or_ns], axis=1)
        self.data["synonymity"] = synonymity


class ConsensusSequences(NucleotidesEngine, AminoAcidsEngine):
    def __init__(self, args={}, p=progress, r=run):
        self.args = args
        self.run = r
        self.progress = p

        A = lambda x, t: t(args.__dict__[x]) if x in args.__dict__ else None
        null = lambda x: x
        self.engine = A('engine', null)
        self.contigs_mode = A('contigs_mode', null)
        self.genes_of_interest_path = A('genes_of_interest', null)
        self.gene_caller_ids = A('gene_caller_ids', null)
        self.genes_of_interest = A('genes_of_interest', null)
        self.compress_samples = A('compress_samples', null)
        self.tab_delimited_output = A('tab_delimited', null)

        self.sequence_name_key = 'gene_caller_id' if not self.contigs_mode else 'contig_name'

        if not self.engine:
            raise ConfigError("You somehow managed to call the ConsensusSequences class with an args object that does not\
                               contain an engine variable. Not appropriate.")

        if self.engine != 'NT':
            raise ConfigError("Currently the only available variability engine for this is 'NT'. You provided %s" % self.engine)

        if self.compress_samples:
            args.quince_mode = True
            self.run.warning("You supplied --compress-samples, so coverage at each variant position for all sample needs to be\
                              calculated. This will take significantly longer.")
        else:
            args.min_departure_from_reference = 0.5 # if < 0.5, consensus is guaranteed to be reference
                                                    # shortcut only used when not compressing samples

        if self.engine == 'NT':
            NucleotidesEngine.__init__(self, args=args, r=self.run, p=self.progress)
        elif self.engine == 'AA':
            AminoAcidsEngine.__init__(self, args=args, r=self.run, p=self.progress)

        if self.contigs_mode:
            if self.gene_caller_ids or self.genes_of_interest_path:
                raise ConfigError("You can't use --contigs-mode with --gene-caller-ids or --genes-of-interest")
            self.splits_of_interest = set(list(self.splits_basic_info.keys()))

        self.sequence_variants_in_samples_dict = {}


    def populate_seqeunce_variants_in_samples_dict(self):
        """Populates the main dictionary that keeps track of variants for each sample."""
        if self.compress_samples:
            # self data needs to be collapsed
            num_samples = self.data['sample_id'].nunique()
            coverage_columns = self.items + ['coverage']
            not_coverage_columns = ['reference',
                                    'codon_order_in_gene',
                                    'base_pos_in_codon',
                                    'pos_in_contig',
                                    'contig_name',
                                    'sample_id',
                                    'corresponding_gene_call']

            array = self.data[coverage_columns].values
            collapsed_coverage_counts = array.reshape((array.shape[0]//num_samples, num_samples, array.shape[1])).sum(axis=1)
            data_append = self.data[not_coverage_columns].drop_duplicates(subset='pos_in_contig').reset_index(drop=True)
            self.data = pd.DataFrame(collapsed_coverage_counts, columns=coverage_columns).reset_index(drop=True)
            self.data[data_append.columns] = data_append
            self.data['sample_id'] = 'merged'
            self.data['consensus'] = self.data[self.items].idxmax(axis=1)

        # no data no play.
        if not len(self.data):
            raise ConfigError("ConsensusSequences class is upset because it doesn't have any data. There can be two reasons\
                               to this. One, anvi'o variability engines reported nothing (in which case you should have gotten\
                               an error much earler). Two, you are a programmer and failed to call the 'process()' on your\
                               instance from this class. Do you see how the second option is much more likely? :/")

        # learn about the sequences, either contigs or genes
        if self.contigs_mode:
            self.init_contig_sequences()
            sequences = {name: sequence['sequence'].lower() for name, sequence in self.contig_sequences.items()}
        else:
            sequences = {}
            for gene_callers_id in self.genes_of_interest:
                _, d = self.get_sequences_for_gene_callers_ids([gene_callers_id])
                sequences[gene_callers_id] = d[gene_callers_id]['sequence'].lower()

        # here we populate a dictionary with all the right items but witout any real data.
        sample_names = set(self.data['sample_id'])
        for sample_name in sample_names:
            self.sequence_variants_in_samples_dict[sample_name] = {}
            for sequence_name in sequences:
                self.sequence_variants_in_samples_dict[sample_name][sequence_name] = {
                    'sequence_as_list': list(sequences[sequence_name]),
                    'num_changes': 0,
                    self.sequence_name_key: None,
                    'in_pos_0': 0,
                    'in_pos_1': 0,
                    'in_pos_2': 0,
                    'in_pos_3': 0
                }

        # here we will go through every single variant in our data, and correct replace
        # some items in sequences for each sample based on variability infomration.
        self.progress.new('Populating sequence variants in samples data')
        self.progress.update('processing %d variants ...' % len(self.data))
        for idx, entry in self.data.iterrows():
            sample_name = entry['sample_id']
            gene_callers_id = entry['corresponding_gene_call']
            sequence_name = entry['contig_name'] if self.contigs_mode else gene_callers_id

            # the dict item we will be playing with
            d = self.sequence_variants_in_samples_dict[sample_name][sequence_name]

            reference = entry['reference']
            consensus = entry['consensus']

            if reference != consensus:
                codon_order = entry['codon_order_in_gene']
                base_pos_in_codon = entry['base_pos_in_codon']
                if self.contigs_mode:
                    nt_position_to_update = entry['pos_in_contig']
                else:
                    nt_position_to_update = ((codon_order * 3) + base_pos_in_codon) - 1

                # update the entry.
                d[self.sequence_name_key] = sequence_name
                d['sequence_as_list'][nt_position_to_update] = consensus
                d['num_changes'] += 1
                d['in_pos_%d' % base_pos_in_codon] += 1

        self.progress.end()


    def get_formatted_consensus_sequence_entry(self, key, sample_name, sequence_name):
        """Gets a sample name and gene callers id, returns a well-formatted dict for sequence
           entry using the `self.sequence_variants_in_samples_dict`.

           `key` must be unique string identifier."""

        F = self.sequence_variants_in_samples_dict[sample_name][sequence_name]
        return {'key': key,
                'sample_name':  sample_name,
                self.sequence_name_key: sequence_name,
                'num_changes': F['num_changes'],
                'in_pos_0': F['in_pos_0'],
                'in_pos_1': F['in_pos_1'],
                'in_pos_2': F['in_pos_2'],
                'in_pos_3': F['in_pos_3'],
                'sequence': ''.join(F['sequence_as_list'])}


    def report(self):
        if not self.sequence_variants_in_samples_dict:
            self.populate_seqeunce_variants_in_samples_dict()

        self.progress.new('Generating the report')
        self.progress.update('...')

        output_d = {}
        counter = 1
        for sample_name in self.sequence_variants_in_samples_dict:
            for sequence_name in self.sequence_variants_in_samples_dict[sample_name]:
                d = self.get_formatted_consensus_sequence_entry('e%.7d' % (counter), sample_name, sequence_name)
                counter += 1

                if self.tab_delimited_output:
                    output_d[d['key']] = d
                else:
                    key = '|'.join([d['key'],
                                    'sample_name:%s' % d['sample_name'],
                                    '{}:{}'.format(self.sequence_name_key, d[self.sequence_name_key]),
                                    'num_changes:%s' % d['num_changes'],
                                    'in_pos_0:%d' % d['in_pos_0'],
                                    'in_pos_1:%d' % d['in_pos_1'],
                                    'in_pos_2:%d' % d['in_pos_2'],
                                    'in_pos_3:%d' % d['in_pos_3']])
                    output_d[key] = d['sequence']

        if self.tab_delimited_output:
            utils.store_dict_as_TAB_delimited_file(output_d, self.output_file_path, headers=['entry_id', 'sample_name', self.sequence_name_key, 'num_changes', 'in_pos_0', 'in_pos_1', 'in_pos_2', 'in_pos_3', 'sequence'])
        else:
            utils.store_dict_as_FASTA_file(output_d, self.output_file_path)

        self.progress.end()

        self.run.info('Num genes reported', pp(len(self.genes_of_interest)))
        self.run.info('Num sequences reported', pp(len(self.sequence_variants_in_samples_dict)))
        self.run.info('Output File', self.output_file_path)


class VariabilityNetwork:
    def __init__(self, args={}, p=progress, r=run):
        self.args = args

        self.run = r
        self.progress = p

        self.samples = None
        self.samples_information_dict = None
        self.data = None

        A = lambda x, t: t(args.__dict__[x]) if x in args.__dict__ else None
        null = lambda x: x
        self.input_file_path = A('input_file', null)
        self.samples_information_path = A('samples_information', null)
        self.max_num_unique_positions = A('max_num_unique_positions', int)
        self.output_file_path = A('output_file', null)

        filesnpaths.is_output_file_writable(self.output_file_path)

        if self.samples_information_path:
            filesnpaths.is_file_tab_delimited(self.samples_information_path)
            self.samples_information_dict = utils.get_TAB_delimited_file_as_dictionary(self.samples_information_path)
            num_attributes = len(list(self.samples_information_dict.values())[0])

            self.run.info('samples_information', '%d attributes read for %d samples' % (num_attributes, len(self.samples_information_dict)))

        if self.input_file_path:
            filesnpaths.is_file_tab_delimited(self.input_file_path)
            self.progress.new('Reading the input file')
            self.progress.update('...')
            self.data = utils.get_TAB_delimited_file_as_dictionary(self.input_file_path)
            self.progress.end()

            self.run.info('input_file', '%d entries read' % len(self.data))


    def generate(self):
        if not self.data:
            raise ConfigError("There is nothing to report. Either the input file you provided was empty, or you\
                                haven't filled in the variable positions data into the class.")

        if self.max_num_unique_positions < 0:
            raise ConfigError("Max number of unique positions cannot be less than 0.. Obviously :/")

        self.samples = sorted(list(set([e['sample_id'] for e in list(self.data.values())])))
        self.run.info('samples', '%d found: %s.' % (len(self.samples), ', '.join(self.samples)))

        if self.samples_information_dict:
            samples_missing_in_information_dict = [s for s in self.samples if s not in self.samples_information_dict]
            if len(samples_missing_in_information_dict):
                raise ConfigError("The sample names you provided in the samples information data is not a subset of\
                                    sample names found in the variable positions data :/ Essentially, every sample name\
                                    appears in the variability data must be present in the samples information data,\
                                    however, you are missing these ones from your samples information: %s."\
                                                % (', '.join(samples_missing_in_information_dict)))

        self.unique_variable_nt_positions = set([e['unique_pos_identifier'] for e in list(self.data.values())])
        self.run.info('unique_variable_nt_positions', '%d found.' % (len(self.unique_variable_nt_positions)))

        if self.max_num_unique_positions and len(self.unique_variable_nt_positions) > self.max_num_unique_positions:
            self.unique_variable_nt_positions = set(random.sample(self.unique_variable_nt_positions, self.max_num_unique_positions))
            self.run.info('unique_variable_nt_positions', 'Unique positions are subsampled to %d' % self.max_num_unique_positions, mc='red')

        self.progress.new('Samples dict')
        self.progress.update('Creating an empty one ...')
        samples_dict = {}
        for sample_name in self.samples:
            samples_dict[sample_name] = {}
            for unique_variable_position in self.unique_variable_nt_positions:
                samples_dict[sample_name][unique_variable_position] = 0

        self.progress.update('Updating the dictionary with data')
        for entry in list(self.data.values()):
            sample_id = entry['sample_id']
            pos = entry['unique_pos_identifier']
            frequency = entry['departure_from_reference']

            samples_dict[sample_id][pos] = float(frequency)

        self.progress.update('Generating the network file')
        utils.gen_gexf_network_file(sorted(list(self.unique_variable_nt_positions)), samples_dict, self.output_file_path, sample_mapping_dict=self.samples_information_dict)
        self.progress.end()

        self.run.info('network_description', self.output_file_path)


class VariabilityData(NucleotidesEngine, CodonsEngine, AminoAcidsEngine):
    def __init__(self, args={}, p=progress, r=run, dont_process=False):
        self.progress = p
        self.run = r

        self.args = args
        A = lambda x, t: t(args.__dict__[x]) if x in args.__dict__ else None
        self.columns_to_load = A('columns_to_load', list)
        self.variability_table_path = A('variability_profile', str)
        self.engine = A('engine', str)

        if not self.variability_table_path:
            raise ConfigError("VariabilityData :: You must declare a variability table filepath.")

        # determine the engine type of the variability table
        inferred_engine = utils.get_variability_table_engine_type(self.variability_table_path)
        if self.engine and self.engine != inferred_engine:
            raise ConfigError("The engine you requested is {}, but the engine inferred from {} is {}.".\
                               format(self.engine, self.variability_table_path, inferred_engine))
        self.engine = inferred_engine

        if self.engine == 'NT':
            self.items = constants.nucleotides
            self.competing_items = 'competing_nts'
        elif self.engine == 'CDN':
            self.items = constants.codons
            self.competing_items = 'competing_codons'
        elif self.engine == 'AA':
            self.items = constants.amino_acids
            self.competing_items = 'competing_aas'
        else:
            pass

        if not dont_process:
            self.process_external_table()


    def load_data(self):
        """load the variability data (output of anvi-gen-variabliity-profile)"""
        self.data = pd.read_csv(self.variability_table_path, sep="\t", usecols=self.columns_to_load)


    def process_external_table(self):
        # load the data
        self.load_data()

        # init the appropriate engine
        self.args.data = self.data
        self.args.engine = self.engine
        self.args.skip_sanity_check = True
        variability_engines[self.engine].__init__(self, self.args, p=self.progress, r=self.run)

        # load residue info data
        if self.append_structure_residue_info:
            self.load_structure_data()

        self.init_commons()


variability_engines = {'NT': NucleotidesEngine, 'CDN': CodonsEngine, 'AA': AminoAcidsEngine}

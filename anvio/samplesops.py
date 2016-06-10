# -*- coding: utf-8
# pylint: disable=line-too-long
"""Module to make sense of samples information and samples order input"""


import anvio.utils as utils
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths

from anvio.errors import SamplesError


__author__ = "A. Murat Eren"
__copyright__ = "Copyright 2015, The anvio Project"
__credits__ = []
__license__ = "GPL 3.0"
__maintainer__ = "A. Murat Eren"
__email__ = "a.murat.eren@gmail.com"


progress = terminal.Progress()
run = terminal.Run()


class SamplesInformation:
    def __init__(self, run=run, progress=progress, quiet=False):
        self.samples_information_dict = {}
        self.samples_order_dict = {}

        self.aliases_to_attributes_dict = {}

        self.available_orders = set([])

        self.sample_names_in_samples_order_file = None
        self.sample_names_in_samples_information_file = None
        self.samples_information_default_layer_order = None

        self.sample_names = None

        self.run = run
        self.prgress = progress
        self.quiet = quiet


    def process_samples_information_file(self, samples_information_path):
        if not samples_information_path:
            return

        self.sample_names_in_samples_information_file = filesnpaths.is_proper_samples_information_file(samples_information_path)

        self.samples_information_dict, self.aliases_to_attributes_dict = self.convert_samples_information_dict(utils.get_TAB_delimited_file_as_dictionary(samples_information_path))
        self.samples_information_default_layer_order = open(samples_information_path, 'rU').readline().strip().split('\t')[1:]

        self.run.info('Samples information', 'Loaded for %d samples' % len(self.samples_information_dict), quiet=self.quiet)


    def recover_samples_information_dict(self, samples_information_dict_from_db, aliases_to_attributes_dict):
        samples_information_dict_with_attributes = {}

        for sample_name in samples_information_dict_from_db:
            samples_information_dict_with_attributes[sample_name] = {}

        for alias in aliases_to_attributes_dict:
            attribute = aliases_to_attributes_dict[alias]['attribute']
            for sample_name in samples_information_dict_with_attributes:
                samples_information_dict_with_attributes[sample_name][attribute] = samples_information_dict_from_db[sample_name][alias]

        return samples_information_dict_with_attributes


    def convert_samples_information_dict(self, samples_information_dict_from_file):
        """Create aliases for each user-declared sample attribute.

           It is important to note that each attribute becomes a database field, and
           databases have limitations on those names that leave no room for creativity.
           For instance, anvi'o uses "X;Y;Z" notation to define a field for bar charts.
           Yet this is not a valid database column.
        """
        samples_information_dict_with_aliases = {}
        aliases_to_attributes_dict = {}

        for sample_name in samples_information_dict_from_file:
            samples_information_dict_with_aliases[sample_name] = {}

        alias_index = 1
        for attribute in samples_information_dict_from_file.values()[0]:
            alias = 'attr_%05d' % alias_index
            aliases_to_attributes_dict[alias] = attribute
            alias_index += 1

            for sample_name in samples_information_dict_from_file:
                samples_information_dict_with_aliases[sample_name][alias] = samples_information_dict_from_file[sample_name][attribute]

        return samples_information_dict_with_aliases, aliases_to_attributes_dict


    def process_samples_order_file(self, samples_order_path):
        if not samples_order_path:
            return

        self.sample_names_in_samples_order_file = filesnpaths.is_proper_samples_order_file(samples_order_path)

        self.samples_order_dict = utils.get_TAB_delimited_file_as_dictionary(samples_order_path)

        self.available_orders = set(self.samples_order_dict.keys())

        if not self.samples_information_default_layer_order:
            pass

        self.run.info('Samples order', 'Loaded for %d attributes' % len(self.samples_order_dict), quiet=self.quiet)


    def update_samples_order_dict(self):
        """Some attributes in the samples information dict may also be used as orders"""

        def F(v):
            if not v:
                return 0.0
            try:
                return float(v)
            except:
                return v

        for sample_attribute_tuples in [[(F(self.samples_information_dict[sample][attribute]), sample, attribute) \
                                            for sample in self.samples_information_dict] \
                                            for attribute in self.aliases_to_attributes_dict]:
            # skip bar charts:
            if ';' in str(sample_attribute_tuples[0][0]):
                continue

            attribute = self.aliases_to_attributes_dict[sample_attribute_tuples[0][2]]
            if attribute not in self.samples_order_dict:
                self.samples_order_dict['>> ' + attribute] = {'newick': '', 'basic': ','.join([t[1] for t in sorted(sample_attribute_tuples)])}
                self.samples_order_dict['>> ' + attribute + ' (reverse)'] = {'newick': '', 'basic': ','.join([t[1] for t in sorted(sample_attribute_tuples, reverse=True)])}


    def populate_from_input_files(self, samples_information_path=None, samples_order_path=None):
        self.process_samples_information_file(samples_information_path)
        self.process_samples_order_file(samples_order_path)
        self.update_samples_order_dict()

        self.sanity_check()

        self.sample_names = self.sample_names_in_samples_information_file or self.sample_names_in_samples_order_file


    def sanity_check(self):
        if self.sample_names_in_samples_information_file and self.sample_names_in_samples_order_file:
            if sorted(self.sample_names_in_samples_information_file) != sorted(self.sample_names_in_samples_order_file):
                raise SamplesError, 'OK. Samples described in the information file and order file are not identical :/\
                                     Here are the %d sample names in the information file: "%s", versus the %d sample\
                                     names in the orders file: "%s". And here is the difference: "%s".'\
                                                            % (len(self.sample_names_in_samples_information_file),
                                                               self.sample_names_in_samples_information_file,
                                                               len(self.sample_names_in_samples_order_file),
                                                               self.sample_names_in_samples_order_file,
                                                               list(set(self.sample_names_in_samples_information_file) - set(self.sample_names_in_samples_order_file)))

        if not self.samples_information_default_layer_order:
            # we still don't have a default order. we will try to recover from that here
            # by looking into what we have in the samples order informaiton
            if not len(self.samples_order_dict):
                raise SamplesError, "Something is missing. Anvi'o is having hard time coming up with a default samples\
                                    order for the samples database."

            a_basic_order = [o['basic'] for o in self.samples_order_dict.values()][0]
            a_tree_order = utils.get_names_order_from_newick_tree([o['newick'] for o in self.samples_order_dict.values()][0])

            self.samples_information_default_layer_order = a_basic_order or a_tree_order

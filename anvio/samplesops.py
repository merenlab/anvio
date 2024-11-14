# -*- coding: utf-8
# pylint: disable=line-too-long
"""Module to make sense of samples information and samples order input"""

import os

import anvio.utils as utils
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths

from anvio.errors import SamplesError


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
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

        self.sample_names_in_samples_information_file = (
            filesnpaths.is_proper_samples_information_file(samples_information_path)
        )

        self.samples_information_dict, self.aliases_to_attributes_dict = (
            self.convert_samples_information_dict(
                utils.get_TAB_delimited_file_as_dictionary(samples_information_path)
            )
        )
        self.samples_information_default_layer_order = (
            open(samples_information_path, "r").readline().strip().split("\t")[1:]
        )

        self.run.info(
            "Samples information",
            "Loaded for %d samples" % len(self.samples_information_dict),
            quiet=self.quiet,
        )

    def recover_samples_information_dict(
        self, samples_information_dict_from_db, aliases_to_attributes_dict
    ):
        samples_information_dict_with_attributes = {}

        for sample_name in samples_information_dict_from_db:
            samples_information_dict_with_attributes[sample_name] = {}

        for alias in aliases_to_attributes_dict:
            attribute = aliases_to_attributes_dict[alias]["attribute"]
            for sample_name in samples_information_dict_with_attributes:
                samples_information_dict_with_attributes[sample_name][attribute] = (
                    samples_information_dict_from_db[sample_name][alias]
                )

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
        for attribute in list(samples_information_dict_from_file.values())[0]:
            alias = "attr_%05d" % alias_index
            aliases_to_attributes_dict[alias] = attribute
            alias_index += 1

            for sample_name in samples_information_dict_from_file:
                samples_information_dict_with_aliases[sample_name][alias] = (
                    samples_information_dict_from_file[sample_name][attribute]
                )

        return samples_information_dict_with_aliases, aliases_to_attributes_dict

    def process_samples_order_file(self, samples_order_path):
        if not samples_order_path:
            return

        self.sample_names_in_samples_order_file = (
            filesnpaths.is_proper_samples_order_file(samples_order_path)
        )

        self.samples_order_dict = utils.get_TAB_delimited_file_as_dictionary(
            samples_order_path
        )

        self.available_orders = set(self.samples_order_dict.keys())

        if not self.samples_information_default_layer_order:
            pass

        self.run.info(
            "Samples order",
            "Loaded for %d attributes" % len(self.samples_order_dict),
            quiet=self.quiet,
        )

    def process_single_order_data(self, single_order_path, single_order_name):
        """Just inject a single order into the `self.samples_order_dict`"""

        if not single_order_path:
            return

        if not single_order_name:
            raise SamplesError(
                "You provided a file for a single order, but not a name for it. This is a no no :/"
            )

        filesnpaths.is_file_plain_text(single_order_path)

        single_order_file_content = [
            l.strip("\n") for l in open(single_order_path, "r").readlines()
        ]

        if len(single_order_file_content) != 1:
            raise SamplesError(
                "The single order file should contain a single line of information. It can't have nothing,\
                                it can't have too much. Just a single newick tree, or a comma-separated list of sample\
                                names."
            )

        _order = single_order_file_content.pop()

        # if you are reading this line, please brace yourself to possibly one of the silliest
        # bunch of lines in the anvi'o codebase. the reason we are doing this this way is quite
        # a long story, and deserves a FIXME, but in order to utilize the excellent function
        # in the filesnpaths module to check the contents of the samples order dict rigirously,
        # we need to have this information in a file. a better way could have been implementing
        # a filesnpaths.is_proper_samples_order_content function next to the currently available
        # filesnpaths.is_proper_samples_order_file (the latter would call the former with a dict
        # and it would be much more flexible), but we can't import utils form within filesnpaths.
        # without utils we don't have a get_TAB_delimited_file_as_dictionary function, and we are
        # definitely not going to implement it in two places :( recovering from a poor design by
        # doing something even poorer? couldn't have we fixed this once and for all instead of
        # writing this paragraph? well. just remember that you are thinking about a rethorical
        # question in a comment section. so sometimes we do things that are not quite productive.
        temp_samples_order_file_path = filesnpaths.get_temp_file_path()
        temp_samples_order_file = open(temp_samples_order_file_path, "w")
        temp_samples_order_file.write(
            "\t".join(["attributes", "basic", "newick"]) + "\n"
        )

        if filesnpaths.is_proper_newick(_order, dont_raise=True):
            temp_samples_order_file.write(
                "\t".join([single_order_name, "", _order]) + "\n"
            )
            self.samples_order_dict[single_order_name] = {
                "newick": _order,
                "basic": None,
            }
        else:
            temp_samples_order_file.write(
                "\t".join([single_order_name, _order, ""]) + "\n"
            )
            self.samples_order_dict[single_order_name] = {
                "basic": _order,
                "newick": None,
            }

        temp_samples_order_file.close()

        sample_names_in_samples_order_file = filesnpaths.is_proper_samples_order_file(
            temp_samples_order_file_path
        )
        os.remove(temp_samples_order_file_path)

        if not self.sample_names_in_samples_information_file:
            self.sample_names_in_samples_order_file = sample_names_in_samples_order_file

        self.available_orders.add(single_order_name)

        self.run.info(
            "Samples order",
            "A single order for '%s' is also loaded" % single_order_name,
            quiet=self.quiet,
        )

    def update_samples_order_dict(self):
        """Some attributes in the samples information dict may also be used as orders"""

        def F(v):
            if isinstance(v, type(None)):
                return ""

            if not v:
                return 0.0

            try:
                return float(v)
            except:
                return v

        for sample_attribute_tuples in [
            [
                (F(self.samples_information_dict[sample][attribute]), sample, attribute)
                for sample in self.samples_information_dict
            ]
            for attribute in self.aliases_to_attributes_dict
        ]:
            # skip bar charts:
            if ";" in str(sample_attribute_tuples[0][0]):
                continue

            attribute = self.aliases_to_attributes_dict[sample_attribute_tuples[0][2]]
            if attribute not in self.samples_order_dict:
                try:
                    self.samples_order_dict[">> " + attribute] = {
                        "newick": "",
                        "basic": ",".join(
                            [t[1] for t in sorted(sample_attribute_tuples)]
                        ),
                    }
                    self.samples_order_dict[">> " + attribute + " (reverse)"] = {
                        "newick": "",
                        "basic": ",".join(
                            [
                                t[1]
                                for t in sorted(sample_attribute_tuples, reverse=True)
                            ]
                        ),
                    }
                except TypeError:
                    raise SamplesError(
                        "OK. Anvi'o has good and bad news. The bad news is that your samples information "
                        "is kaput, because one of the columns in it has mixed data types (not everything has the "
                        "same type). The good news is that we know what column is that: it is the column '%s'. "
                        "Please take a look." % attribute
                    )

    def populate_from_input_files(
        self,
        samples_information_path=None,
        samples_order_path=None,
        single_order_path=None,
        single_order_name=None,
    ):
        if (
            not samples_information_path
            and not samples_order_path
            and not single_order_path
        ):
            raise SamplesError(
                "At least one of the input files must be declared to create or to update an "
                "anvi'o samples information database :/ But maybe not. Maybe anvi'o should be "
                "able to create an empty samples information database, too. Do you need this? "
                "Write to us!"
            )

        self.process_samples_information_file(samples_information_path)
        self.process_samples_order_file(samples_order_path)
        self.process_single_order_data(single_order_path, single_order_name)
        self.update_samples_order_dict()

        self.sanity_check()

        self.sample_names = (
            self.sample_names_in_samples_information_file
            or self.sample_names_in_samples_order_file
        )

    def sanity_check(self):
        if (
            self.sample_names_in_samples_information_file
            and self.sample_names_in_samples_order_file
        ):
            if sorted(self.sample_names_in_samples_information_file) != sorted(
                self.sample_names_in_samples_order_file
            ):
                raise SamplesError(
                    "OK. Samples described in the information file and order file are not identical :/ "
                    'Here are the %d sample names in the information file: "%s", versus the %d sample '
                    'names in the orders file: "%s". And here is the difference: "%s".'
                    % (
                        len(self.sample_names_in_samples_information_file),
                        self.sample_names_in_samples_information_file,
                        len(self.sample_names_in_samples_order_file),
                        self.sample_names_in_samples_order_file,
                        list(
                            set(self.sample_names_in_samples_information_file)
                            - set(self.sample_names_in_samples_order_file)
                        ),
                    )
                )

        if not self.samples_information_default_layer_order:
            # we still don't have a default order. we will try to recover from that here
            # by looking into what we have in the samples order informaiton
            if not len(self.samples_order_dict):
                raise SamplesError(
                    "Something is missing. Anvi'o is having hard time coming up with a default samples "
                    "order for the samples database."
                )

            a_basic_order = [
                o["basic"].split(",") if o["basic"] else None
                for o in list(self.samples_order_dict.values())
            ][0]
            a_tree_order = utils.get_names_order_from_newick_tree(
                [
                    o["newick"] if o["newick"] else None
                    for o in list(self.samples_order_dict.values())
                ][0]
            )

            self.samples_information_default_layer_order = a_basic_order or a_tree_order

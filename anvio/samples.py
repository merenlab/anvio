# -*- coding: utf-8
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


class Samples:
    def __init__(self, run = run, progress = progress):
        self.samples_information_dict = {}
        self.samples_order_dict = {}

        self.run = run
        self.prgress = progress


    def process_samples_information_file(self, samples_information_path):
        if not samples_information_path:
            return

        filesnpaths.is_proper_samples_information_file(samples_information_path)

        self.samples_information_dict = utils.get_TAB_delimited_file_as_dictionary(samples_information_path)

        self.run.info('Samples information', 'Loaded for %d samples' % len(self.samples_information_dict))


    def process_samples_order_file(self, samples_order_path):
        if not samples_order_path:
            return

        filesnpaths.is_proper_samples_order_file(samples_order_path)

        self.samples_order_dict = utils.get_TAB_delimited_file_as_dictionary(samples_order_path)

        self.run.info('Samples order', 'Loaded for %d attributes' % len(self.samples_order_dict))


    def populate_from_input_files(self, samples_information_path = None, samples_order_path = None):
        self.process_samples_information_file(samples_information_path)
        self.process_samples_order_file(samples_order_path)

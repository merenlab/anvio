# -*- coding: utf-8
# pylint: disable=line-too-long

import numpy

import anvio
import anvio.tables as t
import anvio.utils as utils
import anvio.terminal as terminal


__author__ = "Developers of anvi'o (see AUTHORS.txt)"
__copyright__ = "Copyleft 2015-2018, the Meren Lab (http://merenlab.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "A. Murat Eren"
__email__ = "a.murat.eren@gmail.com"
__status__ = "Development"


run = terminal.Run()
progress = terminal.Progress()
pp = terminal.pretty_print


class TableForNtPositions(object):
    def __init__(self, run=run, progress=progress):
        self.run = run
        self.progress = progress
        self.numpy_data_type = 'uint8'
        self.db_entries = []


    def append(self, contig_name, position_info_list):
        position_info_blob = utils.convert_numpy_array_to_binary_blob(numpy.array(position_info_list, dtype=self.numpy_data_type))
        self.db_entries.append((contig_name, position_info_blob, ))


    def store(self, db):
        db.insert_many(t.nt_position_info_table_name, entries=self.db_entries)
        self.db_entries = []

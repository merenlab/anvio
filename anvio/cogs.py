# -*- coding: utf-8
# pylint: disable=line-too-long
"""
    Making sense of COGs.
"""

import os

import anvio
import anvio.utils as utils
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths


__author__ = "A. Murat Eren"
__copyright__ = "Copyright 2016, The anvio Project"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "A. Murat Eren"
__email__ = "a.murat.eren@gmail.com"


run = terminal.Run()
progress = terminal.Progress()
pp = terminal.pretty_print


COGs_PATH = os.path.join(os.path.dirname(anvio.__file__), 'data/misc/COGs/COGs.txt')
CATEGORIES_PATH = os.path.join(os.path.dirname(anvio.__file__), 'data/misc/COGs/CATEGORIES.txt')


class COGsData:
    """A class to make sense of COG ids and categories"""
    def __init__(self, COGs_path=COGs_PATH, categories_path=CATEGORIES_PATH):
        filesnpaths.is_file_tab_delimited(COGs_path)
        filesnpaths.is_file_tab_delimited(categories_path)

        self.cogs = utils.get_TAB_delimited_file_as_dictionary(COGs_path, no_header=True, column_names=['COG', 'categories', 'annotation'])
        self.categories = utils.get_TAB_delimited_file_as_dictionary(categories_path, no_header=True, column_names=['category', 'description'])

        for cog in self.cogs:
            self.cogs[cog]['categories'] = [c.strip() for c in self.cogs[cog]['categories'].split(',')]

        for cat in self.categories:
            self.categories[cat] = self.categories[cat]['description']

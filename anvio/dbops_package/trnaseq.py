# -*- coding: utf-8 -*-
# pylint: disable=line-too-long

import os

import anvio
import anvio.dbops_package as dbops_package
import anvio.constants as constants
import anvio.tables as t
import anvio.terminal as terminal
import anvio.utils as utils

__author__ = "Developers of anvi'o (see AUTHORS.txt)"
__copyright__ = "Copyleft 2015-2020, the Meren Lab (http://merenlab.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__maintainer__ = "Samuel Miller"
__email__ = "samuelmiller10@gmail.com"
__status__ = "Development"


class TRNASeqDatabase(dbops_package.Database):
    """ Used to initialize an empty tRNA-seq database or access an existing tRNA-seq database"""

    def __init__(self, db_path, run=terminal.Run(), progress=terminal.Progress(), quiet=True):
        args = anvio.EmptyArgs()
        if not os.path.exists(db_path):
            args.db_type = 'trnaseq'
            args.db_version = anvio.__trnaseq__version__
        args.meta_int_keys = [] # metadata to be stored as an int
        args.meta_float_keys = [] # metadata to be stored as a float
        args.table_info = [
            (t.trnaseq_sequences_table_name, t.trnaseq_sequences_table_structure, t.trnaseq_sequences_table_types),
            (t.trnaseq_info_table_name, t.trnaseq_info_table_structure, t.trnaseq_info_table_types),
            (t.trnaseq_features_table_name, t.trnaseq_features_table_structure, t.trnaseq_features_table_types),
            (t.trnaseq_unconserved_table_name, t.trnaseq_unconserved_table_structure, t.trnaseq_unconserved_table_types),
            (t.trnaseq_unpaired_table_name, t.trnaseq_unpaired_table_structure, t.trnaseq_unpaired_table_types),
            (t.trnaseq_trimmed_table_name, t.trnaseq_trimmed_table_structure, t.trnaseq_trimmed_table_types),
            (t.trnaseq_normalized_table_name, t.trnaseq_normalized_table_structure, t.trnaseq_normalized_table_types),
            (t.trnaseq_modified_table_name, t.trnaseq_modified_table_structure, t.trnaseq_modified_table_types)]

        super().__init__(db_path, args=args, run=run, progress=progress, quiet=quiet)
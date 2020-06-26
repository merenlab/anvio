# -*- coding: utf-8 -*-
# pylint: disable=line-too-long

import os

import anvio
import anvio.dbops as dbops
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


THREEPRIME_VARIANT_LIST = constants.trnaseq.THREEPRIME_VARIANT_LIST


class TRNASeqSuperclass:
    def __init__(self, args, run=terminal.Run(), progress=terminal.Progress()):
        pass


class TRNASeqDatabase(dbops.Database):
    def __init__(self, db_path, run=terminal.Run(), progress=terminal.Progress(), quiet=True):
        args = {}
        args['meta_int_keys'] = ['charging_recorded']
        args['meta_float_keys'] = []
        args['table_info'] = [(t.trnaseq_sequences_table_name, t.trnaseq_sequences_table_structure, t.trnaseq_sequences_table_types),
                              (t.trnaseq_acceptor_table_name, t.trnaseq_acceptor_table_structure, t.trnaseq_acceptor_table_types),
                              (t.trnaseq_subsequence_table_name, t.trnaseq_subsequence_table_structure, t.trnaseq_subsequence_table_types),
                              (t.trnaseq_info_table_name, t.trnaseq_info_table_structure, t.trnaseq_info_table_types),
                              (t.trnaseq_features_table_name, t.trnaseq_features_table_structure, t.trnaseq_features_table_types),
                              (t.trnaseq_unconserved_table_name, t.trnaseq_unconserved_table_structure, t.trnaseq_unconserved_table_types),
                              (t.trnaseq_unpaired_table_name, t.trnaseq_unpaired_table_structure, t.trnaseq_unpaired_table_types),
                              (t.trnaseq_long_sequences_table_name, t.trnaseq_long_sequences_table_structure, t.trnaseq_long_sequences_table_types),
                              (t.trnaseq_long_acceptor_table_name, t.trnaseq_long_acceptor_table_structure, t.trnaseq_long_acceptor_table_types),
                              (t.trnaseq_long_subsequence_table_name, t.trnaseq_long_subsequence_table_structure, t.trnaseq_long_subsequence_table_types),
                              (t.trnaseq_long_info_table_name, t.trnaseq_long_info_table_structure, t.trnaseq_long_info_table_types),
                              (t.trnaseq_long_features_table_name, t.trnaseq_long_features_table_structure, t.trnaseq_long_features_table_types),
                              (t.trnaseq_long_unconserved_table_name, t.trnaseq_long_unconserved_table_structure, t.trnaseq_long_unconserved_table_types),
                              (t.trnaseq_long_unpaired_table_name, t.trnaseq_long_unpaired_table_structure, t.trnaseq_long_unpaired_table_types)]

        super().__init__(db_path, args=args, run=run, progress=progress, quiet=quiet)
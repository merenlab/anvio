# -*- coding: utf-8
# pylint: disable=line-too-long
"""Gather data from SSMs"""

import os
import glob

import anvio.utils as u
import anvio.terminal as terminal
import anvio.constants as constants

from anvio.errors import ConfigError

run = terminal.Run()

engines = ['AA', 'CDN', 'NT']

def get(engine, run=run):
    data = {}

    if engine not in engines:
        raise ConfigError("Anvi'o was about to populate the SSMs, but it does not know about the engine '%s'." % engine)

    dir_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), engine)
    substitution_matrix_paths = [s for s in glob.glob(os.path.join(dir_path, '*')) if s.endswith('.txt')]

    for matrix_path in substitution_matrix_paths:
        matrix_id = os.path.basename(matrix_path).split('.txt')[0]

        matrix_rows = u.get_column_data_from_TAB_delim_file(matrix_path, column_indices=[0])[0][1:]
        matrix_columns = u.get_columns_of_TAB_delim_file(matrix_path, include_first_column=False)

        if sorted(matrix_columns) != sorted(matrix_rows):
            raise ConfigError("Anvi'o found a substitution scoring matrix named '%s'. However, it doesn't look like\
                               a nicely done matrix. Substitution scoring matrices must contain the same row and column\
                               names (i.e., a square matrix that is equal to its transpose). Well. This one does not :/" \
                                                    % (os.path.basename(matrix_path)))

        if engine == 'AA':
            expected_items = set(list(constants.amino_acids))
        elif engine == 'NT':
            expected_items = set(list(constants.nucleotides))
        elif engine == 'CDN':
            expected_items = set(list(constants.codons))

        unexpected_items_in_matrix = [item for item in matrix_columns if item not in expected_items]
        if len(unexpected_items_in_matrix):
            raise ConfigError("It seems you have a poorly done substitution scoring matrix named '%s' in the data directory.\
                               Anvi'o expects an %s substitution matrix to describe one or more of these %d guys: '%s'. But\
                               the matrix %s had stuff anvi'o is not familiar with: '%s'." % \
                                            (matrix_id, engine, len(expected_items), ', '.join(expected_items),
                                             matrix_id, ', '.join(unexpected_items_in_matrix)))

        matrix_data = u.get_TAB_delimited_file_as_dictionary(matrix_path, column_mapping = [str] + [float] * len(expected_items))
        data[matrix_id] = matrix_data

    if len(data):
        run.warning('%d matri%s been loaded: "%s".' % \
                                    (len(data),
                                     'ces have' if len(data) > 1 else 'x has',
                                     ', '.join(list(data.keys()))),
                    header='%s substitution scoring matrices' % engine,
                    lc="green")

    return data

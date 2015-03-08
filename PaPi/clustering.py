#!/usr/bin/env python
# -*- coding: utf-8

# Copyright (C) 2014, A. Murat Eren
#
# This program is free software; you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free
# Software Foundation; either version 2 of the License, or (at your option)
# any later version.
#
# Please read the COPYING file.

import numpy as np
import PaPi.utils as utils
import PaPi.terminal as terminal


def set_num_components_for_each_matrix(config):
    denominator = float(sum([r['ratio'] for r in config.matrices_dict.values()]))
    for matrix in config.matrices:
        m = config.matrices_dict[matrix]
        num_components_for_ratio = int(round(config.num_components * (m['ratio'] / denominator)))
        m['num_components'] = num_components_for_ratio

    return config


def set_null_ratios_for_matrices(config):
    for matrix in config.matrices:
        config.matrices_dict[matrix]['ratio'] = 1
    return config


def order_contigs_simple(config, progress = terminal.Progress(verbose=False), run = terminal.Run(), debug = False):
    if not config.matrices_dict[config.matrices[0]]['ratio']:
        config = set_null_ratios_for_matrices(config)

    for matrix in config.matrices:
        m = config.matrices_dict[matrix]
        if m['normalize']:
            m['scaled_vectors'] = utils.get_normalized_vectors(m['vectors'])
        else:
            m['scaled_vectors'] = np.array(m['vectors'])

    progress.new('Vectors from %d matrices' % config.num_matrices)
    progress.update('Combining ...')
    config.combined_vectors = []
    config.combined_id_to_sample = {}

    for i in range(0, len(config.master_rows)):
        row = config.master_rows[i]
        config.combined_id_to_sample[i] = config.master_rows[i]
        combined_scaled_vectors_for_row = [m['scaled_vectors'][m['sample_to_id'][row]] for m in config.matrices_dict.values()]
        config.combined_vectors.append(np.concatenate(combined_scaled_vectors_for_row))

    if debug:
        for i in range(0, 10):
            print config.combined_vectors[i]

    progress.update('Clustering ...')
    tree = utils.get_clustering_as_tree(config.combined_vectors, progress = progress)
    newick = utils.get_tree_object_in_newick(tree, config.combined_id_to_sample)
    progress.end()

    if config.output_file_path:
        open(config.output_file_path, 'w').write(newick + '\n')

    return newick


def order_contigs_experimental(config, progress = terminal.Progress(verbose=False), run = terminal.Run(), debug = False):
    if not config.multiple_matrices:
        # there is one matrix. could be coverage, could be tnf. we don't care.
        # we do what we gotta do: skip scaling and perform clustering using all
        # dimensions.
        m = config.matrices_dict[config.matrices[0]]

        progress.new('Single matrix (%s)' % m['alias'])
        progress.update('Performing cluster analysis ...')
        tree = utils.get_clustering_as_tree(m['vectors'], progress = progress)
        newick = utils.get_tree_object_in_newick(tree, m['id_to_sample'])
        progress.end()

        if config.output_file_path:
            open(config.output_file_path, 'w').write(newick + '\n')

        return newick

    else:
        # FIXME: this part needs to be parallelized.

        # ok. there is more than one matrix, so there will be a mixture of scaled vectors prior to clustering.
        # we first will determine whether ratios were set in the config file. if ratios were not set the simplest
        # thing to do is to equally distributing num_components across all matrices; so we will set ratios to 1.
        # a heuristic that handles the initial config file before calling this function can determine what ratios
        # would be appropriate considering the number of samples in the experiment and/or other experiment-specific
        # properties
        if not config.matrices_dict[config.matrices[0]]['ratio']:
            config = set_null_ratios_for_matrices(config)

        # at this point the ratios are set one way or another. it is time to find out about the distribution of
        # components across matrices. note here we introduce a new member that was not in the original config class,
        # "num_components" per matrix.
        config = set_num_components_for_each_matrix(config)

        

        # now we know the exact number of components for each matrix. we can scale them to the expected number of
        # dimensions now.
        for matrix in config.matrices:
            m = config.matrices_dict[matrix]

            progress.new('Scaling matrix %d of %d (%s), for %d components' % (config.matrices.index(matrix) + 1,
                                                                              config.num_matrices,
                                                                              m['alias'],
                                                                              m['num_components']))

            m['scaled_vectors'] = utils.get_scaled_vectors(m['vectors'],
                                                           user_seed = config.seed,
                                                           n_components = m['num_components'],
                                                           normalize = m['normalize'],
                                                           progress=progress)

            progress.update('Normalizing scaled vectors ...')
            m['scaled_vectors'] = utils.get_normalized_vectors(m['scaled_vectors'])
            progress.end()


        # scaled vectors are in place. it is time to combine them to generate the input for final clustering
        progress.new('Scaled vectors for %d matrices' % config.num_matrices)
        progress.update('Combining ...')
        config.combined_vectors = []
        config.combined_id_to_sample = {}
        for i in range(0, len(config.master_rows)):
            row = config.master_rows[i]
            config.combined_id_to_sample[i] = config.master_rows[i]
            combined_scaled_vectors_for_row = [m['scaled_vectors'][m['sample_to_id'][row]] for m in config.matrices_dict.values()]
            config.combined_vectors.append(np.concatenate(combined_scaled_vectors_for_row))

        progress.update('Clustering ...')
        tree = utils.get_clustering_as_tree(config.combined_vectors, progress = progress)
        newick = utils.get_tree_object_in_newick(tree, config.combined_id_to_sample)
        progress.end()

        if config.output_file_path:
            open(config.output_file_path, 'w').write(newick + '\n')

        return newick

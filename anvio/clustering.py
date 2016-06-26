# -*- coding: utf-8
# pylint: disable=line-too-long
"""Clustering operations and helper functions"""

import os
import numpy as np
from sklearn import manifold
from sklearn import preprocessing
from scipy.cluster import hierarchy

import anvio
import anvio.utils as utils
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths

from anvio.errors import ConfigError
with terminal.SuppressAllOutput():
    from ete2 import Tree


__author__ = "A. Murat Eren"
__copyright__ = "Copyright 2015, The anvio Project"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "A. Murat Eren"
__email__ = "a.murat.eren@gmail.com"
__status__ = "Development"

run = terminal.Run()
progress = terminal.Progress()
progress.verbose = False


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


def depth(current_root):
    """
    current_root: represent for root, used in for loops.
    diff: distance between farthest_node and current_root.
    first_dist : in the beginning value for compare distances
    dist : distance between current node and root node in for loop.

    """

    diff = 0

    if current_root.get_sisters():
        node_sister = current_root.get_sisters()[0]
        farthest_node = current_root

        if current_root.get_children():
            farthest_node = current_root.get_children()[0]

        first_dist = farthest_node.get_distance(current_root, topology_only=True)

        # node variable is current root
        for cr_node in current_root.traverse("preorder"):
            # it's mean, all the cr_node children visited.
            if cr_node == node_sister:
                break

            # finding farthest_node to current_root
            else:
                dist = cr_node.get_distance(current_root, topology_only=True)

                if dist > first_dist:
                    farthest_node = cr_node
                    first_dist = dist

        diff = farthest_node.get_distance(current_root)

    return diff


def synchronize(current_root, control):
        # If the node modified in previous step. modified by node's sister.
        if not current_root.name.endswith(control):
            # this control catch root node.
            if current_root.get_sisters():
                new_root_sister = current_root.get_sisters()[0]

                current_node_depth = depth(current_root)
                sister_node_depth = depth(new_root_sister)

                if current_node_depth > sister_node_depth:
                    diff = current_node_depth - sister_node_depth
                    new_root_sister.dist += diff

                elif current_node_depth < sister_node_depth:
                    diff = sister_node_depth - current_node_depth
                    current_root.dist += diff

                # "_!$!" added into node names, cause this function modify the couples
                # but for loop in main function, send nodes one by one, not both of.
                current_root.name += control
                new_root_sister.name += control


def get_newick_tree_data_for_dict(d):
    matrix_file = filesnpaths.get_temp_file_path()
    utils.store_dict_as_TAB_delimited_file(d, matrix_file, ['items'] + d[d.keys()[0]].keys())
    newick = get_newick_tree_data(matrix_file)
    os.remove(matrix_file)
    return newick


def get_newick_tree_data(observation_matrix_path, output_file_name=None, method='ward', metric='euclidean',
                         norm='l1', progress=progress, transpose=False):
    filesnpaths.is_file_exists(observation_matrix_path)
    filesnpaths.is_file_tab_delimited(observation_matrix_path)

    if output_file_name:
        output_file_name = os.path.abspath(output_file_name)
        output_directory = os.path.dirname(output_file_name)
        if not os.access(output_directory, os.W_OK):
            raise ConfigError, "You do not have write permission for the output directory: '%s'" % output_directory

    id_to_sample_dict, sample_to_id_dict, header, vectors = utils.get_vectors_from_TAB_delim_matrix(observation_matrix_path, transpose=transpose)

    vectors = np.array(vectors)

    # normalize vectors:
    vectors = get_normalized_vectors(vectors, norm=norm, progress=progress)

    tree = get_clustering_as_tree(vectors, method, metric, progress)
    newick = get_tree_object_in_newick(tree, id_to_sample_dict)

    if output_file_name:
        open(output_file_name, 'w').write(newick.strip() + '\n')

    return newick


def get_scaled_vectors(vectors, user_seed=None, n_components=12, normalize=True, progress=progress):
    if user_seed:
        seed = np.random.RandomState(seed=user_seed)
    else:
        seed = np.random.RandomState()

    # FIXME: Make this optional:
    from sklearn.metrics.pairwise import euclidean_distances as d

    vectors = get_normalized_vectors(np.array(vectors)) if normalize else np.array(vectors)

    # compute similarities based on d
    progress.update('Computing similarity matrix')
    similarities = d(vectors)

    progress.update('Scaling using %d components' % n_components)
    mds = manifold.MDS(n_components=n_components, max_iter=300, eps=1e-10, random_state=seed,
                       dissimilarity="precomputed", n_jobs=1)

    progress.update('Fitting')
    scaled_vectors = mds.fit(similarities).embedding_

    return scaled_vectors


def get_normalized_vectors(vectors, norm='l1', progress=progress, pad_zeros=True):
    progress.update('Normalizing vectors using "%s" norm' % norm)
    vectors = np.array(vectors, dtype=np.float64)
    if pad_zeros:
        vectors += 0.0000001
    normalizer = preprocessing.Normalizer(norm=norm)
    return normalizer.fit_transform(vectors)


def get_clustering_as_tree_obsolete(vectors, ward=True, clustering_distance='euclidean', clustering_method='complete', progress=progress):
    try:
        import hcluster
    except ImportError:
        raise ConfigError, "This is an obsolete function requires an obsolete module. If you\
                            still want to play with it, you can. But you've got to get hcluster\
                            installed. Run this in your terminal, and make sure you get no errors:\
                            python -c 'import hcluster'"

    if ward:
        progress.update('Clustering data with Ward linkage and euclidean distances')
        clustering_result = hcluster.ward(vectors)
    else:
        progress.update('Computing distance matrix using "%s" distance' % clustering_distance)
        distance_matrix = hcluster.pdist(vectors, clustering_distance)
        progress.update('Clustering data with "%s" linkage' % clustering_method)
        clustering_result = hcluster.linkage(distance_matrix, method=clustering_method)

    progress.update('Returning results')
    return hcluster.to_tree(clustering_result)


def get_tree_object_in_newick_obsolete(tree, id_to_sample_dict=None):
    """Take a tree object, and create a newick formatted representation of it"""

    root = Tree()
    root.dist = 0
    root.name = "root"
    item2node = {tree: root}

    to_visit = [tree]
    while to_visit:
        node = to_visit.pop()
        cl_dist = node.dist / 2.0
        for ch_node in [node.left, node.right]:
            if ch_node:
                ch = Tree()
                ch.dist = cl_dist

                if ch_node.is_leaf():
                    if id_to_sample_dict:
                        ch.name = id_to_sample_dict[ch_node.id]
                    else:
                        ch.name = ch_node.id
                else:
                    ch.name = 'Int' + str(ch_node.id)

                item2node[node].add_child(ch)
                item2node[ch_node] = ch
                to_visit.append(ch_node)

    return root.write(format=1)


def get_clustering_as_tree(vectors, method='ward', metric='euclidean', progress=progress):
    progress.update('Clustering data with "%s" linkage using "%s" distance' % (method, metric))
    linkage = hierarchy.linkage(vectors, method=method, metric=metric)

    progress.update('Recovering the tree from the clustering result')
    tree = hierarchy.to_tree(linkage, rd=False)

    return tree


def get_tree_object_in_newick(tree, id_to_sample_dict=None):
    """Take a tree object, and create a newick formatted representation of it"""

    new_tree = Tree()
    new_tree.dist = 0
    new_tree.name = "root"

    node_id = 0
    node_id_to_node_in_old_tree = {node_id: tree}
    node_id_to_node_in_new_tree = {node_id: new_tree}

    node_ids_to_visit_in_old_tree = [node_id]

    while node_ids_to_visit_in_old_tree:
        node_id_in_old_tree = node_ids_to_visit_in_old_tree.pop()
        node_in_old_tree = node_id_to_node_in_old_tree[node_id_in_old_tree]
        cl_dist = node_in_old_tree.dist / 2.0

        for ch_node_in_old_tree in [node_in_old_tree.left, node_in_old_tree.right]:
            if ch_node_in_old_tree:
                ch_for_new_tree = Tree()
                ch_for_new_tree.dist = cl_dist

                node_id += 1
                node_id_to_node_in_new_tree[node_id] = ch_for_new_tree

                if ch_node_in_old_tree.is_leaf():
                    if id_to_sample_dict:
                        ch_for_new_tree.name = id_to_sample_dict[ch_node_in_old_tree.id]
                    else:
                        ch_for_new_tree.name = ch_node_in_old_tree.id
                else:
                    ch_for_new_tree.name = 'Int' + str(ch_node_in_old_tree.id)

                node_id_to_node_in_new_tree[node_id_in_old_tree].add_child(ch_for_new_tree)
                node_id_to_node_in_old_tree[node_id] = ch_node_in_old_tree
                node_ids_to_visit_in_old_tree.append(node_id)

    return new_tree.write(format=1)


def order_contigs_simple(config, progress=progress, run=run, debug=False):
    if not config.matrices_dict[config.matrices[0]]['ratio']:
        config = set_null_ratios_for_matrices(config)

    if debug:
        run.info_single('Peak at the first 5 items in the first 5 rows in matrices:', mc='green', nl_before=2)

    for matrix in config.matrices:
        m = config.matrices_dict[matrix]

        m['scaled_vectors'] = np.array(m['vectors'], dtype=np.float64)

        if m['normalize']:
            m['scaled_vectors'] = get_normalized_vectors(m['scaled_vectors'])

        if m['log']:
            m['scaled_vectors'] = np.log10(m['scaled_vectors'] + 1)

        if debug:
            summary = '\n'.join(['%s (...)' % m['scaled_vectors'][i][0:5] for i in range(0, 5)])
            run.warning(summary, 'Vectors for "%s" (%d by %d)' % (matrix, len(m['scaled_vectors']), len(m['scaled_vectors'][0])), lc='crimson', raw=True)

    progress.new('Vectors from %d matrices' % config.num_matrices)
    progress.update('Combining ...')
    config.combined_vectors = []
    config.combined_id_to_sample = {}

    for i in range(0, len(config.master_rows)):
        row = config.master_rows[i]
        config.combined_id_to_sample[i] = config.master_rows[i]
        combined_scaled_vectors_for_row = [m['scaled_vectors'][m['sample_to_id'][row]] for m in config.matrices_dict.values()]
        config.combined_vectors.append(np.concatenate(combined_scaled_vectors_for_row))

    progress.update('Clustering ...')
    tree = get_clustering_as_tree(config.combined_vectors, progress=progress)
    newick = get_tree_object_in_newick(tree, config.combined_id_to_sample)
    progress.end()

    if config.output_file_path:
        open(config.output_file_path, 'w').write(newick + '\n')

    return newick


def order_contigs_experimental(config, progress=progress, run=run, debug=False):
    if not config.multiple_matrices:
        # there is one matrix. could be coverage, could be tnf. we don't care.
        # we do what we gotta do: skip scaling and perform clustering using all
        # dimensions.
        m = config.matrices_dict[config.matrices[0]]

        progress.new('Single matrix (%s)' % m['alias'])
        progress.update('Performing cluster analysis ...')
        tree = get_clustering_as_tree(m['vectors'], progress=progress)
        newick = get_tree_object_in_newick(tree, m['id_to_sample'])
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

            m['scaled_vectors'] = get_scaled_vectors(m['vectors'],
                                                           user_seed=config.seed,
                                                           n_components=m['num_components'],
                                                           normalize=m['normalize'],
                                                           progress=progress)

            progress.update('Normalizing scaled vectors ...')
            m['scaled_vectors'] = get_normalized_vectors(m['scaled_vectors'])
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
        tree = get_clustering_as_tree(config.combined_vectors, progress=progress)
        newick = get_tree_object_in_newick(tree, config.combined_id_to_sample)
        progress.end()

        if config.output_file_path:
            open(config.output_file_path, 'w').write(newick + '\n')

        return newick

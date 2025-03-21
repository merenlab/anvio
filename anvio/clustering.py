# -*- coding: utf-8
# pylint: disable=line-too-long
"""Clustering operations and helper functions"""

import numpy as np
import pandas as pd

from sklearn import manifold
from sklearn import preprocessing
from scipy.cluster import hierarchy
from scipy.spatial import distance as scipy_distance
from scipy.spatial.distance import pdist, squareform

import anvio
import anvio.utils as utils
import anvio.terminal as terminal
import anvio.constants as constants
import anvio.filesnpaths as filesnpaths

from anvio.errors import ConfigError
with terminal.SuppressAllOutput():
    from ete3 import Tree


distance_metrics = ['euclidean', 'cityblock', 'sqeuclidean', 'cosine', 'correlation', 'hamming',\
                    'hamming', 'jaccard', 'chebyshev', 'canberra', 'braycurtis', 'yule', 'matching',\
                    'dice', 'kulsinski', 'rogerstanimoto', 'russellrao', 'braycurtis', 'yule',\
                    'matching', 'dice', 'kulsinski', 'rogerstanimoto', 'russellrao', 'sokalmichener',\
                    'sokalsneath', 'minkowski']

linkage_methods = ['single', 'complete', 'average', 'weighted', 'centroid', 'median', 'ward']


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
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
    denominator = float(sum([r['ratio'] for r in list(config.matrices_dict.values())]))
    for matrix in config.matrices:
        m = config.matrices_dict[matrix]
        num_components_for_ratio = int(round(config.num_components * (m['ratio'] / denominator)))
        m['num_components'] = num_components_for_ratio

    return config


def set_null_ratios_for_matrices(config):
    for matrix in config.matrices:
        config.matrices_dict[matrix]['ratio'] = 1
    return config


def is_distance_and_linkage_compatible(distance, linkage):
    is_linkage_method_OK(linkage)
    is_distance_metric_OK(distance)

    if distance == 'yule' and linkage != 'single':
        raise ConfigError("The distance metric 'yule' will only work with the linkage 'single' :/")

    try:
        hierarchy.linkage([(1, 0), (0, 1), (1, 1)], metric=distance, method=linkage)
    except Exception as exception:
        raise ConfigError("Someone is upset here: %s" % exception)


def is_linkage_method_OK(linkage):
    if linkage not in linkage_methods:
        raise ConfigError("Linkage '%s' is not one of the linkage methods anvi'o recognizes :/ Here "
                           "is a list of all the available ones: %s" % (linkage, ', '.join(linkage_methods)))


def is_distance_metric_OK(distance):
    if distance not in distance_metrics:
        raise ConfigError("Distance '%s' is not one of the metrics anvi'o recognizes :/ Here "
                           "is a list of all the available ones: %s" % (distance, ', '.join(distance_metrics)))


def get_newick_tree_data_for_dict(d, transpose=False, linkage=constants.linkage_method_default, distance=constants.distance_metric_default, norm='l1', zero_fill_missing=False):
    is_distance_and_linkage_compatible(distance, linkage)

    if zero_fill_missing:
        vectors = pd.DataFrame.from_dict(d, orient='index').fillna(0)
    else:
        vectors = pd.DataFrame.from_dict(d, orient='index')

    id_to_sample_dict = dict([(i, vectors.index[i]) for i in range(len(vectors.index))])

    if transpose:
        id_to_sample_dict = dict([(i, vectors.columns[i]) for i in range(len(vectors.columns))])

    newick = get_newick_from_matrix(vectors, distance, linkage, norm, id_to_sample_dict, transpose=transpose)

    return newick


def get_newick(node, parent_dist, leaf_names, newick=''):
    """
    Modified from the solution at https://stackoverflow.com/questions/28222179/save-dendrogram-to-newick-format
    """
    if node.is_leaf():
        return "%s:%.2f%s" % (leaf_names[node.id], parent_dist - node.dist, newick)
    else:
        if len(newick) > 0:
            newick = "):%.2f%s" % (parent_dist - node.dist, newick)
        else:
            newick = ");"
        newick = get_newick(node.get_left(), node.dist, leaf_names, newick=newick)
        newick = get_newick(node.get_right(), node.dist, leaf_names, newick=",%s" % (newick))
        newick = "(%s" % (newick)
        return newick


def get_vectors_for_vectors_with_missing_data(vectors):
    """Get a distance matrix for vectors with missing data.

    Modified from the solution at https://stackoverflow.com/questions/31420912/python-hierarchical-clustering-with-missing-values
    """

    vectors[vectors == None] = np.nan

    Npat = vectors.shape[0]
    dist = np.ndarray(shape=(Npat,Npat))
    dist.fill(0)
    for ix in range(0,Npat):
        x = vectors[ix,]
        if ix >0:
            for iy in range(0,ix):
                y = vectors[iy,]
                dist[ix,iy] = np.nanmean((x - y) ** 2, dtype='float32')
                dist[iy,ix] = dist[ix,iy]

    return scipy_distance.squareform(dist)


def get_newick_from_matrix(vectors, distance, linkage, norm, id_to_sample_dict, transpose=False):
    is_distance_and_linkage_compatible(distance, linkage)

    if transpose:
        vectors = vectors.transpose()

    # if there are missing values in data, we will assume that
    # it was normalized by the user.
    if None in vectors:
        vectors = get_vectors_for_vectors_with_missing_data(vectors)
    else:
        # normalize vectors:
        vectors = get_normalized_vectors(vectors, norm=norm, progress=progress)

    tree = get_clustering_as_tree(vectors, linkage, distance, progress)
    newick = get_tree_object_in_newick(tree, id_to_sample_dict)

    return newick


def create_newick_file_from_matrix_file(observation_matrix_path, output_file_path, linkage=constants.linkage_method_default,
                         distance=constants.distance_metric_default, norm='l1', progress=progress, transpose=False,
                         items_order_file_path=None, pad_with_zeros=False, distance_matrix_output_path=None):
    is_distance_and_linkage_compatible(distance, linkage)
    filesnpaths.is_file_exists(observation_matrix_path)
    filesnpaths.is_file_tab_delimited(observation_matrix_path)
    filesnpaths.is_output_file_writable(output_file_path)

    if distance_matrix_output_path:
        filesnpaths.is_output_file_writable(distance_matrix_output_path)

    if items_order_file_path:
        filesnpaths.is_output_file_writable(items_order_file_path)

    id_to_sample_dict, sample_to_id_dict, header, vectors = utils.get_vectors_from_TAB_delim_matrix(observation_matrix_path, transpose=transpose, pad_with_zeros=pad_with_zeros)

    vectors = np.array(vectors)

    newick = get_newick_from_matrix(vectors, distance, linkage, norm, id_to_sample_dict)

    if distance_matrix_output_path:
        # the user has requested a distance matrix calculated from the
        # same data that yields the newick tree to also be reported.
        distance_matrix = get_distance_matrix(vectors, distance=distance)
        report_distance_matrix(distance_matrix, id_to_sample_dict, distance_matrix_output_path)

    if output_file_path:
        open(output_file_path, 'w').write(newick.strip() + '\n')

    if items_order_file_path:
        open(items_order_file_path, 'w').write('\n'.join(utils.get_names_order_from_newick_tree(newick)) + '\n')


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


def report_distance_matrix(distance_matrix, id_to_sample_dict, distance_matrix_output_path):
    with open(distance_matrix_output_path, 'w') as output:
        line = ['items'] + [id_to_sample_dict[i] for i in range(0, len(distance_matrix))]
        output.write('\t'.join(line) + '\n')
        for i in range(0, len(distance_matrix)):
            line = [id_to_sample_dict[i]] + [str(v) for v in distance_matrix[i].tolist()]
            output.write('\t'.join(line) + '\n')

    run.info('Distance matrix', distance_matrix_output_path)


def get_distance_matrix(vectors, distance=constants.distance_metric_default, progress=progress):
    progress.update("Calculating the distance marix")

    distances = pdist(vectors, metric=distance)
    matrix = squareform(distances)

    progress.clear()

    return matrix


def get_clustering_as_tree(vectors, linkage=constants.linkage_method_default, distance=constants.distance_metric_default, progress=progress):
    is_distance_and_linkage_compatible(distance, linkage)

    progress.update('Clustering data with "%s" linkage using "%s" distance' % (linkage, distance))
    linkage = hierarchy.linkage(vectors, metric=distance, method=linkage)

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
                    # we used to export our trees with internal node labels so we could
                    # do various interface operations more easily:
                    #
                    #    ch_for_new_tree.name = 'Int' + str(ch_node_in_old_tree.id)
                    #
                    # but our new interface design does not require such addditions to
                    # dendrograms. Although here we add 0 branch support for our
                    # dendrograms since we wish to use a standard format to export these
                    # data as a tree.
                    ch_for_new_tree.support = 0.0

                node_id_to_node_in_new_tree[node_id_in_old_tree].add_child(ch_for_new_tree)
                node_id_to_node_in_old_tree[node_id] = ch_node_in_old_tree
                node_ids_to_visit_in_old_tree.append(node_id)

    for node in new_tree.traverse("preorder"):
        if node.is_leaf():
            continue

        has_child_with_dist_or_int = False

        for child in node.get_children():
            if not child.is_leaf() or child.dist > 0:
                has_child_with_dist_or_int = True
                break

        if has_child_with_dist_or_int:
            continue

        # swap childs alphabetically
        node.children = sorted(node.get_children(), key=lambda x:x.name, reverse=True)

    return new_tree.write(format=2)


def order_contigs_simple(config, distance=None, linkage=None, progress=progress, run=run, debug=False):
    """An anvi'o clustering config comes in, a (clustering_id, newick) tuple goes out.

       By default the `linkage` and `distance` is set to the system defaults, constants.linkage_method_default
       and constants.distance_metric_default. If the `config` has either of them defined, the system defaults
       are overwritten with the preference in the config file. If the function gets `linkage` or `distance` as
       parameter, they overwrite both system defaults and config preferences.
    """

    if not config.matrices_dict[config.matrices[0]]['ratio']:
        config = set_null_ratios_for_matrices(config)

    distance = distance if distance else (config.distance or constants.distance_metric_default)
    linkage = linkage if linkage else (config.linkage or constants.linkage_method_default)
    clustering_id = ':'.join([config.name, distance, linkage])

    if len(config.master_rows) == 1:
        # there is a single item to cluster. which means there is nothing to cluster really.
        # return that single item in a newick format:
        return (clustering_id, '(%s);' % config.master_rows[0])

    if debug or anvio.DEBUG:
        run.info_single('Peak at the first 5 items in the first 5 rows in matrices:', mc='green', nl_before=2)

    for matrix in config.matrices:
        m = config.matrices_dict[matrix]

        m['scaled_vectors'] = np.array(m['vectors'], dtype=np.float64)

        if m['normalize']:
            m['scaled_vectors'] = get_normalized_vectors(m['scaled_vectors'])

        if m['log']:
            m['scaled_vectors'] = np.log10(m['scaled_vectors'] + 1)

        if debug or anvio.DEBUG:
            summary = '\n'.join(['%s (...)' % m['scaled_vectors'][i][0:5] for i in range(0, 5)])
            run.warning(summary, 'Vectors for "%s" (%d by %d)' % (matrix, len(m['scaled_vectors']), len(m['scaled_vectors'][0])), lc='crimson', raw=True)

    progress.new('Vectors from %d matrices' % config.num_matrices)
    progress.update('Combining ...')
    config.combined_vectors = []
    config.combined_id_to_sample = {}

    for i in range(0, len(config.master_rows)):
        row = config.master_rows[i]
        config.combined_id_to_sample[i] = config.master_rows[i]
        combined_scaled_vectors_for_row = [m['scaled_vectors'][m['sample_to_id'][row]] for m in list(config.matrices_dict.values())]
        config.combined_vectors.append(np.concatenate(combined_scaled_vectors_for_row))

    progress.update('Clustering ...')

    tree = get_clustering_as_tree(config.combined_vectors, linkage, distance, progress=progress)
    newick = get_tree_object_in_newick(tree, config.combined_id_to_sample)

    progress.end()

    if config.output_file_path:
        open(config.output_file_path, 'w').write(newick + '\n')

    return (clustering_id, newick)


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
            combined_scaled_vectors_for_row = [m['scaled_vectors'][m['sample_to_id'][row]] for m in list(config.matrices_dict.values())]
            config.combined_vectors.append(np.concatenate(combined_scaled_vectors_for_row))

        progress.update('Clustering ...')
        tree = get_clustering_as_tree(config.combined_vectors, progress=progress)
        newick = get_tree_object_in_newick(tree, config.combined_id_to_sample)
        progress.end()

        if config.output_file_path:
            open(config.output_file_path, 'w').write(newick + '\n')

        return newick

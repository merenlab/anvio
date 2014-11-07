# -*- coding: utf-8
#
# Copyright (C) 2010 - 2012, A. Murat Eren
#
# This program is free software; you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free
# Software Foundation; either version 2 of the License, or (at your option)
# any later version.
#
# Please read the COPYING file.

import os
import sys
import time
import json
import socket
import shutil
import hcluster
import textwrap
import itertools
import subprocess
import multiprocessing

from ete2 import Tree
import numpy as np
from sklearn import manifold
from sklearn import preprocessing

from PaPi.constants import IS_ESSENTIAL_FIELD, IS_AUXILIARY_FIELD, allowed_chars, complements, levels_of_taxonomy
import PaPi.fastalib as u
import PaPi.filesnpaths as filesnpaths
from PaPi.contig import Composition

from PaPi.terminal import Progress

# Mock progress object that will not report anything, for general clarity.
progress = Progress()
progress.verbose = False


# absolute path anonymous:
ABS = lambda x: x if x.startswith('/') else os.path.join(os.getcwd(), x)


def rev_comp(seq):
    return seq.translate(complements)[::-1]


class Multiprocessing:
    def __init__(self, target_function, num_thread = None):
        self.cpu_count = multiprocessing.cpu_count()
        self.num_thread = num_thread or (self.cpu_count - (int(round(self.cpu_count / 10.0)) or 1))
        self.target_function = target_function
        self.processes = []
        self.manager = multiprocessing.Manager()


    def get_data_chunks(self, data_array, spiral = False):
        data_chunk_size = (len(data_array) / self.num_thread) or 1
        data_chunks = []
        
        if len(data_array) <= self.num_thread:
            return [[chunk] for chunk in data_array]

        if spiral:
            for i in range(0, self.num_thread):
                data_chunks.append([data_array[j] for j in range(i, len(data_array), self.num_thread)])
            
            return data_chunks
        else:
            for i in range(0, self.num_thread):
                if i == self.num_thread - 1:
                    data_chunks.append(data_array[i * data_chunk_size:])
                else:
                    data_chunks.append(data_array[i * data_chunk_size:i * data_chunk_size + data_chunk_size])

        return data_chunks

                
    def run(self, args, name = None):
        t = multiprocessing.Process(name = name,
                                    target = self.target_function,
                                    args = args)
        self.processes.append(t)
        t.start()


    def get_empty_shared_array(self):
        return self.manager.list()


    def get_empty_shared_dict(self):
        return self.manager.dict()

    
    def get_shared_integer(self):
        return self.manager.Value('i', 0)


    def run_processes(self, processes_to_run, progress = Progress(verbose=False)):
        tot_num_processes = len(processes_to_run)
        sent_to_run = 0
        while 1:
            NumRunningProceses = lambda: len([p for p in self.processes if p.is_alive()])
            
            if NumRunningProceses() < self.num_thread and processes_to_run:
                for i in range(0, self.num_thread - NumRunningProceses()):
                    if len(processes_to_run):
                        sent_to_run += 1
                        self.run(processes_to_run.pop())

            if not NumRunningProceses() and not processes_to_run:
                # let the blastn program finish writing all output files.
                # FIXME: this is ridiculous. find a better solution.
                time.sleep(1)
                break

            progress_obj.update('%d of %d done in %d threads (currently running processes: %d)'\
                                                         % (sent_to_run - NumRunningProceses(),
                                                            tot_num_processes,
                                                            self.num_thread,
                                                            NumRunningProceses()))
            time.sleep(1)


def get_available_port_num(start = 8080, look_upto_next_num_ports = 100):
    """Starts from 'start' and incrementally looks for an available port
       until 'start + look_upto_next_num_ports', and returns the first
       available one."""
    for p in range(start, start + look_upto_next_num_ports):
        if not is_port_in_use(p):
            return p

    return None


def is_port_in_use(port):
    in_use = False
    sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    result = sock.connect_ex(('127.0.0.1', port))

    if(result == 0) :
        in_use = True

    sock.close()
    return in_use


def run_command(cmdline):
    try:
        if subprocess.call(cmdline, shell = True) < 0:
            raise ConfigError, "command was terminated by signal: %d" % (-retcode)
    except OSError, e:
        raise ConfigError, "command was failed for the following reason: '%s' ('%s')" % (e, cmdline)


def store_dict_as_TAB_delimited_file(d, output_path, headers):
    filesnpaths.is_output_file_writable(output_path)
    f = open(output_path, 'w')
    f.write('%s\n' % '\t'.join(headers))
    for k in d.keys():
        line = [k]
        for header in headers[1:]:
            line.append(d[k][header])
        f.write('%s\n' % '\t'.join(line))
    return output_path


def is_all_columns_present_in_TAB_delim_file(columns, file_path):
    header_fields = get_header_fields_of_TAB_delim_file(file_path)
    return False if len([False for c in columns if c not in header_fields]) else True


def get_header_fields_of_TAB_delim_file(file_path):
    return open(file_path).readline().strip('\n').split('\t')[1:]


def get_json_obj_from_TAB_delim_metadata(input_file):
    return json.dumps([line.strip('\n').split('\t') for line in open(input_file).readlines()])


def get_vectors_from_TAB_delim_matrix(file_path, fields_to_return=None):
    filesnpaths.is_file_exists(file_path)
    filesnpaths.is_file_tab_delimited(file_path)

    vectors = []
    id_to_sample_dict = {}

    input_matrix = open(file_path)
    header = input_matrix.readline().strip().split('\t')[1:]

    fields_of_interest = []
    if fields_to_return:
        fields_of_interest = [f for f in range(0, len(header)) if header[f] in fields_to_return and IS_ESSENTIAL_FIELD(header[f])]
    else:
        fields_of_interest = [f for f in range(0, len(header)) if IS_ESSENTIAL_FIELD(header[f])]

    # update header:
    header = [header[i] for i in range(0, len(header)) if i in fields_of_interest]


    if not len(header):
        raise ConfigError, "Only a subset (%d) of fields were requested by the caller, but none of them was found\
                            in the matrix (%s) :/" % (len(fields_to_return), file_path)

    line_counter = 0
    for line in input_matrix.readlines():
        id_to_sample_dict[line_counter] = line.strip().split('\t')[0]
        fields = line.strip().split('\t')[1:]

        if fields_of_interest:
            vector = [float(fields[i]) for i in range(0, len(fields)) if i in fields_of_interest]
        else:
            vector = [float(f) for f in fields]

        vectors.append(vector)

        line_counter += 1

    return id_to_sample_dict, header, vectors


def get_all_ids_from_fasta(input_file):
    fasta = u.SequenceSource(input_file)
    ids = []
    
    while fasta.next():
        ids.append(fasta.id) 

    return ids


def get_GC_content_for_FASTA_entries(file_path):
    filesnpaths.is_file_exists(file_path)
    filesnpaths.is_file_fasta_formatted(file_path)

    GC_content_dict = {}

    fasta = u.SequenceSource(file_path)
    while fasta.next():
        GC_content_dict[fasta.id] = get_GC_content_for_sequence(fasta.seq)

    return GC_content_dict


def get_GC_content_for_sequence(sequence):
    return Composition(sequence).GC_content


def get_cmd_line():
    c_argv = []
    for i in sys.argv:
        if ' ' in i:
            c_argv.append('"%s"' % i)
        else:
            c_argv.append(i)
    return ' '.join(c_argv)


def concatenate_files(dest_file, file_list):
    if not dest_file:
        raise ConfigError, "Destination cannot be empty."
    if not len(file_list):
        raise ConfigError, "File list cannot be empty."
    for f in file_list:
        filesnpaths.is_file_exists(f)
    filesnpaths.is_output_file_writable(dest_file)

    dest_file_obj = open(dest_file, 'w')
    for chunk_path in file_list:
        for line in open(chunk_path):
            dest_file_obj.write(line)

    dest_file_obj.close()
    return dest_file


class ConfigError(Exception):
    def __init__(self, e = None):
        Exception.__init__(self)
        while 1:
            if e.find("  ") > -1:
                e = e.replace("  ", " ")
            else:
                break
        self.e = e
        return
    def __str__(self):
        error_type = 'Config Error'

        max_len = max([len(l) for l in textwrap.fill(self.e, 80).split('\n')])
        error_lines = ['\033[0;30m\033[46m%s%s\033[0m' % (l, ' ' * (max_len - len(l)))\
                                         for l in textwrap.fill(self.e, 80).split('\n')]

        error_message = ['%s: %s' % (error_type, error_lines[0])]
        for error_line in error_lines[1:]:
            error_message.append('%s%s' % (' ' * (len(error_type) + 2), error_line))

        return '\n'.join(error_message)


def get_chunks(contig_length, desired_length):
    num_chunks = contig_length / desired_length

    if num_chunks < 2:
        return [(0, contig_length)]

    chunks = []
    for i in range(0, num_chunks):
        chunks.append((i * desired_length, (i + 1) * desired_length),)
    chunks.append(((i + 1) * desired_length, contig_length),)

    if (chunks[-1][1] - chunks[-1][0]) < (desired_length / 2):
        # last chunk is too small :/ merge it to the previous one.
        last_tuple = (chunks[-2][0], contig_length)
        chunks.pop()
        chunks.pop()
        chunks.append(last_tuple)

    return chunks


def check_sample_id(sample_id):
    if sample_id:
        if len([c for c in sample_id if c not in allowed_chars]):
            raise ConfigError, "Project name ('%s') contains characters that PaPi does not like. Please\
                                limit the characters that make up the project name to ASCII letters,\
                                digits, '_' and '-' (if you had not declared a project name and PaPi made\
                                up one for you, please specify with '-p' parameter specifically)." % sample_id


def check_contig_names(contig_names):
    all_characters_in_contig_names = set(''.join(contig_names))
    characters_PaPi_doesnt_like = [c for c in all_characters_in_contig_names if c not in allowed_chars]
    if len(characters_PaPi_doesnt_like):
        raise ConfigError, "The name of at least one contig in your BAM file %s PaPi does not\
                            like (%s). Please go back to your original files and make sure that\
                            the characters in contig names are limited to to ASCII letters,\
                            digits. Names can also contain underscore ('_'), dash ('-') and dot ('.')\
                            characters. PaPi knows how much work this may require for you to go back and\
                            re-generate your BAM files and is very sorry for asking you to do that, however,\
                            it is critical for later steps in the analysis." \
                                % ("contains multiple characters" if len(characters_PaPi_doesnt_like) > 1 else "contains a character",
                                   ", ".join(['"%s"' % c for c in characters_PaPi_doesnt_like]))


def get_TAB_delimited_file_as_dictionary(file_path, expected_fields = None):
    filesnpaths.is_file_tab_delimited(file_path)

    f = open(file_path)

    header_fields = f.readline().strip('\n').split('\t')

    if expected_fields:
        for field in expected_fields:
            if field not in header_fields:
                raise ConfigError, "The file '%s' does not contain the right type of header. It was expected\
                                    to have these: '%s', however it had these: '%s'" % (file_path,
                                                                                        ', '.join(expected_fields),
                                                                                        ', '.join(header_fields[1:]))

    d = {}

    for line in f.readlines():
        line_fields = line.strip('\n').split('\t')

        d[line_fields[0]] = {}
        e = d[line_fields[0]]
        for i in range(1, len(header_fields)):
            e[header_fields[i]] = line_fields[i]

    return d


def get_newick_tree_data(observation_matrix_path, output_file_name = None, clustering_distance='euclidean',
                         clustering_method = 'complete', norm = 'l1', progress = Progress(verbose=False)):
    filesnpaths.is_file_exists(observation_matrix_path)
    filesnpaths.is_file_tab_delimited(observation_matrix_path)

    if output_file_name:
        output_file_name = ABS(output_file_name)
        output_directory = os.path.dirname(output_file_name)
        if not os.access(output_directory, os.W_OK):
            raise ConfigError, "You do not have write permission for the output directory: '%s'" % output_directory
    
    id_to_sample_dict, header, vectors = get_vectors_from_TAB_delim_matrix(observation_matrix_path)

    vectors = np.array(vectors)

    # normalize vectors:
    vectors = get_normalized_vectors(vectors, norm=norm, progress=progress)

    tree = get_clustering_as_tree(vectors, clustering_distance, clustering_method, progress)
    newick = get_tree_object_in_newick(tree, id_to_sample_dict)
   
    if output_file_name:
        open(output_file_name, 'w').write(newick.strip() + '\n')

    return newick


def get_scaled_vectors(vectors, user_seed = None, n_components = 12, normalize=True, progress = Progress(verbose=False)):
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


def get_normalized_vectors(vectors, norm='l1', progress = Progress(verbose=False)):
    progress.update('Normalizing vectors using "%s" norm' % norm)
    normalizer = preprocessing.Normalizer(norm=norm)
    return normalizer.fit_transform(vectors)


def get_clustering_as_tree(vectors, clustering_distance='euclidean', clustering_method = 'complete', progress = Progress(verbose=False)):
    progress.update('Computing distance matrix using "%s" distance' % clustering_distance)
    distance_matrix = hcluster.pdist(vectors, clustering_distance)
    progress.update('Clustering data with "%s" linkage' % clustering_method)
    clustering_result = hcluster.linkage(distance_matrix, method = clustering_method)
    progress.update('Returning results')
    return hcluster.to_tree(clustering_result)


def get_tree_object_in_newick(tree, id_to_sample_dict):
    """i.e., tree = hcluster.to_tree(c_res)"""

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
                    ch.name = id_to_sample_dict[ch_node.id]
                else:
                    ch.name = 'Int' + str(ch_node.id)

                item2node[node].add_child(ch)
                item2node[ch_node] = ch
                to_visit.append(ch_node)

    return root.write(format=1) 


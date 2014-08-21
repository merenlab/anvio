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
import fcntl
import socket
import string
import struct
import cPickle
import termios 
import hcluster
import textwrap
import tempfile
import itertools
import subprocess
import multiprocessing

from ete2 import Tree

from PaPi.constants import pretty_names
import PaPi.fastalib as u


complements = string.maketrans('acgtrymkbdhvACGTRYMKBDHV',\
                               'tgcayrkmvhdbTGCAYRKMVHDB')
levels_of_taxonomy = ["phylum", "class", "order", "family", "genus", "species"]

# absolute path anonymous:
ABS = lambda x: x if x.startswith('/') else os.path.join(os.getcwd(), x)

def rev_comp(seq):
    return seq.translate(complements)[::-1]


class KMers:
    def __init__(self, k = 4):
        self.kmers = {}
        self.k = k
        
        self.get_kmers()

    def get_kmers(self):
        k = self.k
        arg = ['ATCG'] * k
        kmers = set()
        
        for item in itertools.product(*arg):
            kmer = ''.join(item)
            if rev_comp(kmer) not in kmers:
                kmers.add(kmer)
        
        self.kmers[k] = kmers


    def get_kmer_frequency(self, sequence, dist_metric_safe = True):
        k = self.k
        sequence = sequence.upper()

        if len(sequence) < k:
            return None

        if not self.kmers.has_key(k):
            self.get_kmers(k)
        
        kmers = self.kmers[k]
        frequencies = dict(zip(kmers, [0] * len(kmers)))
        
        for i in range(0, len(sequence) - (k - 1)):
            kmer = sequence[i:i + k]
            
            # FIXME: this can be faster/better
            if len([n for n in kmer if n not in 'ATCG']):
                continue

            if frequencies.has_key(kmer):
                frequencies[kmer] += 1
            else:
                frequencies[rev_comp(kmer)] += 1

        if dist_metric_safe:
            # we don't want all kmer freq values to be zero. so the distance
            # metrics wouldn't go crazy. instead we fill it with 1. which
            # doesn't affect relative distances.
            if sum(frequencies.values()) == 0:
                words = self.kmers[self.k]
                frequencies = dict(zip(words, [1] * len(words)))

        return frequencies


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


    def run_processes(self, processes_to_run, progress_obj = None):
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

            if progress_obj:
                progress_obj.update('%d of %d done in %d threads (currently running processes: %d)'\
                                                         % (sent_to_run - NumRunningProceses(),
                                                            tot_num_processes,
                                                            self.num_thread,
                                                            NumRunningProceses()))
            time.sleep(1)


class Progress:
    def __init__(self):
        self.pid = None
        self.verbose = True
        self.terminal_width = None

        self.get_terminal_width()
        self.color_prefix = '\033[0;30m\033[46m'
        self.color_postfix = '\033[0m'
        
        self.currently_shown = None


    def get_terminal_width(self):
        try:
            self.terminal_width = get_terminal_size()[0]
        except:
            self.terminal_width = 80 


    def new(self, pid):
        if self.pid:
            raise ConfigError, "Progress.new() can't be called before ending the previous one (Existing: '%s', Competing: '%s')." % (self.pid, pid)

        if not self.verbose:
            return

        self.pid = '%s %s' % (get_date(), pid)
        self.get_terminal_width()
        self.currently_shown = None


    def write(self, c):
        surpass = self.terminal_width - len(c)
        
        if surpass < 0:
            c = c[0:-(-surpass + 4)] + ' (...)'
        else:
            self.currently_shown = c
            c = c + ' ' * surpass

        if self.verbose:
            sys.stderr.write(self.color_prefix + c + self.color_postfix)
            sys.stderr.flush()


    def reset(self):
        self.clear()

    def clear(self):
        if not self.verbose:
            return
        null = '\r' + ' ' * (self.terminal_width) 
        sys.stderr.write(null)
        sys.stderr.write('\r')
        sys.stderr.flush()
        self.currently_shown = None


    def append(self, msg):
        if not self.verbose:
            return
        self.write('%s%s' % (self.currently_shown, msg))


    def update(self, msg):
        if not self.verbose:
            return
        self.clear()
        self.write('\r[%s] %s' % (self.pid, msg))


    def end(self):
        self.pid = None
        if not self.verbose:
            return
        self.clear()


class Run:
    def __init__(self, info_file_path = None, verbose = True):
        if info_file_path:
            self.init_info_file_obj(info_file_path)
        else:
            self.info_file_obj = None

        self.info_dict = {}
        self.verbose = verbose


    def init_info_file_obj(self, info_file_path):
            self.info_file_obj = open(info_file_path, 'w')


    def info(self, key, value, quiet = False):
        self.info_dict[key] = value
        
        if quiet:
            return True
        
        if type(value) == int:
            value = pretty_print(value)

        label = get_pretty_name(key)

        info_line = "%s %s: %s\n" % (label, '.' * (65 - len(label)), str(value))
        if self.info_file_obj:
            self.info_file_obj.write(info_line)

        if self.verbose:
            sys.stderr.write(info_line)


    def store_info_dict(self, destination):
        cPickle.dump(self.info_dict, open(destination, 'w'))


    def quit(self):
        if self.info_file_obj:
            self.info_file_obj.close()


def pretty_print(n):
    """Pretty print function for very big integers"""
    if type(n) != int:
        return n

    ret = []
    n = str(n)
    for i in range(len(n) - 1, -1, -1):
        ret.append(n[i])
        if (len(n) - i) % 3 == 0:
            ret.append(',')
    ret.reverse()
    return ''.join(ret[1:]) if ret[0] == ',' else ''.join(ret)


def get_date():
    return time.strftime("%d %b %y %H:%M:%S", time.localtime())


def get_terminal_size():
    """function was taken from http://stackoverflow.com/a/566752"""
    def ioctl_GWINSZ(fd):
        try:
            cr = struct.unpack('hh', fcntl.ioctl(fd, termios.TIOCGWINSZ,
        '1234'))
        except:
            return None
        return cr
    cr = ioctl_GWINSZ(0) or ioctl_GWINSZ(1) or ioctl_GWINSZ(2)
    if not cr:
        try:
            fd = os.open(os.ctermid(), os.O_RDONLY)
            cr = ioctl_GWINSZ(fd)
            os.close(fd)
        except:
            pass
    if not cr:
        try:
            cr = (os.environ['LINES'], os.environ['COLUMNS'])
        except:
            cr = (25, 80)
    return int(cr[1]), int(cr[0])


def get_pretty_name(key):
    if pretty_names.has_key(key):
        return pretty_names[key]
    else:
        return key


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


def is_file_exists(file_path):
    if not file_path:
        raise ConfigError, "No input file is declared..."
    if not os.path.exists(file_path):
        raise ConfigError, "No such file: '%s' :/" % file_path
    return True


def is_output_file_writable(file_path):
    if not file_path:
        raise ConfigError, "No output file is declared..."
    if not os.access(os.path.dirname(ABS(file_path)), os.W_OK):
        raise ConfigError, "You do not have permission to generate the output file '%s'" % file_path


def is_file_tab_delimited(file_path):
    is_file_exists(file_path)
    f = open(file_path)
    line = f.readline()
    if len(line.split('\t')) == 1:
        raise ConfigError, "File '%s' does not seem to have TAB characters.\
                            Did you export this file on MAC using EXCEL? :(" % file_path

    f.seek(0)
    if len(set([len(line.split('\t')) for line in f.readlines()])) != 1:
        raise ConfigError, "Not all lines in the file '%s' have equal number of fields..." % file_path

    f.close()
    return True


def is_file_json_formatted(file_path):
    try:
        json.load(open(file_path))
    except ValueError, e:
        raise ConfigError, "File '%s' does seem to be a properly formatted JSON\
                            file ('%s', cries the library)." % (file_path, e)

    return True


def is_file_fasta_formatted(file_path):
    try:
        f = u.SequenceSource(file_path)
    except u.FastaLibError, e:
        raise ConfigError, "Someone is not happy with your FASTA file '%s' (this is\
                            what the lib says: '%s'." % (file_path, e)

    f.close()

    return True


def is_program_exists(program):
    """adapted from http://stackoverflow.com/a/377028"""
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return True
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return True

    raise ConfigError, "'%s' is not found" % program


def run_command(cmdline):
    try:
        if subprocess.call(cmdline, shell = True) < 0:
            raise ConfigError, "command was terminated by signal: %d" % (-retcode)
    except OSError, e:
        raise ConfigError, "command was failed for the following reason: '%s' ('%s')" % (e, cmdline)


def store_dict_as_TAB_delimited_file(d, output_path, headers):
    is_output_file_writable(output_path)
    f = open(output_path, 'w')
    f.write('%s\n' % '\t'.join(headers))
    for k in d.keys():
        line = [k]
        for header in headers[1:]:
            line.append(d[k][header])
        f.write('%s\n' % '\t'.join(line))
    return output_path


def get_header_fields_of_TAB_delim_file(file_path):
    return open(file_path).readline().strip('\n').split('\t')

def get_json_obj_from_TAB_delim_metadata(input_file):
    return json.dumps([line.strip('\n').split('\t') for line in open(input_file).readlines()])


def get_all_ids_from_fasta(input_file):
    fasta = u.SequenceSource(input_file)
    ids = []
    
    while fasta.next():
        ids.append(fasta.id) 

    return ids


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
        is_file_exists(f)
    is_output_file_writable(dest_file)

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
        return 'Config Error: %s' % textwrap.fill(self.e, 80)


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


def get_temp_directory_path():
    return tempfile.mkdtemp()


def get_temp_file_path():
    f = tempfile.NamedTemporaryFile(delete=False)
    temp_file_name = f.name
    f.close()
    return temp_file_name


def gen_output_directory(output_directory, progress=None, delete_if_exits = False):
    if os.path.exists(output_directory) and delete_if_exits:
        try:
            os.rmtree(output_directory)
        except:
            if progress:
                progress.end()
            raise ConfigError, "I was instructed to remove this directory, but I failed: '%s' :/" % output_directory

    if not os.path.exists(output_directory):
        try:
            os.makedirs(output_directory)
        except:
            if progress:
                progress.end()
            raise ConfigError, "Output directory does not exist (attempt to create one failed as well): '%s'" % \
                                                            (output_directory)
    if not os.access(output_directory, os.W_OK):
        if progress:
            progress.end()
        raise ConfigError, "You do not have write permission for the output directory: '%s'" % output_directory


def check_project_name(project_name):
    allowed_chars = string.ascii_letters + string.digits + '_' + '-' + '.'
    if project_name:
        if len([c for c in project_name if c not in allowed_chars]):
            raise ConfigError, "Project name ('%s') contains characters that PaPi does not like. Please\
                                limit the characters that make up the project name to ASCII letters,\
                                digits, '_' and '-' (if you had not declared a project name and PaPi made\
                                up one for you, please specify with '-p' parameter specifically)." % project_name


def check_contig_names(contig_names):
    allowed_chars = string.ascii_letters + string.digits + '_' + '-' + '.'
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
    is_file_tab_delimited(file_path)

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


def get_newick_tree_data(observation_matrix_path, output_file_name = None, clustering_distance='euclidean', clustering_method = 'complete'):
    vectors = []
    id_to_sample_dict = {}

    is_file_exists(observation_matrix_path)
    is_file_tab_delimited(observation_matrix_path)

    if output_file_name:
        output_file_name = ABS(output_file_name)
        output_directory = os.path.dirname(output_file_name)
        if not os.access(output_directory, os.W_OK):
            raise ConfigError, "You do not have write permission for the output directory: '%s'" % output_directory
    
    input_matrix = open(observation_matrix_path)
    input_matrix.readline()

    line_counter = 0
    for line in input_matrix.readlines():
        fields = line.strip().split('\t')
        id_to_sample_dict[line_counter] = fields[0]
        vector = [float(x) for x in fields[1:]]
        denominator = sum(vector)
        normalized_vector = [p / denominator for p in vector]
        vectors.append(normalized_vector)
        line_counter += 1

    distance_matrix = hcluster.pdist(vectors, clustering_distance)
    
    #clustering_result = hcluster.ward(distance_matrix)
    clustering_result = hcluster.linkage(distance_matrix, method = clustering_method)
    
    tree = hcluster.to_tree(clustering_result)
    
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
   
    if output_file_name:
        root.write(format=1, outfile=output_file_name) 

    return root.write(format=1) 


def reset_output_dir(runinfo_dict_path, cwd=None):
    try:
        runinfo_dict = cPickle.load(open(runinfo_dict_path))
    except:
        raise PaPi.utils.ConfigError, "The input file ('%s') does not seem to be a cPickle object." % (runinfo_dict_path)

    cwd = os.getcwd()

    new_output_dir = os.path.dirname(os.path.join(cwd, runinfo_dict_path))
    old_output_dir = runinfo_dict['output_dir']

    for key in runinfo_dict.keys():
        if isinstance(runinfo_dict[key], str):
            runinfo_dict[key] = runinfo_dict[key].replace(old_output_dir, new_output_dir)

    cPickle.dump(runinfo_dict, open(runinfo_dict_path, 'w'))
    return runinfo_dict


# -*- coding: utf-8
# pylint: disable=line-too-long
"""File/Path operations"""

import os
import h5py
import mmap
import json
import time
import shutil
import tempfile
import anvio.fastalib as u

import anvio

from anvio.terminal import Run
from anvio.terminal import Progress
from anvio.terminal import SuppressAllOutput
from anvio.errors import FilesNPathsError, SamplesError

with SuppressAllOutput():
    from ete2 import Tree

__author__ = "A. Murat Eren"
__copyright__ = "Copyright 2015, The anvio Project"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "A. Murat Eren"
__email__ = "a.murat.eren@gmail.com"
__status__ = "Development"


def is_proper_newick(newick_data):
    try:
        return Tree(newick_data, format=1)
    except Exception as e:
        raise FilesNPathsError, "Your tree doesn't seem to be properly formatted. Here is what ete2 had\
                                 to say about this: '%s'. Pity :/" % e


def is_proper_hdf5_file(hdf5_file_path):
    is_file_exists(hdf5_file_path)

    try:
        h5py.File(hdf5_file_path, 'r')
    except:
        raise FilesNPathsError, "The file '%s' does not seem to be a properly formatted HDF5 data file. Are you sure\
                                 anvi'o generated this?" % (hdf5_file_path)

    return True


def is_proper_genomes_storage_file(storage_path):
    is_file_exists(storage_path)
    is_proper_hdf5_file(storage_path)

    fp = h5py.File(storage_path, 'r')

    if '/data/genomes' not in fp or '/info/genomes' not in fp:
        raise FilesNPathsError, "The file '%s' does not seem to be a proper genomes storage file. If you are not just\
                                 sending random HDF5 files as parameters to mock anvi'o and you are certain that this\
                                 is a geniune genomes storage file, then something may have gone wrong during the\
                                 process that attempted to create it :/" % (storage_path)

    if len(fp['/info/genomes']) != len(fp['/data/genomes']):
        raise FilesNPathsError, "The file '%s' has different number of genomes for data (%d) and information (%d) sections. This\
                                 would have never ever happened if things had gone properly when you generated the file. With its\
                                 current form, there is nothing anvi'o can do with this file :( Sorry about the cryptic and\
                                 not-quite-helpful error message..." % (storage_path, len(fp['/data/genomes']), len(fp['/info/genomes']))


    return True


def is_proper_samples_information_file(file_path):
    is_file_tab_delimited(file_path)

    f = open(file_path, 'rU')

    # quick checks for the header
    columns = f.readline().strip('\n').split('\t')

    if columns[0] != 'samples':
        raise  SamplesError, "The first column of the first row of an anvi'o samples information file\
                              must say 'samples'."

    if len(columns[1:]) != len(set(columns[1:])):
        raise SamplesError, "Every column name in the anvi'o samples information file must be unique (obviously)."

    # quick checks for the samples described
    sample_names = [l.strip('\n').split('\t')[0] for l in f.readlines()]

    if len(sample_names) != len(set(sample_names)):
        raise SamplesError, "Every sample name in the anvi'o samples information file must be unique :/"

    f.close()

    return sample_names


def is_proper_samples_order_file(file_path):
    is_file_tab_delimited(file_path)

    f = open(file_path, 'rU')

    columns = f.readline().strip().split('\t')

    if len(columns) != 3:
        raise SamplesError, "The number of columns in an anvi'o samples order file must be three.\
                             Yours has %d. Please see the documentation if you are lost." % len(columns)

    if columns[0] != 'attributes':
        raise  SamplesError, "The first column of the first row of an anvi'o samples order file \
                              must say 'attributes'. All these rules... Anvi'o promises that they \
                              are for your own good."

    if columns[1] != 'basic':
        raise  SamplesError, "The second column of the first row of an anvi'o samples order file \
                              must read 'basic'."

    if columns[2] != 'newick':
        raise  SamplesError, "The third column of the first row of an anvi'o samples order file \
                              must read 'basic'."

    num_samples_described_in_basic_organizations = []
    num_samples_described_in_newick_organizations = []
    sample_names_described_by_each_organization = []

    for columns in [l.strip('\n').split('\t') for l in f.readlines()]:
        if len(columns) != 3:
            raise SamplesError, "Each line in the samples order file must contain three columns separated\
                                 from each other by TAB characters. You have at least one with %d columns\
                                 :/" % len(columns)

        attribute, basic, newick = columns

        if basic and newick:
            raise SamplesError, 'For the attribute %s, there is both basic and newick form of organization\
                                in the samples order file. For a given attribute, you can define only one\
                                of them, and the other must be blank.' % attribute
        if not basic and not newick:
            raise SamplesError, 'For the attribute %s, there is no organization defined (neither newick, nor\
                                 basic). Is this a test or something? :/' % attribute

        if newick:
            try:
                tree = is_proper_newick(newick)
            except:
                raise SamplesError, 'The newick entry for the attribute %s deos not seem to be a properly\
                                     formatted newick :/' % attribute
            samples = [n.name for n in tree.get_leaves()]
            num_samples_described_in_newick_organizations.append(len(samples))
            sample_names_described_by_each_organization.append(samples)

        if basic:
            if not basic.count(','):
                raise SamplesError, 'The basic samples organization for attribute %s does not seem to be a\
                                     comma-separated list.'
            samples = [s.strip() for s in basic.split(',')]
            num_samples_described_in_basic_organizations.append(len(samples))
            sample_names_described_by_each_organization.append(samples)

    if num_samples_described_in_basic_organizations and len(set(num_samples_described_in_basic_organizations)) != 1:
        raise SamplesError, 'The number of samples described by each comma-separated basic organization line\
                             must be equal. But that does not seem to be the case with your input :/'

    if num_samples_described_in_newick_organizations and len(set(num_samples_described_in_newick_organizations)) != 1:
        raise SamplesError, 'The number of samples described by each newick-formatted organization \
                             must be equal. But that does not seem to be the case with your input :/'

    unique_list_of_samples = [list(x) for x in set(tuple(sorted(x)) for x in sample_names_described_by_each_organization)]
    if len(unique_list_of_samples) != 1:
        raise SamplesError, "At least one organization in the samples order file differs from the others. Each\
                             order should contain the same sample names. Sorry about the cryptic error message,\
                             but your file is not properly formatted :/"

    return unique_list_of_samples[0]


def is_file_exists(file_path, dont_raise=False):
    if not file_path:
        raise FilesNPathsError, "No input file is declared..."
    if not os.path.exists(file_path):
        if dont_raise:
            return False
        else:
            raise FilesNPathsError, "No such file: '%s' :/" % file_path
    return True


def is_output_file_writable(file_path):
    if not file_path:
        raise FilesNPathsError, "No output file is declared..."
    if not os.access(os.path.dirname(os.path.abspath(file_path)), os.W_OK):
        raise FilesNPathsError, "You do not have permission to generate the output file '%s'" % file_path
    return True


def is_output_dir_writable(dir_path):
    if not dir_path:
        raise FilesNPathsError, "No output directory path is declared..."
    if not os.path.isdir(dir_path):
        raise FilesNPathsError, "'%s' is not a directory..." % dir_path
    if not os.access(os.path.abspath(dir_path), os.W_OK):
        raise FilesNPathsError, "You do not have permission to generate files in '%s'" % dir_path
    return True


def is_dir_empty(dir_path):
    if not dir_path:
        raise FilesNPathsError, "is_dir_empty: No directory path is declared..."
    if not os.path.isdir(dir_path):
        raise FilesNPathsError, "is_dir_empty: '%s' is not a directory..." % dir_path

    return False if len(os.listdir(dir_path)) else True


def is_file_tab_delimited(file_path, separator='\t', expected_number_of_fields=None):
    is_file_exists(file_path)
    f = open(file_path, 'rU')

    while True:
        line = f.readline().strip(' ')
        if line.startswith('#'):
            continue
        else:
            break

    if len(line.split(separator)) == 1 and expected_number_of_fields != 1:
        raise FilesNPathsError, "File '%s' does not seem to have TAB characters.\
                            Did you export this file on MAC using EXCEL? :(" % file_path

    f.seek(0)
    num_fields_set = set([len(line.split(separator)) for line in f.readlines()])
    if len(num_fields_set) != 1:
        raise FilesNPathsError, "Not all lines in the file '%s' have equal number of fields..." % file_path

    if expected_number_of_fields:
        num_fields_in_file = list(num_fields_set)[0]
        if num_fields_in_file != expected_number_of_fields:
            raise FilesNPathsError, "The expected number of fileds for '%s' is %d. Yet, it has %d\
                                     of them :/" % (file_path, expected_number_of_fields, num_fields_in_file)

    f.close()
    return True


def is_file_json_formatted(file_path):
    is_file_exists(file_path)

    try:
        json.load(open(file_path, 'rU'))
    except ValueError as e:
        raise FilesNPathsError, "File '%s' does seem to be a properly formatted JSON\
                            file ('%s', cries the library)." % (file_path, e)

    return True


def is_file_fasta_formatted(file_path):
    is_file_exists(file_path)

    try:
        f = u.SequenceSource(file_path)
    except u.FastaLibError as e:
        raise FilesNPathsError, "Someone is not happy with your FASTA file '%s' (this is\
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

    raise FilesNPathsError, "'%s' is not found" % program


def get_temp_directory_path():
    return tempfile.mkdtemp()


def get_temp_file_path():
    f = tempfile.NamedTemporaryFile(delete=False)
    temp_file_name = f.name
    f.close()
    return temp_file_name


def get_num_lines_in_file(file_path):
    if os.stat(file_path).st_size == 0:
        return 0

    f = open(file_path, "rU+")
    buf = mmap.mmap(f.fileno(), 0)
    num_lines = 0
    readline = buf.readline
    while readline():
        num_lines += 1

    return num_lines


def check_output_directory(output_directory, ok_if_exists=False):
    if not output_directory:
        raise FilesNPathsError, "Sorry. You must declare an output directory path."

    output_directory = os.path.abspath(output_directory)

    if os.path.exists(output_directory) and not ok_if_exists:
        raise FilesNPathsError, "The output directory already exists. anvio does not like overwriting stuff."

    return output_directory


def gen_output_directory(output_directory, progress=Progress(verbose=False), run=Run(), delete_if_exists=False):
    if os.path.exists(output_directory) and delete_if_exists and not is_dir_empty(output_directory):
        try:
            run.warning('filesnpaths::gen_output_directory: the client asked the existing directory \
                         "%s" to be removed.. Just so you know :/ (You have 5 seconds to press\
                         CTRL + C).' % output_directory)
            time.sleep(5)
            shutil.rmtree(output_directory)
        except:
            progress.end()
            raise FilesNPathsError, "I was instructed to remove this directory, but I failed: '%s' :/" % output_directory

    if not os.path.exists(output_directory):
        try:
            os.makedirs(output_directory)
        except:
            progress.end()
            raise FilesNPathsError, "Output directory does not exist (attempt to create one failed as well): '%s'" % \
                                                            (output_directory)
    if not os.access(output_directory, os.W_OK):
        progress.end()
        raise FilesNPathsError, "You do not have write permission for the output directory: '%s'" % output_directory

    return output_directory


def get_name_from_file_path(file_path, postfix_separator="."):
    """Return a decent name for a given file at file at `file_path`"""

    splits = os.path.basename(os.path.abspath(file_path)).split(postfix_separator)

    return splits[0] if len(splits) in [1, 2] else '.'.join(splits[:-1])

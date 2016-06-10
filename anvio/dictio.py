# -*- coding: utf-8
# pylint: disable=line-too-long
"""Module to read and write serialized+compressed anvio objects"""

import gzip
import cPickle

from anvio.errors import DictIOError


__author__ = "A. Murat Eren"
__copyright__ = "Copyright 2015, The anvio Project"
__credits__ = []
__license__ = "GPL 3.0"
__maintainer__ = "A. Murat Eren"
__email__ = "a.murat.eren@gmail.com"
__status__ = "Development"


def write_serialized_object(obj, output_file_path):
    """Write serialized object on disk"""
    with gzip.GzipFile(output_file_path, 'w') as output_file:
        cPickle.dump(obj, output_file)


def read_serialized_object(input_file_path):
    """Read serialized object from disk"""
    try:
        with gzip.open(input_file_path, 'rb') as input_file:
            data = input_file.read()
    except IOError:
        raise DictIOError, "anvio is having very hard time reading '%s' as a dictionary. Maybe you\
                            have an idea why?" % input_file_path

    try:
        return cPickle.loads(data)
    except:
        raise DictIOError, "The input file ('%s') does not seem to be a cPickle object." % (input_file_path)


def strip_prefix_from_dict_values(input_dict, prefix):
    """Remove a given prefix from every item in a dict"""
    for key in input_dict.keys():
        if key in ['output_dir', 'input_bam']:
            continue
        if isinstance(input_dict[key], str) and input_dict[key].startswith(prefix):
            input_dict[key] = input_dict[key][len(prefix):].strip('/')

    return input_dict

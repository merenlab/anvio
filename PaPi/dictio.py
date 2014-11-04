# -*- coding: utf-8
#
# Copyright (C) 2014, A. Murat Eren
#
# This program is free software; you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free
# Software Foundation; either version 2 of the License, or (at your option)
# any later version.
#
# Please read the COPYING file.

import gzip
import cPickle
import textwrap

import PaPi.constants as constants

class DictIO(Exception):
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
        return 'Dict Error: %s' % textwrap.fill(self.e, 80)


def write_serialized_object(obj, output_file_path):
    with gzip.GzipFile(output_file_path, 'w') as output_file:
        cPickle.dump(obj, output_file)


def read_serialized_object(input_file_path):
    with gzip.open(input_file_path, 'rb') as input_file:
        data = input_file.read()

    try:
        return cPickle.loads(data)
    except:
        raise DictIO, "The input file ('%s') does not seem to be a cPickle object." % (runinfo_dict_path)


def strip_prefix_from_dict_values(d, prefix):
    for key in d.keys():
        if key in ['output_dir', 'input_bam']:
            continue
        if isinstance(d[key], str) and d[key].startswith(prefix):
            d[key] = d[key][len(prefix):].strip('/')

    return d


def reset_output_dir(runinfo_dict_path, old_output_dir, new_output_dir):
    runinfo_dict = read_serialized_object(runinfo_dict_path)
    runinfo_dict['output_dir'] = new_output_dir
    write_serialized_object(runinfo_dict, runinfo_dict_path)
    return runinfo_dict
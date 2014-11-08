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

# files, and paths

import os
import shutil
import textwrap
import tempfile
import PaPi.fastalib as u

from PaPi.terminal import Progress

ABS = lambda x: x if x.startswith('/') else os.path.join(os.getcwd(), x)

class FilesNPathsError(Exception):
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
        error_type = 'File/Path Error'

        max_len = max([len(l) for l in textwrap.fill(self.e, 80).split('\n')])
        error_lines = ['\033[0;30m\033[46m%s%s\033[0m' % (l, ' ' * (max_len - len(l)))\
                                         for l in textwrap.fill(self.e, 80).split('\n')]

        error_message = ['%s: %s' % (error_type, error_lines[0])]
        for error_line in error_lines[1:]:
            error_message.append('%s%s' % (' ' * (len(error_type) + 2), error_line))

        return '\n'.join(error_message)


def is_file_exists(file_path):
    if not file_path:
        raise FilesNPathsError, "No input file is declared..."
    if not os.path.exists(file_path):
        raise FilesNPathsError, "No such file: '%s' :/" % file_path
    return True


def is_output_file_writable(file_path):
    if not file_path:
        raise FilesNPathsError, "No output file is declared..."
    if not os.access(os.path.dirname(ABS(file_path)), os.W_OK):
        raise FilesNPathsError, "You do not have permission to generate the output file '%s'" % file_path
    return True


def is_file_tab_delimited(file_path):
    is_file_exists(file_path)
    f = open(file_path)
    line = f.readline()
    if len(line.split('\t')) == 1:
        raise FilesNPathsError, "File '%s' does not seem to have TAB characters.\
                            Did you export this file on MAC using EXCEL? :(" % file_path

    f.seek(0)
    if len(set([len(line.split('\t')) for line in f.readlines()])) != 1:
        raise FilesNPathsError, "Not all lines in the file '%s' have equal number of fields..." % file_path

    f.close()
    return True


def is_file_json_formatted(file_path):
    try:
        json.load(open(file_path))
    except ValueError, e:
        raise FilesNPathsError, "File '%s' does seem to be a properly formatted JSON\
                            file ('%s', cries the library)." % (file_path, e)

    return True


def is_file_fasta_formatted(file_path):
    try:
        f = u.SequenceSource(file_path)
    except u.FastaLibError, e:
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


def gen_output_directory(output_directory, progress=Progress(verbose=False), delete_if_exits = False):
    if os.path.exists(output_directory) and delete_if_exits:
        try:
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
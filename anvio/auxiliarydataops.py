# -*- coding: utf-8
# pylint: disable=line-too-long
"""Module to deal with HDF5 files"""

import h5py
import numpy as np

import anvio
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths

from anvio.errors import HDF5Error


__author__ = "A. Murat Eren"
__copyright__ = "Copyright 2015, The anvio Project"
__credits__ = []
__license__ = "GPL 3.0"
__maintainer__ = "A. Murat Eren"
__email__ = "a.murat.eren@gmail.com"


run = terminal.Run()
progress = terminal.Progress()


class HDF5_IO(object):
    def __init__(self, file_path, unique_hash, create_new = False, ignore_hash = False):
        self.file_path = file_path

        if create_new:
            if ignore_hash:
                raise HDF5Error, "When creating a new database, you can't use the 'ignore_hash'\
                                  parameter."

            if not unique_hash:
                raise HDF5Error, "When creating a new database, the 'unique_hash' cannot be None"

            self.fp = h5py.File(self.file_path, 'w')
            self.fp.attrs['hash'] = unique_hash
            self.fp.attrs['version'] = anvio.__hdf5__version__
        else:
            filesnpaths.is_file_exists(self.file_path)
            self.fp = h5py.File(self.file_path, 'r')

            if [h not in self.fp.attrs for h in ['hash', 'version']].count(True):
                raise HDF5Error, "The database at '%s' is missing one or more essential headers that\
                                  should appear in every anvi'o generated HDF5 file. Sorry!" % self.file_path

            if self.fp.attrs['version'] != anvio.__hdf5__version__:
                raise HDF5Error, "The database at '%s' is at version '%s', however your cliend is at\
                                  version '%s'. Bad news." % (self.file_path, self.fp.attrs['version'], anvio.__hdf5__version__)

            if not ignore_hash and self.fp.attrs['hash'] != unique_hash:
                raise HDF5Error, "The database at '%s' does not seem to be compatible with the client :/\
                                  (i.e., the hash values do not match)." % self.file_path


    def add_integer_list(self, path, l, data_type = 'uint16'):
        """Add an array into the the HDF5 file.
        
            >>> h = HDF5_IO('test.h5')
            >>> l = [1, 2, 3, 4, 5]
            >>> h.add_integer_list('/split_1/sample_x', l)
            >>> h.close()
        """

        new_data_obj = self.fp.create_dataset(path, (len(l),), dtype=np.dtype(data_type))
        new_data_obj[...] = np.array(l)


    def get_integer_list(self, path):
        l = self.fp[path]
        return l.value


    def path_exists(self, path):
        return path in self.fp


    def close(self):
        self.fp.close()


class AuxiliaryDataForSplitCoverages(HDF5_IO):
    """A class to handle HDF5 operations to store and access split coverages"""
    def __init__(self, file_path, db_hash, create_new = False, ignore_hash = False, run=run, progress=progress, quiet = False):
        HDF5_IO.__init__(self, file_path, db_hash, create_new = create_new, ignore_hash = ignore_hash)

        self.quiet = quiet


    def is_known_split(self, split_name):
        if not self.path_exists('/data/coverages/%s' % split_name):
            raise HDF5Error, 'The database at "%s" does not know anything about "%s" :(' % (self.file_path, split_name)


    def append(self, split_name, sample_id, coverage_list):
        self.add_integer_list('/data/coverages/%s/%s' % (split_name, sample_id), coverage_list)


    def get(self, split_name):
        self.is_known_split(split_name)

        d = {}

        sample_names = self.fp['/data/coverages/%s' % split_name].keys()

        for sample_name in sample_names:
            d[sample_name] = self.get_integer_list('/data/coverages/%s/%s' % (split_name, sample_name))

        return d


class AuxiliaryDataForNtPositions(HDF5_IO):
    """A class to handle HDF5 operations to store and access split coverages"""
    def __init__(self, file_path, db_hash, create_new = False, run=run, progress=progress, quiet = False):
        HDF5_IO.__init__(self, file_path, db_hash, create_new = create_new)

        self.quiet = quiet


    def is_known_contig(self, contig_name):
        path = '/data/nt_position_info/%s' % contig_name
        return self.path_exists(path)


    def append(self, contig_name, position_info_list):
        self.add_integer_list('/data/nt_position_info/%s' % contig_name, position_info_list, data_type = 'uint8')


    def get(self, contig_name):
        if not self.is_known_contig(contig_name):
            return []

        return self.get_integer_list('/data/nt_position_info/%s' % contig_name)

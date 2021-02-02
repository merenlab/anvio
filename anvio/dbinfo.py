#! /usr/bin/env python

from anvio.db import DB
import anvio.utils as utils
import anvio.filesnpaths as filesnpaths
from anvio.terminal import Run, Progress
from anvio.errors import ConfigError

import os

from abc import ABC, abstractmethod


class DBInfo(ABC):
    db_type = None

    def __new__(cls, path, dont_raise=False):
        if not cls.is_db(path, dont_raise=dont_raise):
            return

        db_type = cls.get_type(path)

        if db_type in ['contigs', 'profile', 'structure']:
            return super().__new__(StructureDBInfo)
        else:
            raise NotImplementedError(f"db_type {db_type} has not been implemented for DBInfo")


    def __init__(self, path):
        if self.db_type is None:
            raise NotImplementedError(f"{self.__class__} must set a `db_type` attribute")

        self.path = path


    @staticmethod
    def is_db(path, dont_raise=False):
        try:
            with DB(path, None, ignore_version=True) as database:
                if 'self' not in database.get_table_names():
                    return False
        except Exception as e:
            if dont_raise:
                return False
            else:
                raise ConfigError(f"Someone downstream doesn't like your so called database, '{path}'. They say "
                                  f"\"{e}\". Awkward :(")
        return True


    @staticmethod
    def get_type(path):
        with DB(path, None, ignore_version=True) as database:
            return database.get_meta_value('db_type')


    @abstractmethod
    def hash(self):
        """This must be defined in the Child class. @abstractmethod takes no crap and ensures this happens"""
        pass


    @property
    def variant(self):
        with self.load_db() as database:
            return database.get_meta_value('db_variant', return_none_if_not_in_table=True)


    def load_db(self):
        return DB(self.path, None, ignore_version=True)



class ContigsDBInfo(DBInfo):
    db_type = 'contigs'
    def __init__(self, path, dont_raise):
        DBInfo.__init__(self, path)


    @property
    def hash(self):
        with self.load_db() as database:
            return database.get_meta_value('contigs_db_hash', return_none_if_not_in_table=True)


class StructureDBInfo(DBInfo):
    db_type = 'structure'
    def __init__(self, path, dont_raise):
        DBInfo.__init__(self, path)


    @property
    def hash(self):
        with self.load_db() as database:
            return database.get_meta_value('contigs_db_hash', return_none_if_not_in_table=True)

class ProfileDBInfo(DBInfo):
    db_type = 'profile'
    def __init__(self, path, dont_raise):
        DBInfo.__init__(self, path)


    @property
    def hash(self):
        with self.load_db() as database:
            return database.get_meta_value('contigs_db_hash', return_none_if_not_in_table=True)


class FindAnvioDBs(object):
    """A helper class to traverse a directory to find anvi'o databases.

    The promary data structure in this class is `self.anvio_dbs` and is initiated
    upon creating an instance from it.

    Parameters
    ==========
    search_path : str
        The beginning of the search. The search will be limited to this directory
        and files and directoreies underneath it.

    max_files_and_dirs_to_process : int, default 50000
        Stop processing if the number of files and directories processed exceeds
        this.
    """

    def __init__(self, search_path='.', max_files_and_dirs_to_process=50000, run=Run(), progress=Progress()):
        self.run = run
        self.progress = progress

        self.search_path = search_path
        self.max_files_and_dirs_to_process = max_files_and_dirs_to_process

        self.anvio_dbs = {}

        for db_path, level in self.walk(depth=3):
            db_info = DBInfo(db_path, dont_raise=True)

            if db_info is not None:
                if db_info.db_type not in self.anvio_dbs:
                    self.anvio_dbs[db_info.db_type] = []

                # Add a cheeky `level` attribute to db_info for no other reason than to order dbs after the walk
                db_info.level = level
                self.anvio_dbs[db_info.db_type].append(db_info)

        # sort by level, so we know what is closest to the search_path root directory
        for db_type in self.anvio_dbs:
            self.anvio_dbs[db_type] = sorted(self.anvio_dbs[db_type], key=lambda d: d.level)


    def listdir(self, path):
        for filename in os.listdir(path):
            yield os.path.join(path, filename)


    def walk(self, depth=None):
        self.progress.new('Searching files and directories')

        total_file_and_directory_names = 0

        if depth and depth == 1:
            filenames = list(self.listdir(self.search_path))

            total_file_and_directory_names += len(filenames)
            self.progress.update(f"processing {total_file_and_directory_names} ...")

            for filename in [f for f in filenames if f.endswith('.db')]:
                yield (filename, 1)
        else:
            top_pathlen = len(self.search_path) + len(os.path.sep)
            for dirpath, dirnames, filenames in os.walk(self.search_path):
                total_file_and_directory_names += (len(filenames) + len(dirnames))

                if total_file_and_directory_names > self.max_files_and_dirs_to_process:
                    self.progress.end()
                    return

                self.progress.update(f"processing {total_file_and_directory_names} ...")

                dirlevel = dirpath[top_pathlen:].count(os.path.sep)
                if depth and dirlevel >= depth - 1:
                    dirnames[:] = []
                else:
                    for filename in [f for f in filenames if f.endswith('.db')]:
                        file_path = os.path.join(dirpath, filename)
                        yield (file_path, file_path.count('/'))

        self.progress.end()



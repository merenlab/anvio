#! /usr/bin/env python

import os
import json

from abc import ABC, abstractmethod

from anvio.db import DB
from anvio.terminal import Run, Progress
from anvio.tables import versions_for_db_types
from anvio.errors import ConfigError


class DBInfo(ABC):
    """Factory class to instantiate the correct DBInfo class

    Parameters
    ==========
    path : str
        The path of a DB

    dont_raise : bool, default False
        If `dont_raise`, if (1) `path` doesn't exist, (2) `path` isn't an anvi'o DB, or (3) `path`
        is an anvi'o DB of the wrong type (specified with `expecting`), then the class instance will be
        None. Otherwise, an error is raised.

    expecting : str OR iterable, default None
        If you are expecting `path` to point to a specific DB type, you can specify the expected DB
        type or types with this parameter. See examples.

    Examples
    ========
    Load up any DB like so:
    >>> from anvio.dbinfo import DBInfo
    >>> anvio_db = DBInfo('PROFILE.db')

    Probe its attributes:
    >>> anvio_db.db_type
    'profile'
    >>> anvio_db.version
    '35'
    >>> print(anvio_db)
    [
      "current_version: 35",
      "db_type: profile",
      "hash: hashdb41e303",
      "hash_name: contigs_db_hash",
      "path: PROFILE.db",
      "variant: None",
      "version: 35"
    ]

    Try and load something that isn't an anvi'o DB:
    >>> DBInfo('nothing')
    Config Error: Someone downstream doesn't like your so called database, 'nothing'. They say "
              File/Path Error: No such file: 'nothing' :/  ". Awkward :(
    >>> x = DBInfo('nothing', dont_raise=True)
    >>> print(x)
    None

    Specify the expected DB type(s)
    >>> DBInfo('PROFILE.db', expecting='contigs')
    Config Error: Was expecting any of the db types ['contigs'], but 'HIMB083/PROFILE.db' has type 'profile'
    >>> DBInfo('PROFILE.db', expecting=['profile', 'pan']).db_type
    'profile'
    """

    db_type = None
    hash_name = None

    def __new__(cls, path, dont_raise=False, expecting=None):
        if not cls.is_db(path, dont_raise=dont_raise):
            return

        db_type = cls.get_type(path)
        if db_type is None:
            if dont_raise:
                return
            else:
                raise ConfigError(f"The database '{path}' has no 'db_type' row in 'self'")

        if expecting:
            if isinstance(expecting, list):
                pass
            elif isinstance(expecting, str):
                expecting = [expecting]
            else:
                raise ConfigError("DBInfo :: `expecting` must be of type list or str")

            for e in expecting:
                if e not in dbinfo_classes:
                    raise ConfigError(f"You are expecting a DB with db_type '{e}', which is not one of the "
                                      f"possible db_types: {list(dbinfo_classes.keys())}")

        if expecting is not None and db_type not in expecting:
            if dont_raise:
                return
            raise ConfigError(f"Was expecting any of the db types {expecting}, but '{path}' has type '{db_type}'")

        if db_type in dbinfo_classes:
            return super().__new__(dbinfo_classes[db_type]['class'])

        if db_type in ['modules', 'trnaseq']:
            # FIXME I don't have the expertise for these two but no point crashing everytime
            # `anvi-interactive` is run without parameters
            return

        raise NotImplementedError(f"db_type {db_type} has no entry in dbinfo_classes")


    def __init__(self, path):
        if self.db_type is None:
            raise NotImplementedError(f"{self.__class__} must set a `db_type` attribute")
        if self.hash_name is None:
            raise NotImplementedError(f"{self.__class__} must set a `hash_name` attribute")

        self.path = path
        self.current_version = versions_for_db_types[self.db_type]


    def __str__(self):
        return json.dumps([f"{attr}: {self.__getattribute__(attr)}" for attr in dir(self) if not callable(getattr(self, attr)) and not attr.startswith("_")], indent=2)


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
            return database.get_meta_value('db_type', return_none_if_not_in_table=True)


    @property
    def variant(self):
        with self.load_db() as database:
            return database.get_meta_value('db_variant', return_none_if_not_in_table=True)


    @property
    def hash(self):
        with self.load_db() as database:
            return database.get_meta_value(self.hash_name, return_none_if_not_in_table=True)


    @property
    def version(self):
        with self.load_db() as database:
            return database.get_meta_value('version', return_none_if_not_in_table=True)


    def load_db(self):
        return DB(self.path, None, ignore_version=True)



class ContigsDBInfo(DBInfo):
    db_type = 'contigs'
    hash_name = 'contigs_db_hash'
    def __init__(self, path, *args, **kwargs):
        DBInfo.__init__(self, path)


class ProfileDBInfo(DBInfo):
    db_type = 'profile'
    hash_name = 'contigs_db_hash'
    def __init__(self, path, *args, **kwargs):
        DBInfo.__init__(self, path)


class GenesDBInfo(DBInfo):
    db_type = 'genes'
    hash_name = 'contigs_db_hash'
    def __init__(self, path, *args, **kwargs):
        DBInfo.__init__(self, path)


class AuxiliaryDBInfo(DBInfo):
    db_type = 'auxiliary data for coverages'
    hash_name = 'contigs_db_hash'
    def __init__(self, path, *args, **kwargs):
        DBInfo.__init__(self, path)


class StructureDBInfo(DBInfo):
    db_type = 'structure'
    hash_name = 'contigs_db_hash'
    def __init__(self, path, *args, **kwargs):
        DBInfo.__init__(self, path)


class GenomeStorageDBInfo(DBInfo):
    db_type = 'genomestorage'
    hash_name = 'hash'
    def __init__(self, path, *args, **kwargs):
        DBInfo.__init__(self, path)


class PanDBInfo(DBInfo):
    db_type = 'pan'
    hash_name = 'genomes_storage_hash'
    def __init__(self, path, *args, **kwargs):
        DBInfo.__init__(self, path)


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


dbinfo_classes = {
    'contigs': {
        'class': ContigsDBInfo,
    },
    'profile': {
        'class': ProfileDBInfo,
    },
    'auxiliary data for coverages': {
        'class': AuxiliaryDBInfo,
    },
    'genes': {
        'class': GenesDBInfo,
    },
    'structure': {
        'class': StructureDBInfo,
    },
    'genomestorage': {
        'class': GenomeStorageDBInfo,
    },
    'pan': {
        'class': PanDBInfo,
    },
}

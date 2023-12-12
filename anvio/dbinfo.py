#! /usr/bin/env python
# -*- coding: utf-8
"""A module of classes to keep track of anvi'o databases

The module includes to major classes: DBInfo, and FindAnvioDBs.

DBInfo class enables seamless access to databae properties for any given
database file. FindAnvioDBs finds databases that are linked to each
other to initiate interactive jobs automatically, even when no db
parameters are provided by the user.
"""

import os
import json

from abc import ABC

import anvio

from anvio.db import DB
from anvio.errors import ConfigError
from anvio.terminal import Run, Progress
from anvio.tables import versions_for_db_types


__author__ = "Developers of anvi'o (see AUTHORS.txt)"
__copyright__ = "Copyleft 2015-2021, the Meren Lab (http://merenlab.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "Evan Kiefl"
__email__ = "kiefl.evan@gmail.com"


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
    functional_annotation_sources_name = None

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
            if len(expecting) == 1:
                raise ConfigError(f"The database at '{path}' is a {db_type} database but you passed it as a '{expecting[0]}' database :/")
            else:
                raise ConfigError(f"The database at '{path}' is a {db_type} database but your parameters claim that it is of type "
                                  f"either {' or '.join(expecting)} :/")

        if db_type in dbinfo_classes:
            # This is the most important line in this method:
            # Return the respective class that __init__ should be called for
            return super().__new__(dbinfo_classes[db_type])

        raise NotImplementedError(f"Database type `{db_type}` at `{path}` has no entry in dbinfo_classes")


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
        if not path:
            raise ConfigError("A low-level function was expecting a database path, but got `None`. A programmer "
                              "needs to look into this :/ Meanwhile, please check your command line parameters. "
                              "Most likely you need to declare a database path, but you do not. You can see the "
                              "entire traceback if you include the flag `--debug` in your command, which may "
                              "help you figure out where did things start going wrong.")

        if not os.path.exists(path):
            raise ConfigError(f"There is nothing at '{path}' :/")

        if os.path.isdir(path):
            raise ConfigError("But this is no file?! This is ah directory! :(")

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


    @property
    def project_name(self):
        with self.load_db() as database:
            return database.get_meta_value('project_name', return_none_if_not_in_table=True)


    def load_db(self):
        return DB(self.path, None, ignore_version=True)


    def get_self_table(self):
        with DB(self.path, None, ignore_version=True) as database:
            return dict(database.get_table_as_list_of_tuples('self'))
    
    def get_functional_annotation_sources(self):
        """If sources are defined in the database, we return them as a list. Otherwise, we return None."""
        if self.functional_annotation_sources_name:
            with self.load_db() as database:
                return database.get_meta_value(self.functional_annotation_sources_name, return_none_if_not_in_table=True).split(',')
        else:
            return None


class ContigsDBInfo(DBInfo):
    db_type = 'contigs'
    hash_name = 'contigs_db_hash'
    functional_annotation_sources_name = 'gene_function_sources'
    def __init__(self, path, *args, **kwargs):
        DBInfo.__init__(self, path)


class ProfileDBInfo(DBInfo):
    db_type = 'profile'
    hash_name = 'contigs_db_hash'
    def __init__(self, path, *args, **kwargs):
        DBInfo.__init__(self, path)


    @property
    def blank(self):
        with self.load_db() as database:
            return True if database.get_meta_value('blank') == 1 else False


    @property
    def merged(self):
        with self.load_db() as database:
            return True if database.get_meta_value('merged') == 1 else False


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
    functional_annotation_sources_name = 'gene_function_sources'
    def __init__(self, path, *args, **kwargs):
        DBInfo.__init__(self, path)


class PanDBInfo(DBInfo):
    db_type = 'pan'
    hash_name = 'genomes_storage_hash'
    def __init__(self, path, *args, **kwargs):
        DBInfo.__init__(self, path)


class TRNADBInfo(DBInfo):
    db_type = 'trnaseq'
    hash_name = 'trnaseq_db_hash'
    def __init__(self, path, *args, **kwargs):
        DBInfo.__init__(self, path)


class ModulesDBInfo(DBInfo):
    db_type = 'modules'
    hash_name = 'hash'
    functional_annotation_sources_name = 'annotation_sources'
    def __init__(self, path, *args, **kwargs):
        DBInfo.__init__(self, path)


class GenomeViewDBInfo(DBInfo):
    db_type = 'genome-view'
    hash_name = 'genome_view_db_hash'
    def __init__(self, path, *args, **kwargs):
        DBInfo.__init__(self, path)


class FindAnvioDBs(object):
    """A helper class to traverse a directory to find anvi'o databases.

    The primary data structure in this class is `self.anvio_dbs` and is initiated
    upon creating an instance from it.

    As of 2021, the primary sole client of this class is `PopulateAnvioDBArgs` in the
    anvio.argsparse module to fill in missing databases into the args object.

    Parameters
    ==========
    search_path : str
        The beginning of the search. The search will be limited to this directory
        and files and directories underneath it.

    max_files_and_dirs_to_process : int, default 50000
        Stop processing if the number of files and directories processed exceeds
        this.
    """

    def __init__(self, search_path='.', max_files_and_dirs_to_process=50000, depth=3, run=Run(), progress=Progress()):
        self.run = run
        self.progress = progress

        self.depth = int(depth)
        self.search_path = search_path
        self.max_files_and_dirs_to_process = max_files_and_dirs_to_process

        self.anvio_dbs = {}

        for db_path, level in self.walk():
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

        self.anvio_dbs_found = True


    def listdir(self, path):
        for filename in os.listdir(path):
            yield os.path.join(path, filename)


    def walk(self):
        self.progress.new('Searching files and directories')

        total_file_and_directory_names = 0

        if self.depth and self.depth == 1:
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
                if self.depth and dirlevel >= self.depth - 1:
                    dirnames[:] = []
                else:
                    for filename in [f for f in filenames if f.endswith('.db')]:
                        file_path = os.path.join(dirpath, filename)
                        yield (file_path, file_path.count('/'))

        self.progress.end()


dbinfo_classes = {
    'contigs': ContigsDBInfo,
    'profile': ProfileDBInfo,
    'auxiliary data for coverages': AuxiliaryDBInfo,
    'genes': GenesDBInfo,
    'structure': StructureDBInfo,
    'genomestorage': GenomeStorageDBInfo,
    'pan': PanDBInfo,
    'trnaseq': TRNADBInfo,
    'modules': ModulesDBInfo,
    'genome-view': GenomeViewDBInfo,
}

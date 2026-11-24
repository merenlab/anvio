#!/usr/bin/env python

import os
import sys

import anvio
import anvio.db as db
import anvio.tables as t
import anvio.dbinfo as dbi
import anvio.utils as utils
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths

from anvio.argparse import ArgumentParser
from anvio.errors import ConfigError, FilesNPathsError


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['FlorianTrigodet', 'meren']
__requires__ = ['profile-db']
__can_use__ = ['contigs-db', 'pan-db', 'genes-db']
__provides__ = []
__description__ = ("Create or drop a SQLite index on a column of an anvi'o database table. Indexes make "
                   "`WHERE column = ...` lookups fast on very large databases at the cost of disk space and "
                   "build time, so anvi'o does not create them by default -- you opt into them here")


class TableIndexer:
    """Create or drop an index on a (table, column) of any anvi'o database.

    Anvi'o ships without column indexes because most databases are small enough that an index
    would only add disk space. For the rare large database (e.g. a merged profile.db with
    hundreds of millions of SNV rows), a single index turns a tens-of-minutes full table scan
    into a subsecond lookup. This program lets a user opt into one of the curated indexes
    listed in `anvio.tables.indexable_table_columns`, and reverse it later with `--drop-index`.
    """

    def __init__(self, args, run=terminal.Run(), progress=terminal.Progress()):
        self.run = run
        self.progress = progress

        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        self.db_path = A('db_path')
        self.table = A('table')
        columns = A('column')
        self.columns = [c.strip() for c in columns.split(',')] if columns else None
        self.drop_index = A('drop_index')
        self.list_indexable = A('list')
        self.just_do_it = A('just_do_it')
        self.reclaim_space = A('reclaim_space')


    def process(self):
        if self.reclaim_space and not self.drop_index:
            raise ConfigError("The `--reclaim-space` flag only makes sense together with `--drop-index`: it runs a "
                              "VACUUM to shrink the database file after an index is removed. There is nothing to "
                              "reclaim when you are creating an index, so anvi'o is not sure what you intended here.")

        if self.list_indexable:
            self.list_indexable_columns()
            return

        if not self.db_path:
            raise ConfigError("You need to provide an anvi'o database to index (or use `--list` to see "
                              "which tables and columns are indexable).")

        if not (self.table and self.columns):
            raise ConfigError("You must provide both a `--table` and a `--column` to index (or use `--list` "
                              "to see your options).")

        filesnpaths.is_file_exists(self.db_path)
        self.db_type = dbi.DBInfo(self.db_path).db_type

        database = db.DB(self.db_path, None, ignore_version=True)
        try:
            self.sanity_check(database)

            if self.drop_index:
                self.drop(database)
            else:
                self.create(database)
        finally:
            # `progress.end()` is a no-op if no progress is active, so this safely
            # tidies the terminal even if `_exec` threw mid-build/mid-drop.
            self.progress.end()
            database.disconnect()


    def list_indexable_columns(self):
        """Print the curated whitelist of indexable (table, column) combinations.

        If a database was provided, scope the listing to that database's type and mark which
        indexes already exist; otherwise show everything grouped by db type.
        """
        db_type = None
        existing = set()
        if self.db_path:
            filesnpaths.is_file_exists(self.db_path)
            db_type = dbi.DBInfo(self.db_path).db_type
            database = db.DB(self.db_path, None, ignore_version=True)
            for table, columns in t.indexable_table_columns.get(db_type, []):
                if self.index_exists(database, t.index_name_for(table, columns)):
                    existing.add((table, tuple(columns)))
            database.disconnect()

        self.run.warning(None, header="INDEXABLE TABLES AND COLUMNS", lc="green")
        for this_db_type, entries in t.indexable_table_columns.items():
            if db_type and this_db_type != db_type:
                continue
            self.run.info_single(f"{this_db_type}-db", mc='green', nl_before=1)
            for table, columns in entries:
                tag = "  [already indexed]" if (table, tuple(columns)) in existing else ""
                self.run.info_single(f"--table {table} --column {','.join(columns)}{tag}", level=2)

        self.run.info_single("To build one of the above, rerun this program with the corresponding `--table` "
                             "and `--column`. To index something not on this list, add `--just-do-it`.",
                             nl_before=1, nl_after=1, mc='yellow')


    def sanity_check(self, database):
        # raises a helpful ConfigError if the table is not in the db
        database.is_table_exists(self.table)

        table_columns = database.get_table_structure(self.table)
        missing = [c for c in self.columns if c not in table_columns]
        if missing:
            raise ConfigError(f"The table '{self.table}' has no column(s) named {', '.join(missing)}. "
                              f"The columns it does have are: {', '.join(table_columns)}.")

        whitelisted = (self.table, self.columns) in [(tbl, cols) for tbl, cols in t.indexable_table_columns.get(self.db_type, [])]
        if not whitelisted and not self.drop_index and not self.just_do_it:
            valid = t.indexable_table_columns.get(self.db_type, [])
            valid_str = '; '.join(f"--table {tbl} --column {','.join(cols)}" for tbl, cols in valid) or "(none)"
            raise ConfigError(f"Indexing '{self.table}({','.join(self.columns)})' is not on anvi'o's curated list of "
                              f"useful indexes for a {self.db_type}-db, so anvi'o is not going to do it unless you "
                              f"insist. The combinations anvi'o knows are worth indexing for this db type are: "
                              f"{valid_str}. If you really know that indexing this column will help your use case, "
                              f"rerun with `--just-do-it` and anvi'o will build it anyway. You can always reverse it "
                              f"with `--drop-index`.")


    def create(self, database):
        index_name = t.index_name_for(self.table, self.columns)

        if self.index_exists(database, index_name):
            self.run.info_single(f"An index named '{index_name}' already exists on "
                                 f"{self.table}({','.join(self.columns)}). Nothing to do here.", mc='green',
                                 nl_before=1, nl_after=1)
            return

        size_before = os.path.getsize(self.db_path)

        # Note: we deliberately do NOT report an exact row count here. A COUNT(*) on the
        # unindexed table is itself a full scan -- on an 800M-row table that is its own
        # multi-minute cost, paid just to print a number. The database file size below is a
        # cheap, equally useful proxy for "how big is this going to be".
        self.run.warning(f"Anvi'o is about to build an index on '{self.table}({','.join(self.columns)})'. For "
                         f"typical databases this finishes in seconds. For very large tables (hundreds of millions "
                         f"of rows) it can take many minutes and grow the database file by several gigabytes. The "
                         f"current size of this database is {utils.human_readable_file_size(size_before)}. You can "
                         f"reverse this any time by rerunning this program with the `--drop-index` flag.",
                         header="ABOUT TO BUILD AN INDEX -- THIS MAY TAKE A WHILE", lc='yellow')

        self.progress.new("Indexing")
        self.progress.update(f"Building {index_name} ...")
        column_list = ','.join(f'"{c}"' for c in self.columns)
        database._exec(f'CREATE INDEX IF NOT EXISTS "{index_name}" ON "{self.table}"({column_list})')
        self.progress.end()

        size_after = os.path.getsize(self.db_path)
        self.run.info("Index created", index_name, mc='green', nl_before=1)
        self.run.info("Table", self.table)
        self.run.info("Column(s)", ', '.join(self.columns))
        self.run.info("DB size before", utils.human_readable_file_size(size_before))
        self.run.info("DB size after", utils.human_readable_file_size(size_after), nl_after=1)


    def drop(self, database):
        index_name = t.index_name_for(self.table, self.columns)
        index_present = self.index_exists(database, index_name)

        # A missing index is only "nothing to do" when the user also did not ask to reclaim
        # space. If they did, they are most likely re-running after an earlier `--drop-index`
        # (which removes the index but leaves the file size untouched) precisely to vacuum now,
        # so we must NOT bail here -- otherwise the documented "reclaim later" workflow is a
        # dead end that never reaches the VACUUM below.
        if not index_present and not self.reclaim_space:
            self.run.info_single(f"There is no index named '{index_name}' on {self.table}({','.join(self.columns)}) "
                                 f"to drop. Nothing to do here.", mc='green', nl_before=1, nl_after=1)
            return

        size_before = os.path.getsize(self.db_path)

        if index_present:
            database._exec(f'DROP INDEX IF EXISTS "{index_name}"')
        else:
            # Index already gone (e.g. dropped by an earlier run without `--reclaim-space`). The
            # user is re-running just to shrink the file, which is a valid thing to want.
            self.run.info_single(f"There is no index named '{index_name}' to drop -- it looks like it was already "
                                 f"removed. Anvi'o will go ahead and reclaim disk space, as you asked.",
                                 mc='yellow', nl_before=1)

        # Dropping the index frees its pages inside the file, but SQLite does not shrink the file
        # on disk until a VACUUM. VACUUM is NOT free: it rewrites the entire database into a new
        # file, so it needs temporary free space roughly equal to the current database size and,
        # on the hundreds-of-GB profile dbs this feature exists for, can run for a very long time
        # (and fail outright if the filesystem is tight). That is worse than leaving some slack
        # space behind, so anvi'o does not VACUUM by default -- the user opts in with
        # `--reclaim-space`, after the same kind of heads-up that the create path gives.
        if self.reclaim_space:
            self.run.warning(f"You asked anvi'o to reclaim the disk space the dropped index was using. To do that "
                             f"it will run a VACUUM, which rewrites the *entire* database file. This needs free disk "
                             f"space roughly equal to the current size of the database "
                             f"({utils.human_readable_file_size(size_before)}), and on very large databases it can "
                             f"take a long while. If you would rather skip it, interrupt now -- the index is already "
                             f"dropped; the file will simply shrink the next time the database is vacuumed.",
                             header="ABOUT TO VACUUM -- THIS MAY TAKE A WHILE", lc='yellow')
            self.progress.new("Reclaiming space")
            self.progress.update("Running VACUUM ...")
            database._exec('VACUUM')
            self.progress.end()

        size_after = os.path.getsize(self.db_path)

        if index_present:
            self.run.info("Index dropped", index_name, mc='green', nl_before=1)
        self.run.info("DB size before", utils.human_readable_file_size(size_before))
        self.run.info("DB size after", utils.human_readable_file_size(size_after), nl_after=1)
        if not self.reclaim_space:
            self.run.info_single("The index is gone, but the database file has not shrunk yet: SQLite frees the "
                                 "space inside the file and reuses it later. To physically reclaim it now, rerun "
                                 "this command with `--reclaim-space` (which runs a VACUUM).", mc='yellow', nl_after=1)


    def index_exists(self, database, index_name):
        response = database._exec("SELECT name FROM sqlite_master WHERE type='index' AND name=?", (index_name,))
        return len(response.fetchall()) > 0


def get_args():
    parser = ArgumentParser(description=__description__)

    parser.add_argument('db_path', metavar='DATABASE_PATH', nargs='?', default=None,
                        help="An anvi'o database (contigs-db, profile-db, pan-db, or genes-db). Not required "
                             "if you are only using `--list`.")

    groupA = parser.add_argument_group('Table & column', "The table and column you wish to (un)index.")
    groupA.add_argument('--table', default=None, metavar='TABLE_NAME',
                        help="The table to index. See `--list` for the tables anvi'o knows are worth indexing.")
    groupA.add_argument('--column', default=None, metavar='COLUMN_NAME',
                        help="The column to index. For a compound index, provide a comma-separated list of columns.")

    groupB = parser.add_argument_group('Mode', "What to do.")
    groupB.add_argument('--list', default=False, action='store_true',
                        help="List the (table, column) combinations anvi'o knows are worth indexing and exit. "
                             "If you also provide a database, the list is scoped to its type and marks the "
                             "indexes that already exist.")
    groupB.add_argument('--drop-index', default=False, action='store_true',
                        help="Drop the index instead of creating it. By default this only removes the index; the "
                             "database file keeps its current size until it is vacuumed. Add `--reclaim-space` to "
                             "physically shrink the file right away.")
    groupB.add_argument('--reclaim-space', default=False, action='store_true',
                        help="Only meaningful together with `--drop-index`: after dropping the index, run a VACUUM "
                             "to shrink the database file on disk. VACUUM rewrites the entire file, so it needs "
                             "free space roughly equal to the database size and can take a long time on very large "
                             "databases. Without this flag, the freed space is reused by the database later.")
    groupB.add_argument(*anvio.A('just-do-it'), **anvio.K('just-do-it', {'help':
                        "Build an index even if the (table, column) is not on anvi'o's curated list."}))

    return parser.get_args(parser)


@terminal.time_program
def main():
    args = get_args()

    try:
        TableIndexer(args).process()
    except ConfigError as e:
        print(e)
        sys.exit(-1)
    except FilesNPathsError as e:
        print(e)
        sys.exit(-2)


if __name__ == "__main__":
    main()

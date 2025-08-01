#!/usr/bin/env python
# -*- coding: utf-8

import sys
from anvio.argparse import ArgumentParser

import anvio
import anvio.db as db
import anvio.tables as t
import anvio.dbinfo as dbi
import anvio.utils as utils
import anvio.dbops as dbops
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths

from anvio.errors import ConfigError, FilesNPathsError


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['meren', 'ekiefl']
__requires__ = ["pan-db", "profile-db", "contigs-db", "genomes-storage-db", "structure-db", "genes-db"]
__description__ = "Access self tables, display values, or set new ones totally on your own risk"


# FIXME: The functionality of this class has a lot of overlaps with the newer `DBInfo` class
#        defined in `anvio/dbinfo.py`. Printing out database information, for instance, should
#        be carried over there, and the modification of self key/value pairs should be carried
#        into another class ALSO in `anvio/dbinfo.py`.

class DatabaseInfo:
    """A class to get information about an anvi'o database, and edit some of them"""
    def __init__(self, args, quiet=False, run=terminal.Run(), progress=terminal.Progress()):
        self.args = args
        self.run = run
        self.progress = progress

        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        self.db_path = A("input")
        self.self_key = A("self_key")
        self.self_value = A("self_value")
        self.just_do_it = A("just_do_it")

        filesnpaths.is_file_exists(self.db_path)

        self.get_db_self()

        if self.self_key or self.self_value:
            if (not self.self_key) or (not self.self_value):
                raise ConfigError("You must use `--self-key` and `--self-value` together. You can't "
                                  "just use one of them.")

            self.set_key_value()
            self.get_db_self()

        if not quiet:
            run.warning(None, "DB Info (no touch)", lc="green")
            run.info("Database Path", self.db_path)
            if "description" in self.db_self:
                desc = self.db_self["description"]
                if len(desc) > 30:
                    run.info("description", "[Found, and it is %d characters long]" % len(desc))
                elif not desc:
                    run.info("description", "[Not found, but it's OK]", mc="red")
                else:
                    run.info("Description", desc)
            run.info("db_type", f"{self.db_type} (variant: {self.db_variant})")
            run.info("version", self.db_version, nl_after=1)

            run.warning(None, "DB Info (no touch also)", lc="green")
            for key in [k for k in self.db_self if k not in ["version", "db_type", "description", "db_variant"]]:
                if key == self.self_key or (self.self_key == 'sample_id' and key == 'samples'):
                    run.info('%s [JUST SET]' % key, self.db_self[key], mc='green', lc='green')
                else:
                    run.info(key, self.db_self[key])

            run.info_single("Please remember that it is never a good idea to change these values. But in "
                            "some cases it may be absolutely necessary to update something here, and a "
                            "programmer may ask you to run this program and do it. But even then, you "
                            "should be extremely careful.", mc='red', nl_before=1, nl_after=1)

        # It wouldn't hurt to print out some DB specific information?
        if self.db_type == 'contigs':
            # try-excepting here since this program is supposed to be able to run even on
            # broken databases. so if the following code cause issues because the DB is broken
            # this program shouldn't crash.
            try:
                dbops.ContigsDatabase(self.db_path).list_gene_caller_sources()
            except:
                pass

            try:
                dbops.ContigsDatabase(self.db_path).list_function_sources()
            except:
                pass

            try:
                dbops.ContigsDatabase(self.db_path).list_available_hmm_sources()
            except:
                pass
        elif self.db_type == 'genomestorage':
            database = dbi.DBInfo(self.db_path)
            genome_names = database.load_db().get_single_column_from_table(t.genome_info_table_name, 'genome_name')

            run.warning(None, header=f"GENOMES IN STORAGE (n={len(genome_names)})", lc="green")
            run.info_single(f"{', '.join(genome_names)}", cut_after=120, level=0)
        else:
            # no additional information for you.
            pass


    def get_db_self(self):
        database = dbi.DBInfo(self.db_path)

        self.db_type = database.db_type
        self.db_self = database.get_self_table()
        self.db_version = database.version
        self.db_variant = database.variant


    def set_key_value(self):
        filesnpaths.is_output_file_writable(self.db_path)

        if self.self_key in ["version", "db_type", "description", "samples"]:
            raise ConfigError(f"Well. This is awkward. The key {self.self_key} is one of the 'no touch keys', so "
                              f"you can't change them with this program. That said, you are an adult and "
                              f"anvi'o can't tell you what to do. If you want to be a rebel, feel free to "
                              f"use the sqlite from the command line to create your own Frankenstein.")

        # make sure we did everything we can to avoid user error
        self.self_key = self.self_key.replace(' ', '_').strip()

        if self.self_key in self.db_self and self.self_key == 'sample_id':
            if self.db_type != 'profile' or int(self.db_self['blank']) or int(self.db_self['merged']):
                raise ConfigError("You can change the self-key `sample_id` only in single and non-blank anvi'o "
                                  "profile databases :/")

            if not filesnpaths.is_file_exists(dbops.get_auxiliary_data_path_for_profile_db(self.db_path), dont_raise=True):
                raise ConfigError("You must have the auxiliary database in its proper location if you wish to "
                                  "change the key `sample_id` in a single profile database. If you are tired of "
                                  "anvi'o telling you what you can or cannot do, feel free to tell it to go fly a kite, "
                                  "and open an sqlite terminal to change anything the way you want.")

            utils.check_sample_id(self.self_value)

            self.run.warning("You are trying to change a special key, `sample_id`. This is a special one,\
                              because changing it will require a bunch of other things to change as well\
                              both in other tables in this profile database, and in the auxiliary database\
                              associated with it.")

        if self.self_key in self.db_self:
            self.run.warning(f"The content of the key '{self.self_key}' in the self table of your {self.db_type}database, "
                             f"which currently has a value of '{self.db_self[self.self_key]}', will be replaced with '{self.self_value}'.")
        else:
            self.run.warning(f"You are about to set a new key, '{self.self_key}', in the self table of your {self.db_type} "
                             f"database which will hold the value '{self.self_value}'.")

        if not self.just_do_it:
            self.run.info_single('If you like what you see up there, please use the `--just-do-it` flag to make it happen.')
            sys.exit(0)

        database = db.DB(self.db_path, None, ignore_version=True)
        database.remove_meta_key_value_pair(self.self_key)
        database.set_meta_value(self.self_key, self.self_value)
        database.disconnect()

        if self.self_key == 'sample_id':
            database = db.DB(self.db_path, None, ignore_version=True)
            current_sample_name = database.get_meta_value("samples")

            self.progress.new("Changing '%s' to '%s'" % (current_sample_name, self.self_value), progress_total_items=5)
            self.progress.increment()
            self.progress.update('In self ...')
            database.remove_meta_key_value_pair('samples')
            database.set_meta_value('samples', self.self_value)

            self.progress.increment()
            self.progress.update('In layer additional table')
            database._exec("""UPDATE %s SET "item_name" = ? WHERE "item_name" == ?""" % t.layer_additional_data_table_name, (self.self_value, current_sample_name))

            self.progress.increment()
            self.progress.update('In variable nucleotides table')
            database._exec("""UPDATE %s SET "sample_id" = ? WHERE "sample_id" == ?""" % t.variable_nts_table_name, (self.self_value, current_sample_name))

            self.progress.increment()
            self.progress.update('In variable codons table')
            database._exec("""UPDATE %s SET "sample_id" = ? WHERE "sample_id" == ?""" % t.variable_codons_table_name, (self.self_value, current_sample_name))
            database.disconnect()

            self.progress.increment()
            self.progress.update('In auxiliary data table')
            database = db.DB(dbops.get_auxiliary_data_path_for_profile_db(self.db_path), None, ignore_version=True)
            database._exec("""UPDATE %s SET "sample_name" = ? WHERE "sample_name" == ?""" % t.split_coverages_table_name, (self.self_value, current_sample_name))
            database.disconnect()

            self.progress.end()


def get_args():
    parser = ArgumentParser(description=__description__)

    groupA = parser.add_argument_group('Input', "The database path you wish to access.")
    groupA.add_argument("input", metavar = "DATABASE_PATH", nargs=None, help = "An anvi'o database for pan, profile, contigs, or auxiliary data")

    groupB = parser.add_argument_group('Very dangerous zone', "For power users with extreme self-control and maturity.")
    groupB.add_argument(*anvio.A('self-key'), **anvio.K('self-key'))
    groupB.add_argument(*anvio.A('self-value'), **anvio.K('self-value'))
    groupB.add_argument(*anvio.A('just-do-it'), **anvio.K('just-do-it'))

    return  parser.parse_known_args()


@terminal.time_program
def main():
    args, unknown = get_args()

    try:
        DatabaseInfo(args)
    except ConfigError as e:
        print(e)
        sys.exit(-1)
    except FilesNPathsError as e:
        print(e)
        sys.exit(-2)


if __name__ == "__main__":
    main()


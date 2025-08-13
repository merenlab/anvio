#!/usr/bin/env python
# -*- coding: utf-8

import os
import sys
import copy
from anvio.argparse import ArgumentParser
import shutil

import anvio
import anvio.db as db
import anvio.utils as utils
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths

from anvio.migrations import migration_scripts
from anvio.version import versions_for_db_types
from anvio.errors import ConfigError, FilesNPathsError


__description__ = "Migrates any anvi'o artifact, whether it is a database or a config file, to a newer version. Pure magic."
__authors__ = ['meren', 'ozcan', 'ekiefl', 'ivagljiva', 'semiller10', 'mschecht']
__requires__ = ["contigs-db", "profile-db", "pan-db", "genes-db", "genomes-storage-db", "structure-db", "modules-db", "workflow-config"]
__provides__ = []


class Migrater(object):
    def __init__(self, args):
        self.args = args
        self.artifact_path = args.artifact_path
        self.artifact_type = 'Unknown'
        self.artifact_version = None
        self.target_version = None

        self.run = terminal.Run()
        self.progress = terminal.Progress()

        filesnpaths.is_file_exists(self.artifact_path)

        if not self.args.migrate_safely and not self.args.migrate_quickly:
            raise ConfigError("You must choose either `--migrate-safely` or `--migrate-quickly`. Anvi'o would "
                              "have chosen one for you, but our lawyers suggested that it may put our life savings in "
                              "jeopardy.")

        if self.args.migrate_safely and self.args.migrate_quickly:
            raise ConfigError("Please don't ask anvi'o to migrate your files both safely and quickly :/")

        self.safe_mode = self.args.migrate_safely

        try:
            self.sqlite_version = utils.get_command_output_from_shell("sqlite3 --version")[0].decode("utf-8").split(' ')[0].split('.')
        except:
            self.sqlite_version = None

        if not (int(self.sqlite_version[0]) >= 3 and int(self.sqlite_version[1]) >= 30):
            if not self.args.just_do_it:
                raise ConfigError("Anvi'o migration tool requires the installed SQLite version to be at least "
                                  "v3.30.0. Yours seem to be %s (if this doesn't make any sense, you can take "
                                  "a look at the output of the command 'sqlite3 --version' on your terminal). "
                                  "You can skip this check if you use the argument `--just-do-it` with your "
                                  "command. In which case we strongly recommend you to use `--migrate-safely` "
                                  "flag since it will be likely for you to lose your file if something goes "
                                  "wrong downstream." % '.'.join(self.sqlite_version) if self.sqlite_version else 'Unknown :/')
            else:
                self.run.warning("Anvi'o is skipping the version number check for SQLite. Brace yourself for impact?")


    def print_info(self):
        self.run.info('Input file path', self.artifact_path)
        self.run.info('Input file type', self.artifact_type)
        self.run.info('Current Version', self.artifact_version)
        self.run.info('Target Version', self.target_version)

        if self.args.migrate_safely:
           self.run.info('Migration mode', "Safe", mc="green")
        else:
           self.run.info('Migration mode', "Adventurous", mc="red")

        self.run.info('SQLite Version', '.'.join(self.sqlite_version) if self.sqlite_version else 'Unknown :/', nl_after=1, nl_before=1)


    def is_good_to_go(self):
        """Basic checks to make sure file types and so on looks correct"""
        # let's start with very basics
        if os.path.isdir(self.artifact_path):
            self.run.info('Input file path', self.artifact_path, mc='red')
            self.run.warning("This is a directory, and yet the training of the members of the migration headquarters "
                        "of anvi'o only covers files. MOVING ON.", header='NOPE: LA DIRECTORY', lc='yellow')

        # check file extensions to make sure we're clear with this.
        if not self.args.just_do_it and self.artifact_path.split('.')[-1] not in ['db', 'h5', 'json']:
            self.run.info('Input file path', self.artifact_path)
            self.run.info('Input file type', self.artifact_type)
            self.run.warning(f"This program only works with files that end with `.db` or `.json` extensions. "
                        f"But if you are sure that your input file at `{self.artifact_path}` is indeed an anvi'o "
                        f"database or config file you can use the flag `--just-do-it` at your own risk.",
                        header="NOPE: LA FILE EXTENSION", lc='yellow')

            return False

        # make sure we have the Python module for HDF is present if we have an .h5 file
        if self.artifact_path.endswith('GENOMES.h5'):
            utils.check_h5py_module()

        # if we have a JSON file, we have to make sure a few things to say we're good to go
        if self.artifact_path.endswith('json') or self.artifact_path.endswith('JSON'):
            # does it look like a JSON file?
            filesnpaths.is_file_json_formatted(self.artifact_path)

            # does it look like an anvi'o config file?
            import json
            with open(self.artifact_path) as json_file:
                json_content = json.load(json_file)
                if 'workflow_name' not in json_content:
                    self.run.info('Input file path', self.artifact_path)
                    self.run.info('Input file type', self.artifact_type, mc='red')

                    self.run.info_single("Your input file is a proper JSON file, but it doesn't look like an anvi'o workflow "
                                    "since nowhere in it there is a variable called `workflow_name`. So, there is nothing "
                                    "anvi'o can do with this.", nl_before=1)

                    return False

        return True


    def get_artifact_meta(self):
        if not self.is_good_to_go():
            return

        if not self.args.just_do_it and self.artifact_path.split('.')[-1] not in ['db', 'h5', 'json']:
            raise ConfigError("This program only works with files that end with `.db` or `.json` extensions. "
                              "But if you sure that, this is in fact an anvi'o database or config file "
                              "you can use --just-do-it flag at your own risk.")

        file_is_config = False
        if self.args.just_do_it:
            # check if this a config file by checking if it is a JSON formatted file
            try:
                filesnpaths.is_file_json_formatted(self.artifact_path)
                file_is_config = True
            except:
                pass


        try:
            if self.artifact_path.endswith('GENOMES.h5'):
                import h5py
                fp = h5py.File(self.artifact_path, 'r')

                self.artifact_type = 'genomestorage'
                self.artifact_version = int(fp.attrs['version'])

                fp.close()
            elif self.artifact_path.endswith('json') or file_is_config:
                # the file is a json file, but is it a config file for an anvi'o workflow?
                # here we will test it by looking for `workflow_name` in the JSON data.
                import anvio.workflows as w
                self.artifact_type = 'config'
                anvio.QUIET = True
                workflow_name, config_version = w.get_workflow_name_and_version_from_config(self.artifact_path, dont_raise=True)
                self.artifact_version = int(config_version)
                anvio.QUIET = False
            else:
                db_conn = db.DB(self.artifact_path, None, ignore_version=True)

                self.artifact_type = db_conn.get_meta_value('db_type')
                self.artifact_version = int(db_conn.get_meta_value('version'))

                db_conn.disconnect()
        except Exception as e:
            raise ConfigError(f"Anvi'o had a problem with `{self.artifact_path}` :( Here is the error: '{e}'")

        # now we know about the target anvi'o artifact, we can fill in the version.
        self.get_target_version()


    def get_target_version(self):
        if self.artifact_type in versions_for_db_types:
            version = int(versions_for_db_types[self.artifact_type])
        else:
            raise ConfigError(f"Anvi'o does not have any version information about the type '{self.artifact_type}'.")

        if self.args.target_version:
            target = int(self.args.target_version)

            if target <= self.artifact_version:
                raise ConfigError("Target version ('%s') can not be lower than the artifact version ('%s')." % (target, self.artifact_version))
            elif target > version:
                raise ConfigError("Target version ('%s') can not be higher than highest available version ('%s') for this type." % (target, version))

            version = target

        self.target_version = int(version)


    def is_something_to_migrate(self):
        """Test if there is anything worth migrating."""

        if not self.artifact_version:
            self.get_artifact_meta()

            if not self.artifact_version():
                raise ConfigError("Anvi'o needs an adult :( PS: if you are the programmer, you are the adult, and your code "
                                  "should have never made it here.")

        if self.target_version == self.artifact_version:
            self.run.info('Input file path', self.artifact_path)
            self.run.info('Input file type', self.artifact_type)

            self.run.info_single(f"Your {self.artifact_type} looks good and needs no updatin' ☘️", nl_before=1, mc="green")

            return False

        return True


    def process(self):
        self.run.warning(None, header="NEW MIGRATION TASK", lc='cyan')
        if not self.is_good_to_go():
            return

        self.get_artifact_meta()

        if not self.is_something_to_migrate():
            return

        self.print_info()

        tasks = []

        # if we are in safe mode, we will first copy the artifact file
        if self.safe_mode:
            self.progress.new("Migration preparation")
            self.progress.update("Copying the original file for safety ...")
            artifact_file_name = os.path.basename(self.artifact_path)
            temp_file_path = filesnpaths.get_temp_file_path() + artifact_file_name
            shutil.copy(self.artifact_path, temp_file_path)
            self.progress.end()
        else:
            temp_file_path = None

        for i in range(self.artifact_version, self.target_version):
            script_name = "v%s_to_v%s" % (i, i + 1)

            if not self.artifact_type in migration_scripts or not script_name in migration_scripts[self.artifact_type]:
                raise ConfigError("Anvi'o can not find a migrate script required for this operation. (Artifact Type: %s, Script name: %s) " % (self.artifact_type, script_name))

            tasks.append(script_name)

        for script_name in tasks:
            try:
                migration_scripts[self.artifact_type][script_name].migrate(self.artifact_path)
            except Exception as e:
                if temp_file_path:
                    shutil.move(self.artifact_path, self.artifact_path + '.broken')
                    shutil.copy(temp_file_path, self.artifact_path)
                    os.remove(temp_file_path)

                    self.progress.reset()
                    raise ConfigError("Something went wrong during the migration of your file :( But anvi'o was able to restore your "
                                      "original artifact, and stored the unfinished upgrade in the same directory with the '.broken' "
                                      "in case you may need it for debugging purposes. Please feel free to get in touch with the "
                                      "developers, who may ask you to make available the original file to reproduce the problem. This "
                                      "was the original error message that caused this: '%s'." % e)
                else:
                    self.progress.reset()

                    shutil.move(self.artifact_path, self.artifact_path + '.broken')

                    raise ConfigError("Anvi'o has very bad news for you :( Your migration failed, and anvi'o has no backups to restore your "
                                      "original file. The current artifact is likely in a broken state, and you will unlikely going to be "
                                      "able to use it. So anvi'o renamed it by adding a prefix '.broken' to its file name. We are very sorry "
                                      "for this error (and anvi'o will certainly not put salt on the wound by reminding you that you could "
                                      "have avoided it by using the `--migrate-safely` flag): \"%s\"." % e)

                    print(e)
                    sys.exit(-1)

            # special cases #
            # after this script is done, genome storage changes extension from .h5 to .db
            if self.artifact_type == 'genomestorage' and script_name == 'v4_to_v5':
                self.artifact_path = self.artifact_path[:-3] + '.db'

        if self.safe_mode:
            self.progress.new("Post migration")
            self.progress.update("Removing the backup file ...")
            os.remove(temp_file_path)
            self.progress.end()
        else:
            temp_file_path = None


def main():
    args, unknown = get_args()

    try:
        if len(args.input) > 1 and args.target_version:
            raise ConfigError(f"You have to provide a single file to use with the `--target-version` parameter. "
                              f"You have provided {len(args.input)} :/")

        for artifact_path in args.input:
            args_for_single_artifact = copy.deepcopy(args)
            args_for_single_artifact.artifact_path = artifact_path

            Migrater(args_for_single_artifact).process()

    except ConfigError as e:
        print(e)
        sys.exit(-1)
    except FilesNPathsError as e:
        print(e)
        sys.exit(-2)


def get_args():
    parser = ArgumentParser(description=__description__)

    groupA = parser.add_argument_group('INPUTS', "You will literally give us anything.")
    groupA.add_argument('input', metavar = 'FILE(S)', nargs='+', help = "One or more files to migrate. If you unleash the program on "
                        "many files by running something like `anvi-migrate *.db *.json`, anvi'o will try to migrate "
                        "all of the files that match to your filters.")

    groupB = parser.add_argument_group('SAFETY', "It is up to you. Safe things take much longer and boring. Unsafe things "
                                                 "are fast, fun, and .. well, don't come to us if your computer loses power "
                                                 "or something.")
    groupB.add_argument(*anvio.A('migrate-safely'), **anvio.K('migrate-safely'))
    groupB.add_argument(*anvio.A('migrate-quickly'), **anvio.K('migrate-quickly'))

    groupC = parser.add_argument_group('PARAMETERS OF CONVENIENCE', "This is how anvi'o spoils you.")
    groupC.add_argument(*anvio.A('just-do-it'), **anvio.K('just-do-it'))
    groupC.add_argument(*anvio.A('target-version'), **anvio.K('target-version'))

    return parser.parse_known_args()


if __name__ == '__main__':
    main()

#!/usr/bin/env python
# -*- coding: utf-8

import os
import glob
import argparse

import anvio
import anvio.tables as t
import anvio.terminal as terminal 
import anvio.dbops as dbops
import anvio.db as db
from anvio.errors import ConfigError
from sqlite3 import OperationalError

run = terminal.Run()

parser = argparse.ArgumentParser(description='A simple script to upgrade profile database to version 5 and import all state files into database')
parser.add_argument('-p', '--profile-db', metavar = 'PROFILE_DB',
                    help = 'Profile database')

args = parser.parse_args()

if args.profile_db is None:
	raise ConfigError, "No profile database."

run.info_single("Trying to upgrade profile database...")

try:
    profile_db = db.DB(args.profile_db, '4')
    profile_db.create_table(t.states_table_name, t.states_table_structure, t.states_table_types)   
    profile_db.remove_meta_key_value_pair('version')
    profile_db.set_version('5')
    profile_db.disconnect()
    run.info_single("Database successfully upgraded.")
except OperationalError:
    run.warning("It seems '%s' table already exists in your profile database." % t.states_table_name)
except ConfigError:
    run.warning("It seems your profile database already in version '5'.")

run.info_single("Trying to import all state files in same path with profile databse.")

profile_db_full_path = os.path.join(os.getcwd(), args.profile_db)
dir_name = os.path.dirname(profile_db_full_path)

search_pattern = os.path.join(dir_name, 'state*.json')
state_files = glob.glob(search_pattern)

run.info("State files found in the path", len(state_files))

for state_file_path in state_files:
    path, filename = os.path.split(state_file_path)

    #replace .json and use as key.
    state_id = filename.replace('.json', '')

    run.info_single("Importing '%s'..." % state_id)

    # get json file content
    with open(state_file_path) as f: content = f.read()

    states = dbops.TablesForStates(args.profile_db, '5')
    states.store_state(state_id, content)


#!/bin/bash
source ../00.sh
set -e

get_random_string() {
    echo `python -c 'import string; import random; print "".join(random.choice(string.ascii_lowercase + string.digits) for _ in range(30))'`
}

add_user() {
cat << EOF | sqlite3 test-output/users-data/USERS.db
INSERT INTO "users" VALUES("$1","$2",'$3',"$1@email",'n7YtL2o4bGG6Q',"$4","`get_random_string`",1,NULL,'Some Affiliation','127.0.0.1','user','2015-12-18',NULL);
EOF

mkdir test-output/users-data/userdata/$4

echo "New user: $1 (password: 'test') ..."
}

add_project() {
proj_dir=`get_random_string`
cp -r anvi_server_files/mock_project_directory_01 test-output/users-data/userdata/$2/$proj_dir

cat << EOF | sqlite3 test-output/users-data/USERS.db
INSERT INTO "projects" VALUES("$3", "$proj_dir", "$1", "$4");
EOF

echo "* New project, $3, for user $1 has been created ..."
}

gen_mock_databases() {
cat << EOF | sqlite3 test-output/users-data/USERS.db
PRAGMA foreign_keys=OFF;
BEGIN TRANSACTION;
CREATE TABLE self (key TEXT PRIMARY KEY, value TEXT);
INSERT INTO "self" VALUES('version','1');
CREATE TABLE users (login TEXT PRIMARY KEY, firstname TEXT, lastname TEXT, email TEXT, password TEXT, path TEXT, token TEXT, accepted INTEGER, project TEXT, affiliation TEXT, ip TEXT, clearance TEXT, date TEXT, visit TEXT);
CREATE TABLE projects (name TEXT PRIMARY KEY, path TEXT, user TEXT, description TEXT);
CREATE TABLE metadata (project INTEGER, attribute TEXT, value TEXT);
CREATE TABLE views (name TEXT, user TEXT, project TEXT, public INTEGER, token TEXT);
COMMIT;
EOF

mkdir test-output/users-data/userdata

INFO "Adding users"
add_user "meren" "A. Murat" "Eren" "merens_dir"
add_user "tobi" "Tobias" "Paczian" "tobis_dir"
add_user "ozcan" "Özcan" "Esen" "ozcans_dir"
add_user "tdelmont" "Tom" "Delmont" "toms_dir"
add_user "testuser" "test" "user" "test_dir"

INFO "Adding projects"
add_project "meren" "merens_dir" "m_proj_01" "description_a"
add_project "meren" "merens_dir" "m_proj_02" "description_b"
add_project "tobi" "tobis_dir" "t_proj_01" "description_c"
add_project "ozcan" "ozcans_dir" "o_proj_01" "description_d"
add_project "ozcan" "ozcans_dir" "o_proj_02" "description_e"
add_project "ozcan" "ozcans_dir" "o_proj_03" "description_f"
add_project "ozcan" "ozcans_dir" "o_proj_04" "description_g"
add_project "testuser" "test_dir" "test_proj_01" "description_h"

INFO "Giving user 'tobi' admin credentials"
cat << EOF | sqlite3 test-output/users-data/USERS.db
UPDATE users SET clearance='admin' WHERE login='tobi';
EOF
}




cd ../sandbox

echo $x

if [ $# -eq 0  ]
then
      echo "\

        No arguments supplied. You must call this script with one argument: 'new',
        'mock', or 'continue':

        'new'     : starts the server with an empty users database.
        'mock'    : starts the server from scratch but with some mock users.
        'continue': starts the server with the previously generated data.
        "
        exit -1
fi

if [ $# -gt 1  ]
then
      echo "
        This scripts expect only one argument ('new','continue', or 'mock').
        "
        exit -1
fi


if [ $1 = "new"  ]
then
    INFO "Creating an empty output directory ..."
    rm -rf test-output
    mkdir test-output
    mkdir test-output/users-data
elif [ $1 = "mock" ]
then
    INFO "Creating an empty output directory with mock users and data..."
    rm -rf test-output
    mkdir test-output
    mkdir test-output/users-data
    gen_mock_databases
elif [ $1 = "continue"  ]
then
    if [ ! -d "test-output/users-data"  ]
    then
      echo "
        You asked to continue with the previously generated users directory,
        but none found... Please re-run this script with parameter 'new' or
        'mock' to create one.
        "
        exit -1
    else
        INFO "Attempting to continue with the previously generated directory ..."
    fi
else
      echo "
        Unknown parameter $1 :/ Try 'new','continue', or 'mock'.
        "
        exit -1
fi

INFO "Anvo'o version ..."
anvi-profile --version

INFO "Running anvi'o server ..."
if [ ! -f "smtp_config.ini"  ]
then
    echo "
        SMTP config file is missing! Therefore this script will run the server
        without an SMTP setup. If you would like to run *with* SMTP support,
        please go into the ../sandbox/ directory, create a copy of the file
        'smtp_config_sample.ini' as 'smtp_config.ini' and edit this copy
        according to your SMTP settings before re-running this script.
        "

    anvi-server -U test-output/users-data -P 8080
else
    anvi-server -E smtp_config.ini -U test-output/users-data -P 8080
fi



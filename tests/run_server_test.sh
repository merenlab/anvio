#!/bin/bash
source 00.sh
set -e

cd sandbox

if [ $# -eq 0  ]
then
      echo "\

        No arguments supplied. If you want to start a new server with an emptry
        users database, call this script with parameter 'new', if you want to
        continue with a previously generated users database, use 'continue' as a
        parameter
        "
        exit -1
fi

if [ $# -gt 1  ]
then
      echo "
        This scripts expect only one argument ('new' or 'continue').
        "
        exit -1
fi


if [ $1 = "new"  ]
then
    INFO "Creating an empty output directory ..."
    rm -rf test-output
    mkdir test-output
elif [ $1 = "continue"  ]
then
    if [ ! -d "test-output/users-data"  ]
    then
      echo "
        You asked to continue with the previously generated users directory,
        but none found... Please re-run this script with parameter 'new' to
        create one.
        "
        exit -1
    else
        INFO "Attempting to continue with the previously generated directory ..."
    fi
else
      echo "
        Unknown parameter $1 :/ Try 'new' or 'continue'
        "
        exit -1
fi

if [ ! -f "smtp_config.ini"  ]
then
    echo "
        SMTP config file is missing! Please go into the sandbox/ directory,
        create a copy of 'smtp_config_sample.ini' as 'smtp_config.ini' and
        edit this copy according to your SMTP settings.
        "
fi


INFO "Anvo'o version ..."
anvi-profile --version

INFO "Running anvi'o server ..."
anvi-server -E smtp_config.ini -U test-output/users-data

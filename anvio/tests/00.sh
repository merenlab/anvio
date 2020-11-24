#!/bin/bash

set -e

C() {
    echo -e "\033[1;30m\033[47m$1\033[0m"
}

INFO() {
    echo
    echo
    C ":: $1 ..."
    echo
}

SETUP_WITH_OUTPUT_DIR() {
    if [ -z "$1"  ]
    then
        output_dir="`pwd`/sandbox/test-output"
        rm -rf $output_dir
    else
        output_dir="$1"
    fi

    if [ -z "$2"  ]
    then
        dry_run_controller=""
    else
        dry_run_controller="--dry-run"
    fi

    mkdir -p $output_dir

    files="sandbox"

    INFO "Output directory"
    echo "$output_dir"
    echo

    INFO "Can has interactive displays?"
    if [ -z "$2"  ]
    then
        echo "Yes, anvi'o will try to show you interactive displays that will require you to come back to the terminal and press CTRL+C to continue"
        echo
    else
        echo "No interactive displays"
        echo
    fi


    INFO "Anvo'o version"
    anvi-profile --version
}

SHOW_FILE() {
    echo
    cat $1
}

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
        output_dir="$1/test-output"
    fi

    mkdir -p $output_dir

    files="sandbox"

    INFO "Output directory"
    echo "$output_dir"
    echo

    INFO "Anvo'o version"
    anvi-profile --version
}

#!/bin/bash

set -e

O() {
    echo -e "\033[1;1m\033[41m:: $1 ...\033[0m"
}

C() {
    echo -e "\033[1;30m\033[47m:: $1 ...\033[0m"
}

INFO() {
    echo
    echo
    C "$1"
    echo
}

OUTPUT() {
    echo
    echo
    O "$1"
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

    if [ "$2" == "interactive" ]
    then
        dry_run_controller=""
    else
        dry_run_controller="--dry-run"
    fi

    if [ -z "$3" ]
    then
        thread_controller=""
    else
        thread_controller="--num-threads $3"
    fi

    mkdir -p $output_dir

    files="sandbox"

    INFO "Output directory"
    echo "$output_dir"
    echo

    INFO "Can has interactive displays?"
    if [ "$2" == "interactive"  ]
    then
        echo "Yes, anvi'o will try to show you interactive displays that will require you to come back to the terminal and press CTRL+C to continue"
        echo
    else
        echo "No interactive displays"
        echo
    fi

    INFO "Do we use threadz yo?"
    if [ -z "$3" ]
    then
        echo "NO THREADS FOR YOU! (well, ok. just one.)"
        echo
    else
        echo "According to our sources, we will multi-thread all available commands with $thread_controller"
        echo
    fi


    INFO "Anvi'o version"
    anvi-profile --version
}

SHOW_FILE() {
    OUTPUT $1
    head -n 10 $1 | anvi-script-tabulate
}

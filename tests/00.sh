#!/bin/bash

C() {
    echo -e "\033[1;30m\033[47m$1\033[0m"
}

INFO() { 
    echo
    echo
    C ":: $1 ..."
    echo
} 

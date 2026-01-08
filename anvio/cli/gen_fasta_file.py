#!/usr/bin/env python
# -*- coding: utf-8

import sys
import os

import anvio
import anvio.utils as utils
import anvio.terminal as terminal
from anvio.errors import ConfigError, FilesNPathsError
import anvio.filesnpaths as filesnpaths

__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['ahenoch']
__provides__ = ['fasta-txt']
__requires__ = ['fasta']
__description__ = "Create the fasta.txt file"


def main():
    args = get_args()

    try:  
        run_program()
    except ConfigError as e:
        print(e)
        sys.exit(-1)
    except FilesNPathsError as e:
        print(e)
        sys.exit(-2)


def run_program():
    args = get_args()
    run = terminal.Run()
    
    filesnpaths.is_output_file_writable(args.fasta_text_file, ok_if_exists=False)

    fasta_text_file_dict = {}

    for file in os.listdir(args.input_dir):
        file_path = os.path.abspath(os.path.join(args.input_dir, file))
        if filesnpaths.is_file_fasta_formatted(file_path, dont_raise=True):
            file_name = os.path.basename(file_path).rsplit('.')[0]
            fasta_text_file_dict[file_name] = {'path': file_path}

    utils.store_dict_as_TAB_delimited_file(fasta_text_file_dict, args.fasta_text_file, key_header='name')

    run.info_single("Successfully created fasta-text-file.")


def get_args():
    from anvio.argparse import ArgumentParser
    parser = ArgumentParser(description=__description__)

    parser.add_argument(*anvio.A('input-dir'), **anvio.K('input-dir', {'required': True}))
    parser.add_argument(*anvio.A('fasta-text-file'), **anvio.K('fasta-text-file', {'required': True}))

    return parser.get_args(parser)


if __name__ == '__main__':
    main()
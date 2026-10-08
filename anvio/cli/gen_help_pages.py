#!/usr/bin/env python

import sys
from pathlib import Path

import anvio
import anvio.terminal as terminal

from anvio.argparse import ArgumentParser
from anvio.docs import AnvioDocs, DocumentationData, HelpPagesRenderer
from anvio.errors import ConfigError, FilesNPathsError

__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['meren', 'Jessica-Pan']
__description__ = "Generate a static web site for anvi'o help pages"


@terminal.time_program
def main():
    args = get_args()

    try:
        if args.from_data:
            path = Path(args.from_data)
            docs_dataset = DocumentationData.load(path)
        else:
            docs = AnvioDocs(args)
            docs_dataset = docs.collect()

        HelpPagesRenderer(
            dataset=docs_dataset,
            output_directory=Path(args.output_dir or 'ANVIO-HELP'),
        ).generate()

    except ConfigError as e:
        print(e)
        sys.exit(-1)
    except FilesNPathsError as e:
        print(e)
        sys.exit(-2)


def get_args():
    parser = ArgumentParser(description=__description__)

    parser.add_argument(*anvio.A('output-dir'), **anvio.K('output-dir'))
    parser.add_argument('--from-data', metavar='DOCUMENTATION_JSON',
                        help="Rebuild help pages from an exported documentation.json and its adjacent images directory, "
                             "without reading CLI metadata or documentation sources. Use a separate output directory.")

    return parser.get_args(parser)


if __name__ == '__main__':
    main()

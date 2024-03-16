#!/usr/bin/env python
# -*- coding: utf-8

import os
import sys

import anvio
import anvio.utils as utils
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths

from anvio.errors import ConfigError, FilesNPathsError

__author__ = "Developers of anvi'o (see AUTHORS.txt)"
__copyright__ = "Copyleft 2015-2024, the Meren Lab (http://merenlab.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['Ge0rges']
__provides__ = ["hmm-source"]
__requires__ = ["hmm"]
__description__ = ("You give this program one or more HMM files from `hmmbuild`, and it generates "
                   "an anvi'o compatible HMM directory to be used with `anvi-run-hmms`")


run = terminal.Run()
progress = terminal.Progress()
pp = terminal.pretty_print


def get_attribute_from_hmm_file(file_path, attribute):
    filesnpaths.is_file_exists(file_path)
    value = None
    with open(file_path) as hmm:
        for line in hmm.readlines():
            if line.startswith(attribute):
                value = [f.strip() for f in line.split(attribute) if len(f)][0]
                break

    if value is None:
        raise ConfigError(f"In the HMM file {file_path} we did not find the attribute {attribute}. This is probably an "
                          f"issue with your file.")
    return value


def main(args):
    A = lambda x: args.__dict__[x] if x in args.__dict__ else None
    hmm_files = A('hmm-list')
    hmm_sources = A('hmm-sourcest')

    # Check input exists
    if hmm_files is None:
        raise ConfigError("You must provide at least one path to a HMM file using `--hmm-list` :/")

    hmm_files = [hmm_files.strip() for e in hmm_files]
    for f in hmm_files:
        if not filesnpaths.is_file_exists(f):
            raise ConfigError(f"The path you provided {f} does not point to a file that exists.")

    if len(hmm_files) == len(set(hmm_files)):
        run.warning("You've specified at least one file more than once in `--hmm-list` this is Ok for Anvi'o but you"
                    "probably should restart the program with a corrected list. "
                    "Otherwise, you might get undefined behaviour downstream. We'll wait a second before continuing.")
        time.sleep(1)

    # Check sources
    if hmm_sources is None:
        raise ConfigError("You must provide at least one source for your HMMs using `--hmm-source` :/")

    if len(hmm_sources) > 1 and not len(hmm_files) == len(hmm_sources):
        raise ConfigError("You must provide either only one source for your HMMs or as many sources as there are HMMs "
                          "using `--hmm-source` :/")
    elif len(hmm_sources) == 1:
        hmm_sources = [hmm_sources[0]] * len(hmm_files)

    # Check output
    output_directory_path = A('output_directory')
    filesnpaths.check_output_directory(output_directory_path)

    # Update the user
    run.info('HMM files to work with', ', '.join(hmm_files))
    run.info('The output directory', output_directory_path)

    # Parse all the HMM files
    progress.new("Parsing", progress_total_items=len(hmm_files))
    progress.update('...')
    data_dict = {}
    for hmm_file, hmm_source in zip(hmm_files, hmm_sources):
        progress.update(hmm_file + ' ...', increment=True)

        data_dict[hmm_file] = {}
        data_dict[hmm_file]['ga'] = get_attribute_from_hmm_file(hmm_file, 'GA ')
        data_dict[hmm_file]['gene'] = get_attribute_from_hmm_file(hmm_file, 'NAME')
        data_dict[hmm_file]['accession'] = get_attribute_from_hmm_file(hmm_file, 'ACC')
        data_dict[hmm_file]['source'] = hmm_source
    progress.end()

    # Generate the Anvi'o HMM
    progress.new("Generating the contents of the HMM directory", progress_total_items=len(data_dict))
    progress.update('...')

    # Create the output directory
    filesnpaths.gen_output_directory(output_directory_path)
    J = lambda x: os.path.join(output_directory_path, x)
    W = lambda p, c: open(J(p), 'w').write(f'{c}\n')

    # Concatenate and compress the genes.hmm
    utils.concatenate_files(J('genes.hmm'), [hmm_file for hmm_file in hmm_files])
    utils.gzip_compress_file(J('genes.hmm'))

    # Generate genes output
    with open(J('genes.txt'), 'w') as genestxt:
        genestxt.write("gene\taccession\thmmsource\n")
        for e in data_dict.values():
            genestxt.write(f"{e['gene']}\t{e['accession']}\te['source']\n")

    # kind
    W('kind.txt', os.path.basename(output_directory_path))
    W('noise_cutoff_terms.txt', '--cut_ga')
    W('reference.txt', "Anvi'o User, http://localhost")
    W('target.txt', "AA:GENE")

    progress.end()

    run.info_single(f"Congratulations. Your anvi'o formatted HMM directory for "
                    f"{terminal.pluralize('HMM file', len(hmm_files))} is ready "
                    f"to be used with `anvi-run-hmms` (all you need to do is to provide the path "
                    f"to your new directory using the `--hmm-profile-dir` parameter).", nl_before=1, nl_after=1)


if __name__ == '__main__':
    from anvio.argparse import ArgumentParser
    parser = ArgumentParser(description=__description__)

    parser.add_argument('--hmm-list', nargs='+', metavar='FILES',
                        help="One or more paths to HMM files. These should have been generated using hmmbuild" 
                             " in HMMER3, or otherwise be in the same format.")
    parser.add_argument('--hmm-source', nargs='+', metavar='SOURCES',
                        help="The source of these HMMs (e.g. pfam.xfam.org). Specify either one or as many as"
                             " there are HMM files. Specify this list in the same order as the files in `--hmm-list`")
    parser.add_argument('-o', '--output-directory', metavar='PATH', help="Output directory for the "
                "anvi'o-formatted HMMs. Choose the name wisely as this will be the name that will "
                "appear in the contigs database after you provide it with `-H` flag to `anvi-run-hmms`. "
                "We suggest you to use a name that does not include any special characters or punctuation"
                "(such as space, question mark, comma, etc.).")

    args = parser.get_args(parser)

    try:
        main(args)
    except ConfigError as e:
        print(e)
        sys.exit(-1)
    except FilesNPathsError as e:
        print(e)
        sys.exit(-2)

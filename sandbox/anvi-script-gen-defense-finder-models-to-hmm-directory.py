#!/usr/bin/env python
# -*- coding: utf-8

import os
import sys
import tarfile
import requests

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
__description__ = ("This program generates an anvi'o compatible HMM directory to be used with `anvi-run-hmms` "
                   "from the MDMParis Defense Finder Models.")


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
        raise ConfigError(f"In the HMM file {file_path} we did not find the attribute {attribute}. "
                          f"You should report this issue on the Defense-Finder-Models GitHub.")
    return value


def download_latest_release():
    progress.new("Downloading", progress_total_items=len(data_dict))
    progress.update('...')

    repo_url = "https://api.github.com/repos/mdmparis/defense-finder-models/releases/latest"
    response = requests.get(repo_url)
    data = response.json()
    download_url = [asset["browser_download_url"] for asset in data["assets"] if asset["name"].endswith(".tar.gz")][0]
    response = requests.get(download_url)
    with open("latest_release.tar.gz", "wb") as f:
        f.write(response.content)
    with tarfile.open("latest_release.tar.gz", "r:gz") as tar:
        tar.extractall()
    os.remove("latest_release.tar.gz")

    progress.end()


def main(args):
    A = lambda x: args.__dict__[x] if x in args.__dict__ else None

    # Check output
    output_directory_path = A('output_directory')
    filesnpaths.check_output_directory(output_directory_path)

    # Update the user
    run.info('HMM files to work with', ', '.join(hmm_files))
    run.info('The output directory', output_directory_path)

    # Download Defense Finder Models
    download_latest_release()
    source_folder = "defense-finder-models"
    filesnpaths.is_file_exists("defense-finder-models/metadata.yml")

    # Parse all the HMM files
    progress.new("Parsing", progress_total_items=len(os.listdir(os.path.join(source_folder, "profiles"))))
    progress.update('...')
    data_dict = {}

    for hmm_name in os.listdir(os.path.join(source_folder, "profiles")):
        hmm_path = os.path.join(source_profile_folder, hmm_name)

        progress.update(hmm_name + ' ...', increment=True)

        data_dict[hmm_name] = {}
        data_dict[hmm_name]['ga'] = get_attribute_from_hmm_file(hmm_file, 'GA ')
        data_dict[hmm_name]['gene'] = get_attribute_from_hmm_file(hmm_file, 'NAME')
        data_dict[hmm_name]['accession'] = get_attribute_from_hmm_file(hmm_file, 'ACC')
        data_dict[hmm_name]['source'] = hmm_source
        data_dict[hmm_name]['hmm_fp'] = hmm_path
    progress.end()

    # Extract version number from metadata.yml
    version_number = None
    with open(os.path.join(source_folder, "metadata.yml"), "r") as f:
        for line in f:
            if "vers" in line:
                version_number = line.split()[1]
                break

    # Generate the Anvi'o HMM
    progress.new("Generating the contents of the HMM directory", progress_total_items=len(data_dict))
    progress.update('...')

    # Create the output directory
    filesnpaths.gen_output_directory(output_directory_path)
    J = lambda x: os.path.join(output_directory_path, x)
    W = lambda p, c: open(J(p), 'w').write(f'{c}\n')

    # Concatenate and compress the genes.hmm
    utils.concatenate_files(J('genes.hmm'), [data_dict[p]['hmm_fp'] for p in data_dict], remove_concatenated_files=True)
    utils.gzip_compress_file(J('genes.hmm'))

    # Generate genes output
    with open(J('genes.txt'), 'w') as genestxt:
        genestxt.write("gene\taccession\thmmsource\n")
        for e in data_dict.values():
            genestxt.write(f"{e['gene']}\t{e['accession']}\t{e['source']}\n")

    # kind
    W('kind.txt', 'DefenseFinder')
    W('noise_cutoff_terms.txt', '--cut_ga')
    W('reference.txt', f"DefenseFinder v{version_number}")
    W('target.txt', "AA:GENE")

    # Delete files
    os.removedirs(source_folder)
    progress.end()

    run.info_single(f"Congratulations. Your anvi'o formatted HMM directory for "
                    f"the MDMParis/Defense-Finder-Models is ready "
                    f"to be used with `anvi-run-hmms` (all you need to do is to provide the path "
                    f"to your new directory using the `--hmm-profile-dir` parameter).", nl_before=1, nl_after=1)


if __name__ == '__main__':
    from anvio.argparse import ArgumentParser
    parser = ArgumentParser(description=__description__)

    args = parser.get_args(parser)

    try:
        main(args)
    except ConfigError as e:
        print(e)
        sys.exit(-1)
    except FilesNPathsError as e:
        print(e)
        sys.exit(-2)

#!/usr/bin/env python
# -*- coding: utf-8

import os
import sys
import json
import shutil
import tarfile

import anvio
import anvio.utils as utils
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths

from anvio.errors import ConfigError, FilesNPathsError

__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['Ge0rges']
__provides__ = ["hmm-source"]
__requires__ = ["hmm-file"]
__description__ = ("This program generates an anvi'o compatible HMM directory to be used with `anvi-run-hmms` "
                   "from the MDMParis Defense Finder Models.")




def main():
    try:
        run_program()
    except ConfigError as e:
        print(e)
        sys.exit(-1)
    except FilesNPathsError as e:
        print(e)
        sys.exit(-2)


def run_program():
    _ = get_args()
    run = terminal.Run()
    progress = terminal.Progress()

    # Check output
    output_directory_path = "DefenseFinder_HMM"
    filesnpaths.check_output_directory(output_directory_path)

    # Download Defense Finder Models
    source_folder = "defense-finder-models"
    source_profile_folder = os.path.join(source_folder, "profiles")

    progress.new("Downloading", progress_total_items=0)
    progress.update('...')

    repo_url = "https://api.github.com/repos/mdmparis/defense-finder-models/releases/latest"
    response = utils.get_remote_file_content(repo_url)
    data = json.loads(response)
    download_url = [asset["browser_download_url"] for asset in data["assets"] if asset["name"].endswith(".tar.gz")][0]
    response = utils.download_file(download_url, "latest_release.tar.gz")

    with tarfile.open("latest_release.tar.gz", "r:gz") as tar:
        tar.extractall()
    os.remove("latest_release.tar.gz")

    progress.end()

    filesnpaths.is_file_exists("defense-finder-models/metadata.yml")

    # Parse all the HMM files
    progress.new("Parsing", progress_total_items=len(os.listdir(source_profile_folder)))
    progress.update('...')
    data_dict = {}

    for hmm_name in os.listdir(source_profile_folder):
        hmm_path = os.path.join(source_profile_folder, hmm_name)

        progress.update(hmm_name + ' ...', increment=True)

        try:
            acc = utils.get_attribute_from_hmm_file(hmm_path, 'ACC')
        except ValueError:
            continue

        try:
            acc = utils.get_attribute_from_hmm_file(hmm_path, "NAME") if acc is None else acc
            ga = utils.get_attribute_from_hmm_file(hmm_path, 'GA ')
            name = utils.get_attribute_from_hmm_file(hmm_path, 'NAME')

        except ValueError:
            raise ConfigError(f"In the HMM file {hmm_path} we did not find one of the attributes we expected to find. "
                              f"You should report this issue on the Defense-Finder-Models GitHub.")

        data_dict[hmm_name] = {}
        data_dict[hmm_name]['ga'] = ga
        data_dict[hmm_name]['gene'] = name
        data_dict[hmm_name]['accession'] = acc
        data_dict[hmm_name]['source'] = "https://github.com/mdmparis/defense-finder-models"
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
    shutil.rmtree(source_folder)
    progress.end()

    run.info_single("Congratulations. Your anvi'o formatted HMM directory for "
                    "the MDMParis/Defense-Finder-Models is ready at 'DefenseFinder_HMM' "
                    "to be used with `anvi-run-hmms` (all you need to do is to provide the path "
                    "to your new directory using the `--hmm-profile-dir` parameter).", nl_before=1, nl_after=1)


def get_args():
    from anvio.argparse import ArgumentParser
    parser = ArgumentParser(description=__description__)

    return parser.get_args(parser)


if __name__ == '__main__':
    main()

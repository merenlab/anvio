#!/usr/bin/env python

import os
import sys
import collections

import anvio
import anvio.utils as utils
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths

from anvio.errors import ConfigError, FilesNPathsError
from anvio.metabolism.context import KeggContext

__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['meren', 'ge0rges', 'ivujic']
__requires__ = ["kofam-accession", "kegg-data"]
__provides__ = ["hmm-source"]
__description__ = ("You give this program one or more KOfam accession ids, and it generates "
                   "an anvi'o compatible HMM directory to be used with `anvi-run-hmms` "
                   "by extracting them from your local KEGG setup.")


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
    args = get_args()
    run = terminal.Run()
    progress = terminal.Progress()

    A = lambda x: args.__dict__[x] if x in args.__dict__ else None
    kofam_accessions_list = A('kofam_accessions_list')
    kofam_accessions_file = A('kofam_accessions_file')

    if kofam_accessions_list and kofam_accessions_file:
        raise ConfigError("You should either provide KOfam accession ids through the command line, or "
                          "list them in a file. Doing both is just... unnecessary.")

    if kofam_accessions_list:
        kofam_accession_ids = kofam_accessions_list
    elif kofam_accessions_file:
        filesnpaths.is_file_tab_delimited(kofam_accessions_file, expected_number_of_fields=1)
        kofam_accession_ids = [a.strip() for a in open(kofam_accessions_file).readlines()]
    else:
        raise ConfigError("You should provide *some* KOfam accession ids to this program :/")

    # sanity check
    not_kofam_accession_ids = [k for k in kofam_accession_ids if not k.startswith("K")]
    if len(not_kofam_accession_ids):
        raise ConfigError(f"The following accessions do not appear to be from KOfam because they do not "
                          f"start with \"K\", please double check the following: {','.join(not_kofam_accession_ids)}")

    # Use KeggContext to find the data
    kegg_context = KeggContext(args)
    kegg_context.run = run
    kegg_context.progress = progress
    ko_list_path = kegg_context.ko_list_file_path
    kofam_hmm_path = kegg_context.kofam_hmm_file_path

    if not os.path.exists(ko_list_path):
        raise ConfigError(f"The KOfam ko_list file was not found at '{ko_list_path}'. "
                          "Please make sure you have run `anvi-setup-kegg-data` or "
                          "provided the correct `--kegg-data-dir`.")

    if not os.path.exists(kofam_hmm_path):
        raise ConfigError(f"The KOfam HMM file was not found at '{kofam_hmm_path}'. "
                          "Please make sure you have run `anvi-setup-kegg-data` or "
                          "provided the correct `--kegg-data-dir`.")

    # Load thresholds and definitions using KeggContext methods
    kegg_context.setup_ko_dict(exclude_threshold=False)
    ko_dict = kegg_context.ko_dict

    # Check for stray KOs
    stray_ko_dict = {}
    if os.path.exists(kegg_context.stray_ko_thresholds_file):
        kegg_context.setup_stray_ko_dict()
        stray_ko_dict = kegg_context.stray_ko_dict

    # determine output directory
    kofam_accession_ids = sorted([e.strip() for e in set(kofam_accession_ids)])
    if not A('output_directory'):
        if len(kofam_accession_ids) == 1:
            output_directory_path = os.path.abspath(kofam_accession_ids[0])
        else:
            output_directory_path = os.path.abspath(f"{kofam_accession_ids[0]}_and_{len(kofam_accession_ids) - 1}_others")
    else:
        output_directory_path = os.path.abspath(A('output_directory'))

    failed_accession_ids = set([])

    run.info('KOfam accessions to work with', ', '.join(kofam_accession_ids))
    run.info('The output directory', output_directory_path)

    filesnpaths.check_output_directory(output_directory_path)

    data_dict = {}

    progress.new("Extracting HMMs", progress_total_items=len(kofam_accession_ids))
    for kofam_accession in kofam_accession_ids:
        progress.update(kofam_accession + ' ...', increment=True)

        # Decide which dictionary and HMM file to use
        current_ko_info = None
        current_hmm_file = kofam_hmm_path

        if stray_ko_dict and kofam_accession in stray_ko_dict:
            current_ko_info = stray_ko_dict[kofam_accession]
            current_hmm_file = kegg_context.stray_ko_hmm_file_path
        elif kofam_accession in ko_dict:
            current_ko_info = ko_dict[kofam_accession]
            current_hmm_file = kofam_hmm_path
        else:
            failed_accession_ids.add(kofam_accession)
            continue

        # extract HMM
        fp = filesnpaths.get_temp_file_path()
        log_file = filesnpaths.get_temp_file_path()
        cmd = ['hmmfetch', '-o', fp, current_hmm_file, kofam_accession]
        if utils.run_command(cmd, log_file):
            failed_accession_ids.add(kofam_accession)
            if os.path.exists(fp):
                os.remove(fp)
            if os.path.exists(log_file):
                os.remove(log_file)
            continue

        if os.path.exists(log_file):
            os.remove(log_file)

        # Patch HMM with ACC, DESC, and GA threshold
        threshold = current_ko_info.get('threshold')
        definition = current_ko_info.get('definition', kofam_accession)
        if ' [' in definition:
            definition = definition.split(' [')[0]

        hmm_lines = open(fp).readlines()
        new_hmm_lines = []

        # attributes we want to ensure are in the HMM file
        attributes_to_add = collections.OrderedDict()
        if not any(line.startswith('ACC ') for line in hmm_lines):
            attributes_to_add['ACC '] = f"ACC   {kofam_accession}\n"
        if not any(line.startswith('DESC ') for line in hmm_lines):
            attributes_to_add['DESC '] = f"DESC  {definition}\n"
        if threshold and threshold != '-' and not any(line.startswith('GA ') for line in hmm_lines):
            attributes_to_add['GA '] = f"GA    {threshold} {threshold};\n"

        ga_added = 'GA ' not in attributes_to_add
        for line in hmm_lines:
            new_hmm_lines.append(line)

            if line.startswith('NAME '):
                for attr in ['ACC ', 'DESC ']:
                    if attr in attributes_to_add:
                        new_hmm_lines.append(attributes_to_add[attr])
                        del attributes_to_add[attr]

            if not ga_added and (line.startswith('STATS') or line.startswith('HMM  ')):
                if 'GA ' in attributes_to_add:
                    new_hmm_lines.append(attributes_to_add['GA '])
                    del attributes_to_add['GA ']
                ga_added = True

        with open(fp, 'w') as f:
            f.writelines(new_hmm_lines)

        data_dict[kofam_accession] = {}
        data_dict[kofam_accession]['gene'] = kofam_accession
        data_dict[kofam_accession]['accession'] = kofam_accession
        data_dict[kofam_accession]['temp_file_path'] = fp

    progress.end()

    if len(failed_accession_ids):
        if not len(failed_accession_ids) == len(kofam_accession_ids):
            run.warning(f"Anvi'o couldn't find some of your accession ids :/ But it will continue to build "
                        f"the HMM directory with {len(kofam_accession_ids) - len(failed_accession_ids)} of "
                        f"{len(kofam_accession_ids)} accessions. Here is the list of those that failed FYI: "
                        f"{', '.join(failed_accession_ids)}.")

            for accession in failed_accession_ids:
                kofam_accession_ids.remove(accession)
        else:
            raise ConfigError("Anvi'o couldn't find any of your accession ids in the local KOfam setup. "
                              "Please make sure your accession ids are correct and that your KEGG "
                              "data directory is properly set up.")

    progress.new("Generating the contents of the HMM directory", progress_total_items=len(data_dict))
    progress.update('...')

    # create the output dir
    filesnpaths.gen_output_directory(output_directory_path)
    J = lambda x: os.path.join(output_directory_path, x)
    W = lambda p, c: open(J(p), 'w').write(f'{c}\n')

    # concatenate and compress the genes.hmm
    utils.concatenate_files(J('genes.hmm'), [data_dict[p]['temp_file_path'] for p in data_dict], remove_concatenated_files=True)
    utils.gzip_compress_file(J('genes.hmm'))

    # generate genes output
    with open(J('genes.txt'), 'w') as genestxt:
        genestxt.write("gene\taccession\thmmsource\n")
        for e in data_dict.values():
            genestxt.write(f"{e['gene']}\t{e['accession']}\tKOfam\n")

    # kind
    W('kind.txt', os.path.basename(output_directory_path))
    W('noise_cutoff_terms.txt', '--cut_ga')
    W('reference.txt', "10.1093/bioinformatics/btz859")
    W('target.txt', "AA:GENE")

    progress.end()

    run.info_single(f"Congratulations. Your anvi'o formatted HMM directory for "
                    f"{terminal.pluralize('KOfam accession id', len(kofam_accession_ids))} is ready "
                    f"to be used with `anvi-run-hmms` (all you need to do is to provide the path "
                    f"to your new directory using the `--hmm-profile-dir` parameter).", nl_before=1, nl_after=1)


def get_args():
    from anvio.argparse import ArgumentParser
    parser = ArgumentParser(description=__description__)

    parser.add_argument('--kofam-accessions-list', nargs='+', help="One or more KOfam accession IDs "
                "(such as K00001). If you have multiple accessions, you can separate them from "
                "each other with a space. If you have too many, consider using the "
                "`--kofam-accessions-file` parameter instead.", metavar='KOFAM_ACCESSION')
    parser.add_argument('--kofam-accessions-file', help="A single column text file where each column "
                "is a single KOfam accession ID (such as K00001). You may have as many accession "
                "ids as you like in this file.", metavar='FILE')
    parser.add_argument('--kegg-data-dir', metavar='PATH', help="Path to the KEGG data directory "
                "set up by `anvi-setup-kegg-data`. If not provided, the default location will be used.")
    parser.add_argument('--include-nt-KOs', action='store_true', default=False, help="By default this "
                "program will skip KOfams that do not have a bitscore threshold defined by KEGG. "
                "But you can include them using this flag.")
    parser.add_argument('-O', '--output-directory', metavar='PATH', help="Output directory for the "
                "anvi'o-formatted HMMs. Choose the name wisely as this will be the name that will "
                "appear in the contigs database after you provide it with `-H` flag to `anvi-run-hmms`.")

    return parser.get_args(parser)


if __name__ == '__main__':
    main()

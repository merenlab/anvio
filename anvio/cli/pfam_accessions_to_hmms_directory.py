#!/usr/bin/env python
# -*- coding: utf-8

import os
import sys

import anvio
import anvio.utils as utils
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths

from anvio.errors import ConfigError, FilesNPathsError

__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['meren', 'ge0rges']
__requires__ = ["pfam-accession"]
__provides__ = ["hmm-source"]
__description__ = ("You give this program one or more PFAM accession ids, and it generates "
                   "an anvi'o compatible HMM directory to be used with `anvi-run-hmms`")


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
    pfam_accessions_list = A('pfam_accessions_list')
    pfam_accessions_file = A('pfam_accessions_file')

    if pfam_accessions_list and pfam_accessions_file:
        raise ConfigError("You should either provide PFAM accession ids thorugh the command line, or "
                          "list them in a file. In theory you should be able to do both indeed, but in "
                          "practice why should you, really?")

    if pfam_accessions_list:
        pfam_accession_ids = pfam_accessions_list
    elif pfam_accessions_file:
        filesnpaths.is_file_tab_delimited(pfam_accessions_file, expected_number_of_fields=1)
        pfam_accession_ids = [a.strip() for a in open(pfam_accessions_file).readlines()]
    else:
        raise ConfigError("You should provide *some* PFAM accession ids to this program :/")

    utils.sanity_check_pfam_accessions(pfam_accession_ids)

    output_directory_path = A('output_directory') or os.path.abspath('./UNKNOWN_HMMS_FROM_PFAM')

    pfam_accession_ids = [e.strip() for e in set(pfam_accession_ids)]
    failed_accession_ids = set([])

    run.info('PFAM accessions to work with', ', '.join(pfam_accession_ids))
    run.info('The output directory', output_directory_path)

    filesnpaths.check_output_directory(output_directory_path)

    data_dict = {}

    progress.new("Downloading", progress_total_items=len(data_dict))
    progress.update('...')
    for pfam_accession in pfam_accession_ids:
        progress.update(pfam_accession + ' ...', increment=True)

        fp = filesnpaths.get_temp_file_path()

        try:
            # this downloads a gzipped model,
            utils.download_file(f"https://www.ebi.ac.uk/interpro/wwwapi/entry/pfam/{pfam_accession}?annotation=hmm", fp + '.gz')

            # which we need to decompress before moving on.
            utils.gzip_decompress_file(fp + '.gz')
        except:
            failed_accession_ids.add(pfam_accession)
            continue

        try:
            data_dict[pfam_accession] = {}
            data_dict[pfam_accession]['ga'] = utils.get_attribute_from_hmm_file(fp, 'GA ')
            data_dict[pfam_accession]['gene'] = utils.get_attribute_from_hmm_file(fp, 'NAME')
            data_dict[pfam_accession]['accession'] = utils.get_attribute_from_hmm_file(fp, 'ACC')
            data_dict[pfam_accession]['temp_file_path'] = fp

        except ValueError:
            run.warning(f"The PFAM accession {pfam_accession} could not be included due to a missing attribute. "
                        f"Anvi'o will continue without it")

    progress.end()

    if len(failed_accession_ids):
        if not len(failed_accession_ids) == len(pfam_accession_ids):
            run.warning(f"Anvi'o couldn't download some of your accession ids :/ But it will continue to build "
                        f"the HMM directory with {len(pfam_accession_ids) - len(failed_accession_ids)} of "
                        f"{len(pfam_accession_ids)} accessions. Here is the list of those that failed FYI: "
                        f"{', '.join(failed_accession_ids)}.")

            for accession in failed_accession_ids:
                pfam_accession_ids.remove(accession)
        else:
            raise ConfigError("Anvi'o couldn't download any of your accession ids. This could be due to bad "
                              "accession ids that don't exist in the PFAM database, or a bad connection"
                              "to the server :/")

    # make sure all models have the GA cutoff defined
    pfam_hmms_without_ga = [p for p in pfam_accession_ids if not data_dict[p]['ga']]
    if len(pfam_hmms_without_ga):
        for pfam_accession in data_dict:
            os.remove(data_dict[pfam_accession]['temp_file_path'])

        raise ConfigError(f"Not all PFAM accession ids you are interested in seem to have 'GA' "
                          f"noise cutoff defined :/ This script is only able to setup "
                          f"anvi'o HMM directories from data_dict that have that cutoff defined. "
                          f"So this is kind of the end of the road for us. Here are the PFAM accessions "
                          f"that violate this: {', '.join(pfam_hmms_without_ga)}.")


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
            genestxt.write(f"{e['gene']}\t{e['accession']}\tpfam.xfam.org\n")

    # kind
    W('kind.txt', os.path.basename(output_directory_path))
    W('noise_cutoff_terms.txt', '--cut_ga')
    W('reference.txt', "Anvi'o User, http://localhost")
    W('target.txt', "AA:GENE")

    progress.end()

    run.info_single(f"Congratulations. Your anvi'o formatted HMM directory for "
                    f"{terminal.pluralize('PFAM accession id', len(pfam_accession_ids))} is ready "
                    f"to be used with `anvi-run-hmms` (all you need to do is to provide the path "
                    f"to your new directory using the `--hmm-profile-dir` parameter).", nl_before=1, nl_after=1)


def get_args():
    from anvio.argparse import ArgumentParser
    parser = ArgumentParser(description=__description__)

    parser.add_argument('--pfam-accessions-list', nargs='+', help="One or more PFAM accession IDs "
                "(such as PF14437). If you have multiple accessions, you can separate them from "
                "each other with a space. If you have too many, consider using the "
                "`--pfam-accessions-file` parameter instead.", metavar='PFAM_ACCESSION')
    parser.add_argument('--pfam-accessions-file', help="A single column text file where each column "
                "is a single PFAM accession ID (such as PF14437). You may have as many accession "
                "ids as you like in this file.", metavar='FILE')
    parser.add_argument('-O', '--output-directory', metavar='PATH', help="Output directory for the "
                "anvi'o-formatted HMMs. Choose the name wisely as this will be the name that will "
                "appear in the contigs database after you provide it with `-H` flag to `anvi-run-hmms`. "
                "We suggest you use a name that does not include any of those millenial characters "
                "(like space, question mark, comma, and such, you know).")

    return parser.get_args(parser)


if __name__ == '__main__':
    main()

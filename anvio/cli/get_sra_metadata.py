#!/usr/bin/env python

import sys

import anvio
import anvio.sra as sra
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths

from anvio.errors import ConfigError


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['FlorianTrigodet']
__provides__ = ['sra-metadata-txt']
__requires__ = ['samples-txt']
__description__ = ("Ask NCBI what it knows about a set of SRA run accessions, and keep the answers in a "
                   "TAB-delimited file. Anvi'o workflows that download reads from the SRA need to know "
                   "whether each run is paired-end or single-end, whether it came off a short-read or a "
                   "long-read instrument, and how big it is, all before they can start. They will build this "
                   "file themselves the first time they need it; this program exists so you can build it "
                   "ahead of time (handy if the computer that runs your workflow has no internet access), "
                   "and so you can see what anvi'o thinks before it acts on it")


run = terminal.Run()


def main():
    try:
        run_program()
    except ConfigError as e:
        print(e)
        sys.exit(-1)


def run_program():
    args = get_args()

    accessions = get_accessions(args)

    entries = sra.get_metadata_for_accessions(accessions,
                                              cache_path=args.output_file,
                                              source_of_accessions=describe_where_accessions_came_from(args),
                                              run=run)

    sra.warn_about_long_read_technologies(entries, run=run, cache_path=args.output_file)

    read_types = {}
    for entry in entries.values():
        read_types[entry['read_type']] = read_types.get(entry['read_type'], 0) + 1

    total_gb = sum(sra.predict_peak_disk_usage_in_gb(e) for e in entries.values())

    run.warning(None, header="WHAT ANVI'O LEARNED", lc="green")
    run.info('SRA runs', len(entries))
    for read_type in sorted(read_types):
        run.info(f"  {read_type}", read_types[read_type])
    run.info('Peak disk space if all were downloaded at once', f"{total_gb:.1f} GB")
    run.info('Metadata file', args.output_file, nl_after=1)

    if sra.SINGLE_END_SHORT_READS in read_types:
        run.warning(f"{terminal.pluralize('run', read_types[sra.SINGLE_END_SHORT_READS])} in this set "
                    f"{'are' if read_types[sra.SINGLE_END_SHORT_READS] > 1 else 'is'} single-end short reads, "
                    f"which the anvi'o metagenomics workflow cannot process yet. You will have to leave "
                    f"{'them' if read_types[sra.SINGLE_END_SHORT_READS] > 1 else 'it'} out of your samples-txt.")


def get_accessions(args):
    """Collect accessions from whichever of the two possible inputs was given."""

    if not args.accession_list and not args.samples_txt:
        raise ConfigError("You need to tell anvi'o which SRA runs you are interested in, either with a file of "
                          "accessions (`--accession-list`) or with a samples-txt that has an `sra_accession` "
                          "column in it (`--samples-txt`).")

    if args.accession_list and args.samples_txt:
        raise ConfigError("Please use either `--accession-list` or `--samples-txt`, but not both at once.")

    if args.accession_list:
        filesnpaths.is_file_exists(args.accession_list)
        accessions = [line.strip() for line in open(args.accession_list) if line.strip()]
    else:
        from anvio.artifacts.samples_txt import SamplesTxt
        samples_txt = SamplesTxt(args.samples_txt, expected_format="free")
        accessions = samples_txt.all_sra_accessions()

        if not accessions:
            raise ConfigError(f"There is not a single SRA accession in '{args.samples_txt}'. This program needs a "
                              f"samples-txt with an `sra_accession` column in it.")

    if not accessions:
        raise ConfigError(f"Anvi'o found no accessions in '{args.accession_list}'.")

    return list(dict.fromkeys(accessions))


def describe_where_accessions_came_from(args):
    return f"your accession list ('{args.accession_list}')" if args.accession_list else f"'{args.samples_txt}'"


def get_args():
    from anvio.argparse import ArgumentParser
    parser = ArgumentParser(description=__description__)

    groupI = parser.add_argument_group('INPUT', "Where the accessions come from. Pick one.")
    groupI.add_argument('--accession-list', type=str, default=None, help="A file with one SRA run accession per line.")
    groupI.add_argument('--samples-txt', type=str, default=None, help="A samples-txt file with an `sra_accession` column.")

    groupO = parser.add_argument_group('OUTPUT')
    groupO.add_argument('-o', '--output-file', type=str, default='SRA-METADATA.txt',
                        help="Where to keep what anvi'o learns. If this file already exists, anvi'o only looks up the "
                             "accessions that are missing from it and adds them, so anything you have corrected by "
                             "hand stays as you left it. The default is '%(default)s', which is also where the "
                             "metagenomics workflow looks unless you tell it otherwise.")

    return parser.get_args(parser)


if __name__ == '__main__':
    main()

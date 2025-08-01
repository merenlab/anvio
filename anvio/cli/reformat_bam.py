#!/usr/bin/env python
# -*- coding: utf-8

import os
import sys
import pysam


import anvio
import anvio.utils as utils
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths

from anvio.errors import ConfigError, FilesNPathsError

__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__   = []
__license__   = "GPL 3.0"
__version__   = anvio.__version__
__authors__   = ['telatin']
__requires__  = ["bam-file", "contig-rename-report-txt"]
__provides__  = ["bam-file"]
__description__ =  ("Reformat a BAM file to match the updated sequence names after running "
                    "anvi-script-reformat-fasta. You will need this script to fix your BAM "
                    "file if you run `anvi-script-reformat-fasta` on a FASTA file of "
                    "sequences *after* you already used the previous version of the FASTA "
                    "file for read recruitment.")

OUTPUT_FILE_SUFFIX = ".reformatted.bam"


def main():
    try:
        run_program()
    except ConfigError as e:
        print(e)
        sys.exit(-1)
    except FilesNPathsError as e:
        print(e)
        sys.exit(-2)


def load_dict_from_tsv(path):
    name_dict = { }
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            if not '\t' in line:
                raise Exception(f"ERROR: Invalid line in {path} (missing tab): {line}")
            new, old = line.strip().split("\t")
            utils.is_this_name_OK_for_database('contig name prefix', new)
            old = old.split(" ")[0].split("\t")[0]
            name_dict[old] = new

    return name_dict


def update_read(read, new_header):
    new_read = pysam.AlignedSegment(header=new_header)
    new_read.query_name = read.query_name
    new_read.query_sequence = read.query_sequence
    new_read.flag = read.flag
    new_read.reference_id = read.reference_id
    new_read.reference_start = read.reference_start
    new_read.mapping_quality = read.mapping_quality
    new_read.cigar = read.cigar
    new_read.next_reference_id = read.next_reference_id
    new_read.next_reference_start = read.next_reference_start
    new_read.template_length = read.template_length
    new_read.query_qualities = read.query_qualities
    new_read.tags = read.tags
    return new_read


def run_program():
    args = get_args()
    run = terminal.Run()

    # print dump of `args`
    run.info('Input BAM file', args.bam_file)
    run.info('Rename list', args.report_file)

    # Parameters check
    input_basename = os.path.basename(args.bam_file)
    output_filename = args.output_bam

    if not args.output_bam:
        # Remove .bam and add .reformatted.bam
        output_filename = os.path.join(os.path.dirname(args.bam_file),
                                       input_basename[:-4] + OUTPUT_FILE_SUFFIX if input_basename.endswith(".bam") else input_basename + ".reformatted.bam"
        )
        run.info("Output BAM file", output_filename)

    # Some sanity checks
    filesnpaths.is_file_bam_file(args.bam_file)
    filesnpaths.is_file_exists(args.report_file)

    if args.bam_file == output_filename:
        raise ConfigError("The input file and the output file names can't be the same :/")

    filesnpaths.is_output_file_writable(output_filename, ok_if_exists=False)

    # Load the list (tsv file, two columns, new name TAB old name)
    run.warning(None, header='WHAT WAS THERE', lc="cyan")
    rename_dict = load_dict_from_tsv(args.report_file)
    run.info(f"Sequences in {args.report_file}", f"{len(rename_dict):,}", mc="green")

    # Open the BAM file with pysam and rename all the contigs in header and alignment records
    with pysam.AlignmentFile(args.bam_file, "rb") as bam:
        with pysam.AlignmentFile(args.bam_file, "rb") as bam:
            # Convert header to dictionary
            header_dict = bam.header.to_dict()

            # sanity check for matching number of sequences
            run.info("Sequences in BAM file", len(header_dict['SQ']))
            if len(rename_dict) != len(header_dict['SQ']):
                all_seqs_in_bam = [x['SN'] for x in header_dict['SQ']]
                mismatching_bam = [x for x in all_seqs_in_bam if x not in rename_dict]
                mismatching_rename = [x for x in rename_dict if x not in all_seqs_in_bam]
                raise ConfigError(f"ERROR: The number of sequences in the BAM file ({len(header_dict['SQ'])}) does "
                                  f"not match to the number of sequences in the reformat report ({len(rename_dict)}). "
                                  f"For example, here are the sequences that are only in the BAM file and not in the "
                                  f"reformat report (if any): {', '.join(mismatching_bam)} . And here are the sequences that "
                                  f"are only in the reformat report and not in the BAM file (if any):  {', '.join(mismatching_rename)}")

            # update contig names, and create new header
            for sq in header_dict['SQ']:
                old_name = sq['SN']
                if old_name in rename_dict:
                    sq['SN'] = rename_dict[old_name]
            new_header = pysam.AlignmentHeader.from_dict(header_dict)

        # Create a new BAM file with the new header and updated records
        with pysam.AlignmentFile(args.bam_file, "rb") as bam, pysam.AlignmentFile(output_filename, "wb", header=new_header) as out_bam:
            for read in bam:
                updated_read = update_read(read, new_header)
                out_bam.write(updated_read)

    run.warning(None, header='WHAT WAS DONE', lc="cyan")
    run.info("Output BAM file created", output_filename)


def get_args():
    from anvio.argparse import ArgumentParser

    parser = ArgumentParser(description=__description__)

    groupA = parser.add_argument_group('INPUT', 'The input file you wish to work with.')
    groupA.add_argument('bam_file', metavar='BAM FILE')
    groupA.add_argument('-r', '--report-file', metavar="TXT FILE", help="A TAB-delimited file produced by `anvi-script-reformat-fasta` "
                        "to report sequence name changes due to `--simplify-names` flag.", required=True)

    groupB = parser.add_argument_group('OUTPUT', 'Dealing with the output.')
    groupB.add_argument('-o', '--output-bam', required=False, metavar='FORMATTED BAM',
                        help=f"The output BAM file name. If you don't specify an output file name, the new BAM file will be reported "
                        f"with the suffix `{OUTPUT_FILE_SUFFIX}`.")


    return parser.get_args(parser)


if __name__ == '__main__':
    main()

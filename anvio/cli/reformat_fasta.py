#!/usr/bin/env python
# -*- coding: utf-8

import sys
import shutil

import pandas as pd

import anvio
import anvio.fastalib as u
import anvio.utils as utils
import anvio.terminal as terminal
import anvio.constants as constants
import anvio.filesnpaths as filesnpaths

from anvio.errors import ConfigError, FilesNPathsError

__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['meren', 'ekiefl', 'vinisalazar']
__requires__ = ["fasta"]
__provides__ = ["contigs-fasta", "contig-rename-report-txt"]
__description__ =  ("Reformat FASTA file (remove contigs based on length, or based on a given list of "
                    "deflines, and/or generate an output with simpler names)")


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
    pp = terminal.pretty_print
    P = terminal.pluralize

    if args.keep_ids and args.exclude_ids:
        raise ConfigError("You can't use`--exclude-ids and --keep-ids together :/")

    if args.exact_length and args.min_len:
        raise ConfigError(f"You can't ask your reads to be an exact lenght of {args.exact_length} "
                          f"and longer than {args.min_len} at the same time. It just doesn't make "
                          f"any sense and that's why.")

    if args.output_file == args.contigs_fasta:
        raise ConfigError("You can't set the same path for both your input file name and output "
                          "file name. It makes no sense UNLESS you wish to 'overwrite' your input "
                          "file in-place. But anvi'o not only has a flag for that (`--overwrite-input`) "
                          "but also has means to track those who do not read the help menu for "
                          "programs. Just so you know.")

    if args.output_file and args.overwrite_input:
        raise ConfigError("You can't ask anvi'o to overwrite your input file, and also provide an "
                          "output file name at the same time.")

    if not args.overwrite_input and not args.output_file:
        raise ConfigError("You have not provided an output file name. If you are feeling extra "
                          "adventurous today and would like anvi'o to overwrite your input file "
                          "you must use the `--overwrite` flag.")

    # if the user wants to overwrite the input file, we will do it by first
    # setting the output file path to a temp file, and then moving it on
    # top of the input file, overwriting it forever.
    if args.overwrite_input:
       args.output_file = filesnpaths.get_temp_file_path()
       if args.contigs_fasta.endswith(".gz"):
            args.output_file += ".gz"

    filesnpaths.is_file_fasta_formatted(args.contigs_fasta)
    filesnpaths.is_output_file_writable(args.output_file)

    if not args.ignore_empty_sequences:
        # check for empty sequences
        fasta = u.SequenceSource(args.contigs_fasta)

        sequence_lenghts = []
        while next(fasta):
            sequence_lenghts.append(len(fasta.seq))

        num_empty_sequences = sequence_lenghts.count(0)
        if num_empty_sequences:
            raise ConfigError(f"We have a problem, Houston. Of the total {len(sequence_lenghts)} entries "
                              f"in your FASTA file, {num_empty_sequences} {P('has', num_empty_sequences, alt='have' )} "
                              f"no sequences (i.e., they're blank). You have two options: either (1) use the "
                              f"flag `--ignore-empty-seqeunces` so anvi'o can ignore these FASTA entries, or "
                              f"(2) go back to your FASTA file and figure out why they are empty.")
    else:
        pass

    report_file = open(args.report_file, 'w') if args.report_file and args.simplify_names else None
    prefix = args.prefix if args.prefix else None

    if prefix:
        utils.is_this_name_OK_for_database('contig name prefix', prefix)

    if args.exclude_ids:
        filesnpaths.is_file_exists(args.exclude_ids)
        exclude_ids = set([l.split('\t')[0].strip() for l in open(args.exclude_ids, 'r').readlines()])
        run.info('Input IDs to remove', '%d found' % len(exclude_ids))
    else:
        exclude_ids = None

    if args.keep_ids:
        filesnpaths.is_file_exists(args.keep_ids)
        keep_ids = set([l.split('\t')[0].strip() for l in open(args.keep_ids, 'r').readlines()])
        run.info('Input IDs to consider', '%d found' % len(keep_ids))
    else:
        keep_ids = None

    if args.seq_type is not None:
        replace_chars = True
        if args.seq_type == 'AA':
            acceptable_chars = set(sorted(list(constants.AA_to_single_letter_code.values())))
            acceptable_chars.add('X')
            replacement = 'X'
        else:
            acceptable_chars = set(constants.nucleotides)
            acceptable_chars.add('N')
            replacement = 'N'
    else:
        replace_chars = False

    if args.export_gap_counts_table:
        gaps_info_list = []

    # summary of where we are
    run.info('Input', args.contigs_fasta)

    if args.overwrite_input:
        run.info('Output', "(anvi'o will overwrite your input file)", mc='red')
    else:
        run.info('Output', args.output_file)

    total_num_nucleotides = 0
    total_num_contigs = 0
    total_num_nucleotides_removed = 0
    total_num_nucleotides_modified = 0
    total_num_contigs_removed = 0

    fasta = u.SequenceSource(args.contigs_fasta)
    output = u.FastaOutput(args.output_file)

    while next(fasta):
        l = len(fasta.seq)

        total_num_nucleotides += l
        total_num_contigs += 1

        if not l:
            # this is an entry with a blank sequence
            total_num_contigs_removed += 1
            continue

        if replace_chars:
            seq = []
            for char in fasta.seq:
                if char not in acceptable_chars:
                    seq.append(replacement)
                    total_num_nucleotides_modified += 1
                else:
                    seq.append(char)
            fasta.seq = ''.join(seq)

        if keep_ids and fasta.id.split()[0] not in keep_ids:
            total_num_nucleotides_removed += l
            total_num_contigs_removed += 1
            continue

        if exclude_ids and fasta.id.split()[0] in exclude_ids:
            total_num_nucleotides_removed += l
            total_num_contigs_removed += 1
            continue

        length = len(fasta.seq)
        if args.exact_length and length != args.exact_length:
            total_num_nucleotides_removed += l
            total_num_contigs_removed += 1
            continue
        elif length < args.min_len:
            total_num_nucleotides_removed += l
            total_num_contigs_removed += 1
            continue

        num_gaps = fasta.seq.count('-')
        if args.export_gap_counts_table:
            gaps_info_list.append([fasta.id, num_gaps])

        percentage_of_gaps = num_gaps * 100.0 / l
        if percentage_of_gaps >= args.max_percentage_gaps:
            total_num_nucleotides_removed += l
            total_num_contigs_removed += 1
            continue

        if num_gaps >= args.max_gaps:
            total_num_nucleotides_removed += l
            total_num_contigs_removed += 1
            continue

        if args.simplify_names or prefix:
            if prefix:
                defline = '%s_%012d' % (prefix, fasta.pos)
            else:
                defline = 'c_%012d' % fasta.pos

            output.write_id(defline)
            output.write_seq(fasta.seq, split = False)

            if report_file and args.simplify_names:
                report_file.write('%s\t%s\n' % (defline, fasta.id))
        else:
            output.store(fasta, split = False)

    if args.export_gap_counts_table:
        df = pd.DataFrame(gaps_info_list, columns=["header", "num_gaps"])
        df.to_csv(args.export_gap_counts_table + ".tsv", sep='\t', index=False)

    if report_file:
        report_file.close()

    fasta.close()
    output.close()

    if args.overwrite_input:
        shutil.move(args.output_file, args.contigs_fasta)

    run.warning(None, header='WHAT WAS THERE', lc="cyan")
    run.info('Total num contigs', total_num_contigs)
    run.info('Total num nucleotides', total_num_nucleotides)

    run.warning(None, header='WHAT WAS ASKED', lc="cyan")
    run.info('Simplify deflines?', "Yes" if args.simplify_names else "No")
    run.info('Add prefix to sequence names?', f"Yes, add '{args.prefix}'" if args.prefix else "No")
    run.info('Exact length of contigs to keep', args.exact_length, mc="red") if args.exact_length else run.info("Minimum length of contigs to keep", args.min_len)
    run.info('Max % gaps allowed', '%.2f%%' % args.max_percentage_gaps)
    run.info('Max num gaps allowed', args.max_gaps)
    run.info('Exclude specific sequences?', f"Yes, those listed in {args.exclude_ids}" if args.exclude_ids else "No")
    run.info('Keep specific sequences?', f"Yes, those listed in {args.keep_ids}" if args.keep_ids else "No")
    run.info('Enforce sequence type?', f"Yes, enforce '{args.seq_type}'" if args.seq_type else "No")

    run.warning(None, header='WHAT HAPPENED', lc="cyan")
    if args.ignore_empty_sequences:
        run.info('Entries w/blank sequences discarded', 'Yes', mc='red')
    run.info('Contigs removed', f'{pp(total_num_contigs_removed)} ({total_num_contigs_removed * 100.0 / total_num_contigs:.2f}% of all)', mc='green')
    run.info('Nucleotides removed', f'{pp(total_num_nucleotides_removed)} ({total_num_nucleotides_removed * 100.0 / total_num_nucleotides:.2f}% of all)', mc='green')
    run.info('Nucleotides modified', f'{pp(total_num_nucleotides_modified)} ({total_num_nucleotides_modified * 100.0 / total_num_nucleotides:.5f}% of all)', mc='green')
    run.info('Deflines simplified', args.simplify_names, mc='green')

    if args.overwrite_input:
        run.info_single("The contents of your input file have changed because you used "
                        "the flag `--overwrite-input`.", nl_before=1, nl_after=1)


def get_args():
    from anvio.argparse import ArgumentParser
    parser = ArgumentParser(description=__description__)

    groupA = parser.add_argument_group('INPUT', 'The input file you wish to work with.')
    groupA.add_argument('contigs_fasta', metavar='FASTA FILE')

    groupB = parser.add_argument_group('OUTPUT', 'Dealing with the output.')
    groupB.add_argument('-o', '--output-file', required=False, metavar='FASTA FILE PATH',
                        help="Output file path for a NEW FASTA file. If you want all changes to be "
                             "stored in the same file, see the `--overwrite` flag.")
    groupB.add_argument('--overwrite-input', default=False, action="store_true",
                        help="Instead of reporting a new file, overwrite the input file. For those of "
                             "us who loves to live life at the very edge of it. A case of we can't go "
                             "back, Kate.")
    groupB.add_argument('-r', '--report-file', required=False, metavar='TXT FILE PATH',
                        help="When you run this program with `--simplify-names` flag, all changes to "
                             "deflines will be reported in this file in case you need to go back to this "
                             "information later. It is not mandatory to declare one, but it is a very good "
                             "idea to have it.")

    groupC = parser.add_argument_group('OPERATIONS: RENAME YOUR DEFLINES')
    groupC.add_argument('--simplify-names', default=False, action="store_true",
                        help="Edit deflines to make sure they contigs have simple names.")
    groupC.add_argument('--prefix', default=None, metavar="PREFIX",
                        help="Use this parameter if you would like to add a prefix to your contig\
                              names while simplifying them. The prefix must be a single word (you\
                              can use underscor character, but nothing more!).")

    groupD = parser.add_argument_group('OPERATIONS: TAME YOUR SEQUENCES')
    groupD.add_argument('-l', '--min-len', type=int, default=0, metavar='INT',
                        help="Minimum length of contigs to keep (contigs shorter than this value\
                              will not be included in the output file). The default is %(default)d,\
                              so nothing will be removed if you do not declare a minimum size.")
    groupD.add_argument('--exact-length', type=int, default=None, metavar='INT',
                        help="Exact lenght of the sequences you wish to keep from a FASTA fille. Why \
                              would anyone need that? Anvi'o does not know. But there you have it. This \
                              parameter is indeed incompatible with a lot of others, such as `--min-len`.")
    groupD.add_argument('--max-percentage-gaps', type=float, default=100, metavar='PERCENTAGE',
                        help="Maximum fraction of gaps in a sequence (any sequence with \
                              more gaps will be removed from the output FASTA file). The \
                              default is %(default)f.")
    groupD.add_argument('-M','--max-gaps', type=float, default=1000000,
                        help="Maximum amount of gaps allowed per sequence in the alignment. \
                              Don't know which threshold to pick? Use --export-gap-counts-table \
                              to explore the gap counts per sequence distribution!")
    groupD.add_argument('-i', '--exclude-ids', required=False, metavar='TXT FILE',
                        help="IDs to remove from the FASTA file. You cannot provide both\
                              --keep-ids and --exclude-ids.")
    groupD.add_argument('--export-gap-counts-table', required=False, metavar='TSV FILE',
                        help="Export a table with the number of gaps per sequence. \
                              Please provide a prefix to name the file.")
    groupD.add_argument('-I', '--keep-ids', type=str, required=False, metavar='TXT FILE',
                        help="If provided, all IDs not in this file will be excluded from the\
                              reformatted FASTA file. Any additional filters (such as --min-len)\
                              will still be applied to the IDs in this file. You cannot provide both\
                              --exclude-ids and --keep-ids.")
    groupD.add_argument('--seq-type', default=None, metavar="SEQ TYPE", choices={'AA', 'NT'},
                        help=("Supply either 'NT' or 'AA' (if you want). If 'NT', any characters besides {A,C,T,G} will "
                              "by replaced with 'N'. If 'AA', any characters that are not 1-letter amino acid "
                              "characters will be replaced with 'X'. If you don't supply anything, no charaters will be "
                              "modified."))
    groupD.add_argument('--ignore-empty-sequences', default=False, action="store_true",
                        help=("If your FASTA file contains entries with no sequences, you will either have to ask anvi'o "
                              "to ignore them, or go back to your file and figure out why it is the case (because you will "
                              "certainly get an error somewhere during the process)."))
    return parser.get_args(parser)


if __name__ == '__main__':
    main()

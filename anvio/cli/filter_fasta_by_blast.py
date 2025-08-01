#!/usr/bin/env python
# -*- coding: utf-8

import sys
import pandas as pd

# maybe not all are needed, who is to say?
import anvio
import anvio.fastalib as u
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths

from anvio.errors import ConfigError, FilesNPathsError

__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['ekiefl', 'meren', 'blankenberg']
__provides__ = ["contigs-fasta"]
__requires__ = ["contigs-fasta", "blast-table"]
__description__ = "Filter FASTA file according to BLAST table (remove sequences with bad BLAST alignment)"



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
    pp = terminal.pretty_print

    A = lambda x: args.__dict__[x] if x in args.__dict__ else None
    input_fasta_file_path = A('fasta_file')
    blast_output_file_path = A('blast_output')
    output_file_path = A('output_file')
    threshold = A('threshold') or 30.0
    outfmt = A('outfmt')
    ok_if_exists = A('just_do_it')

    if not output_file_path:
        raise ConfigError("An output file path is required for this to work.")

    if not blast_output_file_path:
        raise ConfigError("You forgot to provide an input BLAST output.")

    if not input_fasta_file_path:
        raise ConfigError("You must declare an input FASTA file path...")

    try:
        threshold = float(threshold)
    except:
        raise ConfigError("Threshold must be of type int or float.")

    filesnpaths.is_output_file_writable(output_file_path, ok_if_exists=ok_if_exists)

    cols = outfmt.split(" ")
    required_cols = ['qseqid', 'bitscore', 'length', 'qlen', 'pident']
    if len([c for c in required_cols if c not in cols]):
        raise ConfigError("According to your `outfmt` parameter, columns in your BLAST output table "
                          "lack one or more important stuff. These columns must be present in your "
                          "output table (see help for details): %s." % ', '.join(required_cols))

    filesnpaths.is_file_tab_delimited(blast_output_file_path, expected_number_of_fields=len(cols))

    # DONE WITH THE INITAL CHECKS. TIME TO GET READY FOR THE ACTUAL STUFF
    # BUT FIRST, REPORTING.

    run.info('Original FASTA file', input_fasta_file_path)
    run.info('BLAST output file', blast_output_file_path)
    run.info('Threshold for filtering', '%.1f%%' % threshold)
    run.info('OUTFMT', "'%s'" % ' '.join(cols))

    progress.new('Processing data')

    progress.update('Reading BLAST output file...')
    blastout = pd.read_csv(blast_output_file_path, sep="\t", names=cols)

    progress.update('Selecting the hits with highest bitscore for each sequence...')
    idx = blastout.groupby("qseqid")["bitscore"].transform(max) == blastout["bitscore"]
    blastout = blastout[idx]

    progress.update('Generating a proper idnetity column in the data frame...')
    # Rarely, the reported align length is greater
    # than the query sequence length as a result of indels. If this ever
    # happens we set the align length = query sequence length.
    blastout["proper_pident"] = blastout["length"]/blastout["qlen"]
    blastout.loc[blastout["proper_pident"]>1, "proper_pident"] = 1
    blastout["proper_pident"] *= blastout["pident"]

    progress.update('Filtering hits based on proper identity scores...')
    blastout = blastout[blastout["proper_pident"] >= args.threshold]
    blastout = blastout.sort_values("proper_pident",ascending=False)
    seqs_to_keep = [str(x) for x in list(blastout.qseqid.unique())]

    # iterate through fasta. If fasta defline is in seqs_to_keep, add to new fasta
    progress.update('Filtering input FASTA file...')
    output = u.FastaOutput(output_file_path)
    fasta = u.SequenceSource(args.fasta_file)
    total_num_seq, filtered_num_seq = 0, 0
    while next(fasta):
        total_num_seq += 1
        if fasta.id in seqs_to_keep:
            filtered_num_seq += 1
            output.write_id(fasta.id)
            output.write_seq(fasta.seq, split = False)

    fasta.close()
    output.close()

    progress.end()

    run.info('Total number of sequences', pp(total_num_seq))
    run.info('Final number of filtered sequences', pp(filtered_num_seq))
    run.info('Filtered FASTA output', output_file_path, mc='green')


def get_args():
    from anvio.argparse import ArgumentParser

    parser = ArgumentParser(description=__description__)

    parser.add_argument(*anvio.A('fasta-file'), **anvio.K('fasta-file', {'required': False}))
    parser.add_argument(*anvio.A('output-file'), **anvio.K('output-file', {'required': False}))
    parser.add_argument("-b", "--blast-output", required = True, help = "BLAST table generated with \
                                blastp. `--outfmt 6` as the output format is assumed.", metavar='TAB DELIMITED FILE')
    parser.add_argument("-s", "--outfmt", required = True, help = "Specify the column ordering of your\
                                BLAST report. We add the following paramter to our BLAST searches\
                                so the output report contains the `qlen` field, which is not included by\
                                default: `-outfmt '6 qseqid sseqid pident length mismatch gapopen qstart\
                                qend sstart send evalue bitscore qlen slen'`. You may have used a different\
                                `-outfmt` paramter, and you should use this parameter to explicitly define the\
                                column names in your output file. For instance, if you had used the parameter\
                                mentioned above, then the correct version of this parameter would be: \"qseqid\
                                sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore\
                                qlen slen\". Regardless of the BLAST output format, your columns MUST contain\
                                the following parameters for this program to work properly: 'qseqid', \
                                'bitscore', 'length', 'qlen', and 'pident'.")
    parser.add_argument("-t", "--threshold", required = True, type = float, default=30.0, help = "What `proper_pident`\
                                threshold do you want to use for filtering out sequences whose top bit-score\
                                matches have `proper_pident`s less than this threshold?  We have defined \
                                `proper_pident` to be the percentage of the query amino acids that both \
                                aligned to and were identical to the corresponding matched amino acid. \
                                Note that the `pident` parameter output by BLAST does not include regions\
                                of the query sequence unaligned to the matched sequence, whereas `proper_pident`\
                                does. For example, a sequence that's only half aligned by a match but with 100%%\
                                identity at matched regions has a `pident` of 100 but a `proper_pident` of 50. The\
                                default is %(default).1f%%.")
    parser.add_argument(*anvio.A('just-do-it'), **anvio.K('just-do-it'))

    return parser.get_args(parser)


if __name__ == '__main__':
    main()

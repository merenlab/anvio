# -*- coding: utf-8
# pylint: disable=line-too-long

"""Pangenome graph preprocessing class"""

import os
import sys
import subprocess
import pandas as pd

from Bio import SearchIO, SeqIO
from Bio.SeqRecord import SeqRecord

import anvio
import anvio.fastalib as f
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths

from anvio.argparse import ArgumentParser
from anvio.errors import ConfigError, FilesNPathsError


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "A. Murat Eren"
__email__ = "a.murat.eren@gmail.com"
__status__ = "Development"


run = terminal.Run()
progress = terminal.Progress()
pp = terminal.pretty_print

class external_genomes_preprocess:
    """Takes in a fasta.txt file
    """

    def __init__(self, args, run=run, progress=progress):
        self.args = args
        self.run = run
        self.progress = progress

        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        self.fasta_text_file = A('fasta_text_file')
        self.output_dir = A('output_dir')

        if A('prioritize_genome_size'):
            self.prioritize = 'genome_size'
        elif A('prioritize_number_of_contigs'):
            self.prioritize = 'number_of_contigs'
        else:
            raise ConfigError("Please pick a method :)")

        self.max_contig_num = 100
        self.do_reorient_contigs = True
        self.max_ani = 0.95
        self.max_genomes = 5
        self.fasta_text_df = pd.read_csv(self.fasta_text_file, header=0, sep='\t', index_col='name')


    def process(self):
        filesnpaths.gen_output_directory(self.output_dir)
        self.sanity_check()
        self.preprocess_fasta()


    def sanity_check(self):
        filesnpaths.is_output_dir_writable(self.output_dir)


    def get_num_contigs_and_genome_length(self, file_path):
        fasta = f.SequenceSource(file_path)

        num_contigs = 0
        length = 0

        while next(fasta):
            num_contigs += 1
            length += len(fasta.seq)

        return(num_contigs, length)


    def preprocess_fasta(self, reference_fasta = None):

        self.run.warning(None, header='Reorient contigs.', lc='green')

        best_length = 0
        best_num_contigs = sys.maxsize
        genome_properties = {}

        for name, row in self.fasta_text_df.iterrows():
            fasta_file = row['path']
            num_contigs, length = self.get_num_contigs_and_genome_length(fasta_file)

            if self.max_contig_num != -1 and num_contigs > self.max_contig_num:
                self.run.info_single(f'Removed {name} due to contig number filter.')
                self.fasta_text_df.drop(name, inplace=True)
            else:
                self.run.info_single(f"Inlcude {name} with {num_contigs} contigs and {length} nt length.)")

                if self.prioritize == 'genome_size' and length > best_length:
                    reference_fasta = name
                    best_length = length
                elif self.prioritize == 'number_of_contigs' and num_contigs < best_num_contigs:
                    reference_fasta = name
                    best_num_contigs = num_contigs
                elif self.prioritize == 'number_of_contigs' and num_contigs == best_num_contigs:
                    if length > best_length:
                        reference_fasta = name
                        best_length = length
                else:
                    pass

                genome_properties[name] = {'path': fasta_file, 'num_contigs': num_contigs, 'length': length}

        if self.do_reorient_contigs and genome_properties:

            self.run.info_single(f"Method to select reference {'Minimum number of contigs' if self.prioritize == 'number_of_contigs' else 'Maximum length'}.")
            self.run.info_single(f"Reference FASTA is {reference_fasta}.")

            if genome_properties[reference_fasta]['num_contigs'] != 1:
                if self.prioritize == 'genome_size':
                    raise ConfigError("Sadly, the reference genome anvi'o chose for you based on geonome size priority seems to have multiple contigs, "
                                      "which suggests that it is not a complete genome :/ The current implementation of this script does not know "
                                      "how to use a genome with multiple contigs as reference. You can try to re-run the program with `--prioritize-number-of-contigs` "
                                      "flag if you have a complete genome.")
                elif self.prioritize == 'number_of_contigs':
                    raise ConfigError("The genome in your collection with the smallest number of contigs is still not a complete genome :( There is nothing "
                                      "this program can do for you at this point. Sorry!")
                else:
                    self.run.info_single('WTF case')

            for name, fasta_dict in genome_properties.items():

                fasta_file = fasta_dict['path']
                num_contigs = fasta_dict['num_contigs']
                length = fasta_dict['length']
                contigs = SeqIO.to_dict(SeqIO.parse(fasta_file, 'fasta'))

                num_hits = 0
                num_reverse = 0

                blast_result_file = os.path.join(self.output_dir, name + '-blast.xml')

                if name != reference_fasta:
                    cmd_line = ["blastn", "-query", genome_properties[reference_fasta]['path'], "-out", blast_result_file, "-outfmt", "5", "-subject", fasta_file]
                    subprocess.run(cmd_line, text=True, capture_output=True)

                    qresult = SearchIO.read(blast_result_file, 'blast-xml')
                    contig_hits = []

                    for hit in qresult.hits:
                        best_hsp = max(hit.hsps, key=lambda hsp: hsp.bitscore)
                        contig_hits.append({
                            'id': hit.id,
                            'query_start': best_hsp.query_start,
                            'query_end': best_hsp.query_end,
                            'hit_start': best_hsp.hit_start,
                            'hit_end': best_hsp.hit_end,
                            'strand': best_hsp.hit_strand,  # 1 or -1
                            'bitscore': best_hsp.bitscore
                        })

                    # Sort by query_start to get correct order
                    contig_hits.sort(key=lambda x: x['query_start'])

                    # Prepare oriented contigs
                    oriented_contigs = []
                    for ch in contig_hits:
                        seq = contigs[ch['id']].seq
                        if ch['strand'] == -1:
                            seq = seq.reverse_complement()
                            num_reverse += 1
                        else:
                            num_hits += 1

                        new_record = SeqRecord(
                            seq=seq,
                            id=ch['id'],
                            description=f"aligned to query {ch['query_start']}-{ch['query_end']} (strand: {ch['strand']})"
                        )
                        oriented_contigs.append(new_record)

                    self.run.info_single(f"Found {num_hits + num_reverse} / Reversed {num_reverse} out of {num_contigs} contigs for fasta {fasta_file}.")

                    output_file = os.path.abspath(os.path.join(self.output_dir, os.path.basename(fasta_file)))
                    self.fasta_text_df.at[name, 'path'] = output_file

                    with open(output_file, "w") as output_handle:
                        SeqIO.write(oriented_contigs, output_handle, 'fasta')

                    self.run.info_single(f"Currated copy of {fasta_file} created as {output_file}.")

                    if num_hits + num_reverse != num_contigs:
                        self.run.info_single("Warning one or more contigs in the file might be created from contamination.")

                    self.fasta_text_df.to_csv(self.fasta_text_file, sep='\t')

            self.run.info_single("Done.")
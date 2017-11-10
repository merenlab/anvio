# -*- coding: utf-8
# pylint: disable=line-too-long
"""
    Classes for gene calling.
"""

import shutil

import anvio
import anvio.utils as utils
import anvio.tables as t
import anvio.terminal as terminal
import anvio.fastalib as u
import anvio.filesnpaths as filesnpaths

from anvio.errors import ConfigError

from anvio.drivers.prodigal import Prodigal


__author__ = "A. Murat Eren"
__copyright__ = "Copyright 2016, The anvio Project"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "A. Murat Eren"
__email__ = "a.murat.eren@gmail.com"


run = terminal.Run()
progress = terminal.Progress()


class GeneCaller(object):
    def __init__(self, fasta_file_path, gene_caller=None, progress=progress, run=run, debug=False):
        filesnpaths.is_file_exists(fasta_file_path)
        filesnpaths.is_file_fasta_formatted(fasta_file_path)

        self.fasta_file_path = fasta_file_path

        self.run = run
        self.progress = progress

        self.debug = debug
        self.tmp_dirs = []

        self.gene_callers = {'prodigal': Prodigal}

        self.gene_caller = gene_caller or 'prodigal'

        if self.gene_caller not in self.gene_callers:
            raise ConfigError("The gene caller you requested ('%s') is not available at this point.\
                                here is a list of what we have: %s." % (self.gene_caller, ', '.join(self.gene_callers)))


    def process(self):
        output_dir = filesnpaths.get_temp_directory_path()
        self.tmp_dirs.append(output_dir)
        gene_caller = self.gene_callers[self.gene_caller]()

        gene_calls_dict, protein_sequences_dict = gene_caller.process(self.fasta_file_path, output_dir)

        if not self.debug:
            self.clean_tmp_dirs()

        return gene_calls_dict, protein_sequences_dict


    def clean_tmp_dirs(self):
        for tmp_dir in self.tmp_dirs:
            shutil.rmtree(tmp_dir)


class ExternalCeneCalls(object):
    def __init__(self, fasta_file_path, gene_calls_file, progress=progress, run=run, debug=False, ignore_internal_stop_codons=False):
        self.fasta_file_path = fasta_file_path
        self.gene_calls_file = gene_calls_file
        self.progress = progress
        self.run = run


    def process(self):
        gene_calls_dict = utils.get_TAB_delimited_file_as_dictionary(gene_calls_file,
                                                                     expected_fields=t.genes_in_contigs_table_structure,
                                                                     only_expected_fields=True,
                                                                     column_mapping=[int, str, int, int, str, int, str, str])


        if not len(gene_calls_dict):
            raise ConfigError("You provided an external gene calls file, but it returned zero gene calls. Assuming that\
                               this is an error, anvi'o will stop here and complain. If this is not an error and you\
                               in fact expected this, the proper way of doing this is to use `--skip-gene-calls` flag,\
                               instead of providing an emtpy external gene calls file. You don't agree? You need this\
                               for some weird step for you weird pipeline? Let us know, and we will consider changing\
                               this.")

        self.run.info("External gene calls", "%d gene calls recovered and will be processed." % len(gene_calls_dict))

        protein_sequences = {}
        contig_sequences = {}
        
        fasta = u.SequenceSource(self.fasta_file_path)
        while next(fasta):
            contig_sequences[fasta.id] = {'sequence': fasta.seq}
        fasta.close()

        num_genes_with_internal_stops = 0
        number_of_impartial_gene_calls = 0
        for gene_callers_id in gene_calls_dict:
            gene_call = gene_calls_dict[gene_callers_id]
            contig_name = gene_call['contig']

            if contig_name not in contig_sequences:
                # remove the partial contigs database so things don't get screwed later
                os.remove(self.db_path)
                raise ConfigError("You are in big trouble :( The contig name '%s' in your external gene callers file\
                                    does not appear to be in the contigs FASTA file. How did this happen?" % contig_name)

            if gene_call['partial']:
                protein_sequences[gene_callers_id] = ''
                number_of_impartial_gene_calls += 1
                continue

            sequence = contig_sequences[contig_name]['sequence'][gene_call['start']:gene_call['stop']]
            if gene_call['direction'] == 'r':
                sequence = utils.rev_comp(sequence)

            protein_sequence = utils.get_DNA_sequence_translated(sequence, gene_callers_id)

            # check if there are any internal stops:
            if protein_sequence.find('*') > -1:
                if ignore_internal_stop_codons:
                    protein_sequence = protein_sequence.replace('*', 'X')
                    num_genes_with_internal_stops += 1
                else:
                    os.remove(self.db_path)
                    raise ConfigError("Oops. Anvi'o run into an amino acid seqeunce (that corresponds to the gene callers id '%s')\
                                       which had an internal stop codon :/ This usually indicates that your external gene calls\
                                       have problems. If you still want to continue, you can ask anvi'o to ignore internal stop\
                                       codons on your own risk. It will probably look very ugly on your screen, but here is the\
                                       DNA sequence for that gene in case you don't trust anvi'o (which only would be fair since\
                                       anvi'o does not trust you either): %s" % (str(gene_callers_id), sequence))

            protein_sequences[gene_callers_id] = protein_sequence

        if num_genes_with_internal_stops:
            percent_genes_with_internal_stops = num_genes_with_internal_stops * 100.0 / len(gene_calls_dict)
            self.run.warning("Please read this carefully: Your external gene calls contained open reading frames with internal\
                              stop codons, and you asked anvi'o to ignore those. Anvi'o replaced internal stop codons with 'X'\
                              characters, and stored them in the contigs database that way. %d of your genes, which corresponded\
                              to %.2f%% of the total %d genes, had internal stop codons. We hope you are happy." % \
                                        (num_genes_with_internal_stops, percent_genes_with_internal_stops, len(gene_calls_dict)))

        if number_of_impartial_gene_calls:
            self.run.warning('%d of your %d gene calls were impartial, hence the translated protein sequences for those\
                              were not stored in the database.' % (number_of_impartial_gene_calls, len(gene_calls_dict)))

        return gene_calls_dict, protein_sequences_dict



def check_gene_calls_dict(self, gene_calls_dict):
    if not isinstance(gene_calls_dict, type({})):
        raise ConfigError("Gene calls dict must be a dict instance :/")

    try:
        [int(g) for g in list(gene_calls_dict.keys())]
    except ValueError:
        raise ConfigError("Keys of a gene calls dict must be integers!")

    if False in [x['direction'] in ['f', 'r'] for x in list(gene_calls_dict.values())]:
        raise ConfigError("The values in 'direction' column can't be anything but 'f' (for forward)\
                            or 'r' (for reverse). You have other stuff, and it is not cool.")

    if False in [x['stop'] > x['start'] for x in list(gene_calls_dict.values())]:
        raise ConfigError("For each gene call, the stop position must be bigger than the start position.\
                            Your gene calls dict does not conform to that. If you have reverse gene calls\
                            you must use the 'direction' column to declare that.")

    if False in [(x['stop'] - float(x['start'])) % 3.0 == 0 for x in list(gene_calls_dict.values())]:
        raise ConfigError("Something is wrong with your gene calls. For every gene call, the (stop - start)\
                            should be multiply of 3. It is not the case for all, which is a deal breaker.")

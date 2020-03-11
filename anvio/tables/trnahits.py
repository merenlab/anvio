# -*- coding: utf-8
# pylint: disable=line-too-long
"""
    The purpose of this module is to help the dealing with tRNA genes
    in contigs.
"""

import os
import shutil

from collections import Counter

import anvio
import anvio.utils as utils
import anvio.terminal as terminal
import anvio.constants as constants
import anvio.filesnpaths as filesnpaths
import anvio.drivers.trnscan_se as trnascan_se

from anvio.tables.hmmhits import TablesForHMMHits
from anvio.tables.genefunctions import TableForGeneFunctions


__author__ = "Developers of anvi'o (see AUTHORS.txt)"
__copyright__ = "Copyleft 2015-2020, the Meren Lab (http://merenlab.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "A. Murat Eren"
__email__ = "a.murat.eren@gmail.com"



class TablesForTransferRNAs:
    """ This class follows the structure of how Ribosomal RNAs are populated in contigs databases. A
        critical point is that ribosomal RNAs are identified in contigs context (rather than gene
        context), JUST LIKE the transfer RNAs will be, but with a major difference: we identify rRNAs
        through HMMs, and so this logic for ribosomal RNAs has ben neatly implemented in anvio/tables/hmmhits.
        But for tRNAs we don't have HMMs, therefore we can't directly benefit from the comprehensive HMM
        infrastructure of anvi'o for tRNAs. BUT WE MUST be able to use it to not implement a lot of
        redundant code to get tRNAs into contigs databases while still benefiting from what is already
        in place. HENCE This class. Which knits a workflow parallel to how ribosomal RNAs are added to
        contigs databases through the HMMs. But instaed of using models, it uses tRNAScan-SE output by
        default.
    """

    def __init__(self, args, run=terminal.Run(), progress=terminal.Progress()):
        self.args = args
        self.run = run
        self.progress = progress

        # just to make sure we have what it takes to continue later:
        trnascandriver = trnascan_se.tRNAScanSE(self.args, skip_sanity_check=True)
        trnascandriver.check_programs(quiet=True)

        self.tmp_directory_path = filesnpaths.get_temp_directory_path()

        P = lambda p: os.path.abspath(os.path.expanduser(p))
        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        self.num_threads = A('num_threads') or 1
        self.hits_file_path = P(A('trna_hits_file') or os.path.join(self.tmp_directory_path, 'hits_file.txt'))
        self.log_file_path = P(A('log_file') or os.path.join(self.tmp_directory_path, 'log.txt'))
        self.cutoff_score = A('trna_cutoff_score') or 20

        self.amino_acids = set([aa for aa in constants.AA_to_codons.keys() if aa != 'STP'])
        self.codons = set([cdn for cdn in constants.codon_to_AA.keys() if constants.codon_to_AA[cdn] in self.amino_acids])

        # the following variable is to meet the requirements of TablesForHMMHits to work with a new HMM
        # source.
        self.source = {'ref': 'Chan and Lowe, https://doi.org/10.1007/978-1-4939-9173-0_1',
                       'kind': 'Transfer_RNAs',
                       'domain': None,
                       'genes': ['%s_%s' % (constants.codon_to_AA[cdn], cdn) for cdn in self.codons],
                       'target': 'RNA:CONTIG',
                       'noise_cutoff_terms': None,
                       'model': None}

        self.source_name = 'Transfer_RNAs'
        self.kind_of_search = self.source['kind']
        self.domain = self.source['domain']
        self.all_genes_searched_against = self.source['genes']
        self.hmm_model = self.source['model']
        self.reference = self.source['ref']
        self.noise_cutoff_terms = self.source['noise_cutoff_terms']


    def run_trnascan_on_FASTA(self, fasta_file_path):
        self.run.info("Temporary dir", self.tmp_directory_path)
        self.run.info("FASTA file for contigs", fasta_file_path)
        self.run.info("tRNA hits output", self.hits_file_path)
        self.run.info("Log file", self.log_file_path)
        self.run.info("Cutoff score", self.cutoff_score)

        self.args.fasta_file = fasta_file_path
        self.args.log_file = self.log_file_path
        self.args.trna_hits_file = self.hits_file_path

        trnascandriver = trnascan_se.tRNAScanSE(self.args)
        results_dict = trnascandriver.process()

        return results_dict


    def populate_search_tables(self, contigs_db_path):
        utils.is_contigs_db(contigs_db_path)
        filesnpaths.is_output_file_writable(contigs_db_path, ok_if_exists=True)

        contig_sequences_fasta_path = os.path.join(self.tmp_directory_path, 'contig_sequences.fa')

        utils.export_sequences_from_contigs_db(contigs_db_path,
                                               contig_sequences_fasta_path)

        search_results_dict = self.run_trnascan_on_FASTA(fasta_file_path=contig_sequences_fasta_path)

        # At this point we need to turn this search_results_dict into one that matches how it is used
        # in HMM operations. Here is an entry from tRNA results dict:
        #
        # {1: {'contig': 'Bfragilis_0100_000000000001',
        #      'trna_no': '1',
        #      'start': 135361,
        #      'stop': 135433,
        #      'amino_acid': 'Thr',
        #      'codon': 'CGT',
        #      'score': 67.6}}
        #
        # and here is one exmple from the rRNA HMMs results dict:
        #
        # {1: {'entry_id': 0,
        #      'gene_name': 'Bacterial_23S_rRNA',
        #      'gene_hmm_id': '-',
        #      'contig_name': 'Bfragilis_0100_000000000001',
        #      'start': 1110877,
        #      'stop': 1113757,
        #      'e_value': 0.0}}
        #
        # so we will have to make the former look like the latter. I have the feeling that the
        # score / e_value will cause issues later :(

        missing_amino_acids = Counter()
        missing_codons = Counter()
        entries_to_remove = set([])
        for entry_id in search_results_dict:
            entry = search_results_dict[entry_id]

            aa, codon = entry['amino_acid'], utils.rev_comp(entry['anticodon'])

            if codon not in self.codons:
                missing_codons[codon] += 1
                entries_to_remove.add(entry_id)
                continue

            if aa not in self.amino_acids:
                missing_amino_acids[aa] += 1
                entries_to_remove.add(entry_id)
                continue

            aa_codon = '%s_%s' % (aa, codon)

            entry['gene_name'] = aa_codon
            entry['e_value'] = entry['score']
            entry['gene_hmm_id'] = '-'
            if entry['stop'] > entry['start']:
                # so we are forward
                entry['start'] = entry['start'] - 1 # setting the pythonic start.
            else:
                # so this one is reverse
                entry['stop'] = entry['stop'] - 1

        for entry_id in entries_to_remove:
            search_results_dict.pop(entry_id)

        self.run.info("Num tRNA genes recovered", len(search_results_dict))

        if len(missing_codons):
            self.run.warning("While anvi'o was trying to parse the output from tRNAScan-SE, it "
                             "became clear that some of the codons the tool identified was not "
                             "known to anvi'o, so we conservatively discareded those entries. "
                             "Here is the list of codons that were discareded and their frequency "
                             "among your contigs: '%s'." % (', '.join(['%s (%d)' % (codon, missing_codons[codon]) for codon in missing_codons])),
                             header="WEIRD CODONS ALERT")

        if len(missing_amino_acids):
            self.run.warning("While anvi'o was trying to parse the output from tRNAScan-SE, it "
                             "run into some amino acid names that were not known to anvi'o. "
                             "All those entries are now gone :/ But here is the list of amino "
                             "acids and their frequencies: '%s'." % (', '.join(['%s (%d)' % (amino_acid, missing_amino_acids[amino_acid]) for amino_acid in missing_amino_acids])),
                             header="WEIRD AMINO ACIDS ALERT")

        search_results_dict = utils.get_pruned_HMM_hits_dict(search_results_dict)

        tables_for_hmm_hits = TablesForHMMHits(contigs_db_path, run=self.run, progress=self.progress)
        search_results_dict = tables_for_hmm_hits.add_new_gene_calls_to_contigs_db_and_update_serach_results_dict(self.kind_of_search,
                                                                                                                  search_results_dict,
                                                                                                                  skip_amino_acid_sequences=True)
        tables_for_hmm_hits.append(self.source_name, self.reference, self.kind_of_search, self.domain, self.all_genes_searched_against, search_results_dict)


        # when the code comes all the way here, the entries in the search results dict already look like
        # this, so we have a gene callers id for the newly generate genes for tRNAs. we will use it
        # to populate a functions dict and submit it to the contigs database as well:
        #
        #     {'contig_name': 'Bfragilis_0100_000000000001',
        #      'trna_no': '1',
        #      'start': 135361,
        #      'stop': 135433,
        #      'amino_acid': 'Thr',
        #      'anticodon': 'CGT',
        #      'score': 67.6,
        #      'gene_name': 'Thr_ACG',
        #      'e_value': 67.6,
        #      'gene_hmm_id': '-',
        #      'gene_callers_id': 4502}
        #
        functions_dict = {}
        for entry_id in search_results_dict:
            entry = search_results_dict[entry_id]

            function_text = 'tRNA gene for amino acid %s (codon: %s; anticodon:%s; score:%.1f; intron_start:%d; intron_end:%d)' \
                                            % (entry['amino_acid'], codon, entry['anticodon'], \
                                               entry['score'], entry['intron_start'], entry['intron_end'])

            functions_dict[entry_id] = {'gene_callers_id': entry['gene_callers_id'],
                                        'source': self.source_name,
                                        'accession': '%s_%d' % (aa_codon, entry['gene_callers_id']),
                                        'function': function_text,
                                        'e_value': 0.0}

        gene_function_calls_table = TableForGeneFunctions(contigs_db_path, run=terminal.Run(verbose=False), progress=terminal.Progress(verbose=False))
        gene_function_calls_table.create(functions_dict)


        if not anvio.DEBUG:
            self.clean_tmp_directory()
            self.run.info_single("Temp directory is now cleaned (if you would like to keep it the "
                                 "next time use the flag `--debug`).", nl_before=1)
        else:
            self.run.info_single("Due to the `--debug` flag, anvi'o did not remove the temoporary files "
                                 "directory (which is still at '%s')." % (self.tmp_directory_path), nl_before=1)


    def clean_tmp_directory(self):
        shutil.rmtree(self.tmp_directory_path)

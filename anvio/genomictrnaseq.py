import os
import tempfile
import pandas as pd

import anvio
import anvio.utils as utils
import anvio.hmmops as hmmops
import anvio.tables as tables
import anvio.fastalib as fastalib
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths
import anvio.ccollections as ccollections

from itertools import combinations, product

from anvio.dbinfo import DBInfo
from anvio.errors import ConfigError
from anvio.drivers.blast import BLAST


__author__ = "Developers of anvi'o (see AUTHORS.txt)"
__copyright__ = "Copyleft 2015-2022, the Meren Lab (http://merenlab.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__maintainer__ = "Samuel Miller"
__email__ = "samuelmiller10@gmail.com"


progress = terminal.Progress()
run = terminal.Run()
pp = terminal.pretty_print

class SeedPermuter(object):
    default_min_nt_frequency = 0.05
    default_max_variable_positions = 5

    def __init__(self, args={}, p=progress, r=run, do_sanity_check=True):
        self.progress = p
        self.run = r

        self.args = args
        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        self.contigs_db_path = A('trnaseq_contigs_db')
        self.modifications_txt_path = A('modifications_txt')
        self.min_nt_frequency = A('min_nt_frequency')
        if self.min_nt_frequency == None:
            self.min_nt_frequency = self.default_min_nt_frequency
        self.max_variable_positions = A('max_variable_positions')
        if self.max_variable_positions == None:
            self.max_variable_positions = self.default_max_variable_positions
        self.permuted_seeds_fasta_path = A('permuted_seeds_fasta')

        self.tmp_dir = None
        if self.permuted_seeds_fasta_path == None:
            self.tmp_dir = anvio.TMP_DIR if anvio.TMP_DIR else tempfile.gettempdir()
            self.permuted_seeds_fasta_path = os.path.join(self.tmp_dir, "permuted_seeds.fa")

        if do_sanity_check:
            self.sanity_check()

        self.contigs_db_info = DBInfo(self.contigs_db_path, expecting='contigs')


    def sanity_check(self):
        contigs_db_info = DBInfo(self.contigs_db_path, expecting='contigs')
        if contigs_db_info.variant != 'trnaseq':
            raise ConfigError(f"The database at '{self.contigs_db_path}' was a '{contigs_db_info.variant}' variant, not the required 'trnaseq' variant.")

        filesnpaths.is_file_exists(self.modifications_txt_path)

        if not (0 <= self.min_nt_frequency < 1):
            raise ConfigError(f"The specified minimum nucleotide frequency for permutation is '{self.min_nt_frequency}', "
                              "which is not in the required range [0, 1).")

        if self.max_variable_positions < 1:
            raise ConfigError(f"The specified maximum number of variable positions in a permuted sequence is '{self.max_variable_positions}', "
                              "but it needs to be an integer with a value of 1 or greater.")

        filesnpaths.is_output_file_writable(self.permuted_seeds_fasta_path)


    def go(self):
        with self.contigs_db_info.load_db() as contigs_db:
            seed_contig_names_seqs = contigs_db.get_table_as_list_of_tuples('contig_sequences')

        variable_nts_gb = self.get_variable_nts_gb()

        output_fasta = fastalib.FastaOutput(self.permuted_seeds_fasta_path)
        self.run.info("FASTA file of permuted seeds", self.permuted_seeds_fasta_path)

        seed_count = 0
        seed_total = len(seed_contig_names_seqs)
        permuted_seq_count = 0
        max_permuted_seq_count = 0
        self.progress.new("Permuting seed sequences", progress_total_items=seed_total)
        for seed_contig_name, seed_seq in seed_contig_names_seqs:
            self.progress.update(f"{seed_count}/{seed_total}")
            seed_count += 1
            self.progress.increment()

            try:
                variable_nts_df = variable_nts_gb.get_group(seed_contig_name)
            except KeyError:
                # No modification positions were predicted in the seed.
                output_fasta.write_id(seed_contig_name)
                output_fasta.write_seq(seed_seq)
                continue

            variable_nts_dict = self.get_variable_nts_dict(variable_nts_df)
            if not variable_nts_dict:
                # No nts beside the most abundant had relative frequencies above the minimum threshold
                # at any predicted modification positions in the seed.
                output_fasta.write_id(seed_contig_name)
                output_fasta.write_seq(seed_seq)
                continue

            permuted_seed_info = self.get_permuted_seed_info(variable_nts_dict, seed_seq)
            self.remove_original_seed_seq(permuted_seed_info, seed_seq)

            permuted_seq_count += len(permuted_seed_info)
            max_permuted_seq_count = max(len(permuted_seed_info), max_permuted_seq_count)

            output_fasta.write_id(seed_contig_name)
            output_fasta.write_seq(seed_seq)
            for permuted_seq, permuted_positions, permuted_nts in permuted_seed_info:
                seq_id = f"{seed_contig_name}|{'_'.join([str(permuted_position) + permuted_nt for permuted_position, permuted_nt in zip(permuted_positions, permuted_nts)])}"
                output_fasta.write_id(seq_id)
                output_fasta.write_seq(permuted_seq)
        self.progress.end()

        self.run.info("Mean permuted seqs per seed", round(permuted_seq_count / seed_total, 1))
        self.run.info("Max permuted seqs from a seed", pp(max_permuted_seq_count))


    def get_variable_nts_gb(self):
        modifications_df = pd.read_csv(self.modifications_txt_path, sep='\t', header=0, usecols=['contig_name', 'seed_position', 'A', 'C', 'G', 'T'])
        modifications_df = modifications_df.set_index('contig_name')
        modifications_df = modifications_df.dropna()
        variable_nts_gb = modifications_df.groupby('contig_name')

        return variable_nts_gb


    def get_variable_nts_dict(self, variable_nts_df):
        variable_nts_dict = {}
        for position, position_variable_nts_df in variable_nts_df.set_index('seed_position').groupby('seed_position', sort=False):
            total_nts_df = position_variable_nts_df.sum(axis=0)
            total_nts_frequency_series = total_nts_df / total_nts_df.sum()
            total_nts_series = total_nts_df[total_nts_frequency_series > self.min_nt_frequency]
            total_nts_series = total_nts_series.loc[~total_nts_series.index.isin([total_nts_series.idxmax()])]

            if len(total_nts_series) == 0:
                continue

            variable_nts_dict[position] = tuple(total_nts_series.index)
        return variable_nts_dict


    def get_permuted_seed_info(self, variable_nts_dict, seed_seq):
        permuted_seed_info = []
        for num_variable_positions in range(1, min(len(variable_nts_dict), self.max_variable_positions) + 1):
            permutation_combinations = combinations(variable_nts_dict, num_variable_positions)
            for permuted_positions in permutation_combinations:
                for permuted_nts in product(*[variable_nts_dict[position] for position in permuted_positions]):
                    permuted_seed_seq = seed_seq
                    for position, nt in zip(permuted_positions, permuted_nts):
                        permuted_seed_seq = permuted_seed_seq[: position] + nt + permuted_seed_seq[position + 1: ]
                    permuted_seed_info.append((permuted_seed_seq, tuple(permuted_positions), tuple(permuted_nts)))
        return permuted_seed_info


    def remove_original_seed_seq(self, permuted_seed_info, seed_seq):
        original_index = None
        for i, seq_info in enumerate(permuted_seed_info):
            if seq_info[0] == seed_seq:
                original_index = i
                break
        if original_index == None:
            return False
        permuted_seed_info.pop(original_index)
        return True


class Integrator(object):
    default_max_mismatches = 3
    blast_search_output_cols = ['qseqid', 'sseqid', 'mismatch', 'qstart', 'qlen', 'sstart', 'send', 'slen', 'bitscore']

    def __init__(self, args={}, p=progress, r=run, do_sanity_check=True):
        self.progress = p
        self.run = r

        self.args = args
        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        self.trnaseq_contigs_db_path = A('trnaseq_contigs_db')
        self.seeds_specific_txt_path = A('seeds_specific_txt')
        self.modifications_txt_path = A('modifications_txt')
        self.genomic_contigs_db_path = A('contigs_db')
        self.genomic_profile_db_path = A('profile_db')
        self.collection_name = A('collection_name')
        self.num_threads = A('num_threads') or anvio.K('num_threads')['default']
        self.max_mismatches = A('max_mismatches')
        if self.max_mismatches == None:
            self.max_mismatches = self.default_max_mismatches
        self.use_full_length_seeds = A('use_full_length_seeds')
        if self.use_full_length_seeds == None:
            self.use_full_length_seeds = False
        self.blast_dir = A('blast_dir')
        self.permuted_seeds_fasta_path = A('permuted_seeds_fasta')
        self.just_do_it = A('just_do_it')

        if self.blast_dir == None:
            self.blast_dir = anvio.TMP_DIR if anvio.TMP_DIR else tempfile.gettempdir()
        self.trna_genes_fasta_path = os.path.join(self.blast_dir, 'trna_genes.fa')

        self.overwrite_table = self.just_do_it

        if do_sanity_check:
            self.sanity_check()

        self.bin_conscious = True if self.genomic_profile_db_path else False

        self.trnaseq_contigs_db_info = DBInfo(self.trnaseq_contigs_db_path)
        self.genomic_contigs_db_info = DBInfo(self.genomic_contigs_db_path)
        self.genomic_profile_db_info = DBInfo(self.genomic_profile_db_path) if self.genomic_profile_db_path else None


    def sanity_check(self, check_permuted_seeds_fasta=False):
        trnaseq_contigs_db_info = DBInfo(self.trnaseq_contigs_db_path, expecting='contigs')
        if trnaseq_contigs_db_info.variant != 'trnaseq':
            raise ConfigError(f"The database at '{self.trnaseq_contigs_db_path}' was a '{trnaseq_contigs_db_info.variant}' variant, not the required 'trnaseq' variant.")

        # The tRNA-seq contigs db version must be up-to-date to update the tRNA gene hits table.
        required_version = utils.get_required_version_for_db(self.trnaseq_contigs_db_path)
        if str(trnaseq_contigs_db_info.version) != required_version:
            raise ConfigError(f"The database at '{self.trnaseq_contigs_db_path}' is outdated (this database is v{trnaseq_contigs_db_info.version} "
                              f"and your anvi'o installation wants to work with v{required_version}). You can migrate your database without losing "
                              "any data using the program `anvi-migrate` with either of the flags `--migrate-dbs-safely` or `--migrate-dbs-quickly`.")

        with trnaseq_contigs_db_info.load_db() as database:
            stored_genomic_contigs_db_original_path = database.get_meta_value('genomic_contigs_db_original_path')
            stored_genomic_contigs_db_project_name = database.get_meta_value('genomic_contigs_db_project_name')
            stored_genomic_contigs_db_hash = database.get_meta_value('genomic_contigs_db_hash')
            stored_genomic_contigs_db_creation_hmd, stored_genomic_contigs_db_creation_mdy = utils.parse_epoch_time(database.get_meta_value('genomic_contigs_db_creation_date'))
            stored_genomic_profile_db_original_path = database.get_meta_value('genomic_profile_db_original_path')
            stored_genomic_profile_db_creation_hmd, stored_genomic_profile_db_creation_mdy = utils.parse_epoch_time(database.get_meta_value('genomic_profile_db_creation_date'))
            stored_genomic_profile_db_collection_name = database.get_meta_value('genomic_profile_db_collection_name')
        if stored_genomic_contigs_db_original_path and not self.overwrite_table:
            if stored_genomic_profile_db_original_path:
                additional_message = (f" Hits were restricted to the bin collection '{stored_genomic_profile_db_collection_name}' from the associated profile database "
                                      f"originally created at {stored_genomic_profile_db_creation_hmd} on {stored_genomic_profile_db_creation_mdy} (UTC).")
            else:
                additional_message = ""
            raise ConfigError(f"Use the flag `--just-do-it` to overwrite existing tRNA gene hits in the tRNA-seq database at '{self.genomic_contigs_db_path}'. "
                              f"The database currently has hits to a (meta)genomic contigs database with project name '{stored_genomic_contigs_db_project_name}' "
                              f"originally created at {stored_genomic_contigs_db_creation_hmd} on {stored_genomic_contigs_db_creation_mdy} (UTC) "
                              f"at the path '{stored_genomic_contigs_db_original_path}' "
                              f"and with the identifying hash '{stored_genomic_contigs_db_hash}'."
                              f"{additional_message}")

        filesnpaths.is_file_exists(self.seeds_specific_txt_path)
        filesnpaths.is_file_exists(self.modifications_txt_path)

        genomic_contigs_db_info = DBInfo(self.genomic_contigs_db_path, expecting='contigs')
        if genomic_contigs_db_info.variant != 'unknown':
            raise ConfigError(f"The database at '{self.genomic_contigs_db_path}' was a '{genomic_contigs_db_info.variant}' variant. "
                              "This should be a normal (meta)genomic contigs database, technically of an 'unknown' variant, produced by `anvi-gen-contigs-database`.")

        if 'Transfer_RNAs' not in hmmops.SequencesForHMMHits(self.genomic_contigs_db_path).hmm_hits_info:
            raise ConfigError(f"It appears that tRNA genes have not been annotated in the (meta)genomic contigs database, '{self.genomic_contigs_db_path}'. "
                              "Please run `anvi-scan-trnas` on the database and try again.")

        if (self.genomic_profile_db_path and not self.collection_name) or (self.collection_name and not self.genomic_profile_db_path):
            raise ConfigError("A profile database cannot be provided without a collection name, "
                              "nor can a collection name be provided without a profile database.")
        if self.genomic_profile_db_path:
            utils.is_profile_db_and_contigs_db_compatible(self.genomic_profile_db_path, self.genomic_contigs_db_path)

            collections = ccollections.Collections()
            collections.populate_collections_dict(self.genomic_profile_db_path)
            if self.collection_name not in collections.collections_dict:
                raise ConfigError(f"The desired collection, '{self.collection_name}', does not exist "
                                  f"in the (meta)genomic profile database, '{self.genomic_profile_db_path}'.")

        if self.num_threads < 1:
            raise ConfigError(f"The number of threads (used by BLAST) must be a positive integer, not the provided value of {self.num_threads}")

        if self.max_mismatches < 0:
            raise ConfigError("The maximum number of mismatches allowed in a seed-gene alignment "
                              f"must be a non-negative integer, not the provided value of {self.max_mismatches}")

        filesnpaths.is_output_dir_writable(self.blast_dir)

        # Ignore this sanity check when using `genomictrnaseq.Permuter` FASTA output.
        if check_permuted_seeds_fasta:
            filesnpaths.is_file_fasta_formatted(self.permuted_seeds_fasta_path)


    def go(self):
        self.write_trna_genes_fasta()
        search_output_path = self.blast()
        hits_df = self.filter_hits(search_output_path)
        unmodified_nt_df = self.find_unmodified_nucleotides(hits_df)
        self.update_trnaseq_contigs_database(hits_df, unmodified_nt_df)


    def write_trna_genes_fasta(self):
        # This method is based on `anvi-get-sequences-for-hmm-hits`.
        trna_gene_info = hmmops.SequencesForHMMHits(self.genomic_contigs_db_path, sources=set(['Transfer_RNAs']))
        splits_dict = {self.genomic_contigs_db_info.project_name: list(self.genomic_contigs_db_info.load_db().smart_get(tables.splits_info_table_name, 'split').keys())}
        hmm_seqs_dict = trna_gene_info.get_sequences_dict_for_hmm_hits_in_splits(splits_dict)

        if not len(hmm_seqs_dict):
            raise ConfigError("Unfortunately it appears that no tRNA genes were found in the contigs database.")

        trna_gene_info.store_hmm_sequences_into_FASTA(hmm_seqs_dict, self.trna_genes_fasta_path, header_section_separator='|', sequence_in_header=True)


    def blast(self):
        blast = BLAST(self.permuted_seeds_fasta_path, self.trna_genes_fasta_path, search_program='blastn', run=self.run, progress=self.progress, num_threads=self.num_threads)
        blast.tmp_dir = self.blast_dir
        blast.search_output_path = os.path.join(self.blast_dir, 'blast-search-results.txt')
        blast.log_file_path = os.path.join(self.blast_dir, 'blast-log.txt')
        blast.makedb(dbtype='nucl')
        blast.blast(outputfmt='6 ' + ' '.join(self.blast_search_output_cols), ungapped=True)

        return blast.search_output_path


    def filter_hits(self, search_output_path):
        hits_df = pd.read_csv(search_output_path, sep='\t', header=None, names=self.blast_search_output_cols)

        hits_df = hits_df[hits_df['mismatch'] <= self.max_mismatches]
        if self.use_full_length_seeds:
            hits_df = hits_df[hits_df['sstart'] == 1]

        hits_df[['seed_contig_name', 'seed_permutation']] = hits_df['qseqid'].str.split('|', expand=True)
        hits_df['seed_permutation'] = hits_df['seed_permutation'].fillna('')
        hits_df = hits_df.drop('qseqid', axis=1)

        contig_bin_id_dict = {}
        if self.bin_conscious:
            bin_contig_names_dict = ccollections.GetSplitNamesInBins(self.args).get_dict()
            for bin_id, split_names in bin_contig_names_dict.items():
                for split_name in split_names:
                    contig_bin_id_dict[split_name.split('_split_')[0]] = bin_id

        decoded_amino_acids = []
        anticodons = []
        bin_ids = []
        trnascan_scores = []
        gene_contig_names = []
        gene_callers_ids = []
        gene_starts = []
        gene_stops = []
        gene_sequences = []
        for sseqid in hits_df['sseqid']:
            split_sseqid = sseqid.split('|')

            decoded_amino_acid, anticodon = split_sseqid[0].split('_')[: 2]
            decoded_amino_acids.append(decoded_amino_acid)
            anticodons.append(anticodon)

            gene_contig_name = split_sseqid[4].split('contig:')[1]
            if self.bin_conscious:
                try:
                    bin_id = contig_bin_id_dict[gene_contig_name]
                except KeyError:
                    bin_id = ''
            else:
                bin_id = ''
            gene_contig_names.append(gene_contig_name)
            bin_ids.append(bin_id)

            trnascan_scores.append(split_sseqid[3].split('e_value:')[1])
            gene_callers_ids.append(split_sseqid[5].split('gene_callers_id:')[1])
            gene_starts.append(split_sseqid[6].split('start:')[1])
            gene_stops.append(split_sseqid[7].split('stop:')[1])
            gene_sequences.append(split_sseqid[9].split('sequence:')[1])
        hits_df['decoded_amino_acid'] = decoded_amino_acids
        hits_df['anticodon'] = anticodons
        hits_df['bin_id'] = bin_ids
        hits_df['trnascan_score'] = trnascan_scores
        hits_df['gene_contig_name'] = gene_contig_names
        hits_df['gene_callers_id'] = gene_callers_ids
        hits_df['gene_start_in_contig'] = gene_starts
        hits_df['gene_stop_in_contig'] = gene_stops
        hits_df['gene_sequence'] = gene_sequences

        # Filter individual alignments.
        retained_indices = []
        for index, decoded_amino_acid, seed_alignment_start, seed_length, gene_alignment_start, gene_alignment_end, gene_sequence, gene_length in zip(hits_df.index, hits_df['decoded_amino_acid'], hits_df['qstart'], hits_df['qlen'], hits_df['sstart'], hits_df['send'], hits_df['gene_sequence'], hits_df['slen']):
            if (gene_length - gene_alignment_end == 0) or ((gene_length - gene_alignment_end == 3) and gene_sequence[-3: ] == 'CCA'):
                # The alignment ends at the end of the gene or just short of a 3'-CCA acceptor in the gene (the seed should never contain the 3'-CCA acceptor).
                if (seed_alignment_start == 1) and (gene_alignment_end - gene_alignment_start == seed_length - 1):
                    # The alignment starts at the beginning of the gene and spans the entire query.
                    retained_indices.append(index)
                elif (decoded_amino_acid == 'His') and (seed_alignment_start == 2) and (gene_alignment_end - gene_alignment_start == seed_length - 2):
                    # The alignment starts at the second position of the tRNA-His seed sequence, which has a post-transcriptional G at the 5' end, and spans the remaining length of the query.
                    retained_indices.append(index)
        hits_df = hits_df.loc[retained_indices]

        # Retain each seed's top-scoring hits.
        hits_df = hits_df[hits_df.groupby('seed_contig_name')['bitscore'].transform('max') == hits_df['bitscore']]

        # If a collection of bins is provided, retain seeds unique to a single bin: disregard seeds
        # that equally strongly hit genes inside and outside a single bin, and disregard seeds that
        # hit unbinned contigs.
        if self.bin_conscious:
            hits_df = hits_df.groupby('seed_contig_name').filter(lambda seed_df: seed_df['bin_id'].nunique() == 1)
            hits_df = hits_df[hits_df['bin_id'] != '']

        # Multiple permutations of the same seed may be retained after filtering by score. There are
        # two and possibly more ways that this can occur. (1) The unmodified nucleotide at a
        # modified position has a very low frequency and so is not used in the permuted sequences.
        # The permuted sequences, none of which contain the correct nucleotide, match this
        # nucleotide in the gene equally well. (2) A permutation is introduced at a predicted
        # modification position that is actually a single nucleotide variant, and different versions
        # of the SNV occur in different (meta)genomic contigs. The following mechanism resolves both
        # of these possibilities, with the last step being the one that resolves the first
        # possibility. (1) Choose the permutation with the most hits. (2) If not resolved, choose
        # the permutation with the fewest permuted positions. (3) If not resolved, break the tie by
        # choosing the first permutation in the table, which will favor permutations toward the 5'
        # end.
        are_permutations_unresolved = True
        if hits_df.groupby('seed_contig_name').ngroups == hits_df.groupby(['seed_contig_name', 'seed_permutation']).ngroups:
            are_permutations_unresolved = False
        if are_permutations_unresolved:
            hits_df['count'] = hits_df.groupby(['seed_contig_name', 'seed_permutation'], as_index=False)['seed_contig_name'].transform(len)
            hits_df = hits_df[hits_df['count'] == hits_df.groupby('seed_contig_name')['count'].transform('max')]
            hits_df = hits_df.drop('count', axis=1)
            if hits_df.groupby('seed_contig_name').ngroups == hits_df.groupby(['seed_contig_name', 'seed_permutation']).ngroups:
                are_permutations_unresolved = False
        if are_permutations_unresolved:
            hits_df['num_permuted_positions'] = hits_df['seed_permutation'].apply(lambda p: p.count('_'))
            hits_df = hits_df[hits_df['num_permuted_positions'] == hits_df.groupby('seed_contig_name')['num_permuted_positions'].transform('min')]
            hits_df = hits_df.drop('num_permuted_positions', axis=1)
            if hits_df.groupby('seed_contig_name').ngroups == hits_df.groupby(['seed_contig_name', 'seed_permutation']).ngroups:
                are_permutations_unresolved = False
        if are_permutations_unresolved:
            hits_df = hits_df[hits_df['seed_permutation'] == hits_df.groupby('seed_contig_name')['seed_permutation'].transform('first')]

        # Seeds can be artifacts of the anvi'o de novo workflow, especially in relatively deeply
        # sequenced samples with high coverages, such as tRNA-seq libraries of pure cultures.
        # `anvi-merge-trnaseq` reports up to the number of seeds set by the user. If the user asks
        # for 1,000 seeds from a bacterial isolate experiment, then ~25-50 of these seeds will be
        # true tRNA sequences and up to ~950-975 will be artifacts (with unaccounted
        # modification-induced indels, nontemplated nucleotides, sequence errors, etc.), typically
        # at low frequency, that could not be resolved as non-tRNA by the tRNA-seq workflow. To
        # remove these artifact seeds, hits to the same gene are sorted by number of mismatches in
        # the alignment and seed abundance, and only the lowest mismatch/highest seed abundance hit
        # is retained. Seed abundance is taken as the average of relative abundance in each sample
        # based on 3' (discriminator nucleotide) coverage of the seed. For example, if there are two
        # tRNA-seq samples in the experiment, and two seeds hit the same gene each with one
        # mismatch, but one seed has relative 3' abundances of 0.02 and 0.03 in the two samples and
        # the other seed has abundances of 0.0006 and 0.00008, then the hit to the former seed will
        # be the only one retained for this gene.
        coverage_df = pd.read_csv(self.seeds_specific_txt_path, sep='\t', header=0, skiprows=[1, 2])
        coverage_df = coverage_df[['contig_name', 'sample_name', 'relative_discriminator_coverage']]
        coverage_df = coverage_df.rename({'contig_name': 'seed_contig_name'}, axis=1)
        seed_contig_names = hits_df['seed_contig_name'].unique()
        coverage_df = coverage_df[coverage_df['seed_contig_name'].isin(seed_contig_names)]
        coverage_df = hits_df.merge(coverage_df, on='seed_contig_name')
        def filter_multiple_hits_to_gene(gene_df):
            min_mismatch_df = gene_df[gene_df['mismatch'] == gene_df['mismatch'].min()]
            if min_mismatch_df['seed_contig_name'].nunique() > 1:
                most_abund_seed_contig_name = min_mismatch_df.groupby('seed_contig_name')['relative_discriminator_coverage'].mean().sort_values(ascending=False).index[0]
                return min_mismatch_df[min_mismatch_df['seed_contig_name'] == most_abund_seed_contig_name]
            return min_mismatch_df
        coverage_df = coverage_df.groupby('gene_callers_id', group_keys=False).apply(filter_multiple_hits_to_gene)
        hits_df = hits_df[hits_df['seed_contig_name'].isin(coverage_df['seed_contig_name'].unique())]

        # Spruce up the columns.
        hits_df['seed_alignment_start'] = hits_df['qstart'] - 1
        hits_df['gene_alignment_start'] = hits_df['sstart'] - 1
        hits_df = hits_df.drop(['qstart', 'qlen', 'sseqid', 'sstart'], axis=1)
        hits_df = hits_df.rename({'send': 'gene_alignment_stop'}, axis=1)
        hits_df = hits_df[['seed_contig_name', 'seed_permutation', # seed info
                           'gene_contig_name', 'gene_callers_id', 'gene_start_in_contig', 'gene_stop_in_contig', 'bin_id', 'trnascan_score', 'decoded_amino_acid', 'anticodon', 'gene_sequence', # gene info
                           'mismatch', 'bitscore', 'seed_alignment_start', 'gene_alignment_start', 'gene_alignment_stop']] # hit info

        return hits_df


    def find_unmodified_nucleotides(self, hits_df):
        modifications_df = pd.read_csv(self.modifications_txt_path, sep='\t', header=0, usecols=['contig_name', 'seed_position'])
        modifications_df = modifications_df.rename({'contig_name': 'seed_contig_name'}, axis=1)
        modifications_df = modifications_df.drop_duplicates()
        modifications_df = modifications_df[modifications_df['seed_contig_name'].isin(hits_df['seed_contig_name'].unique())]

        # Find the unmodified nucleotides at predicted modification positions in the seeds using the
        # matching tRNA gene sequences. If a seed matches multiple genes and the nucleotides at a
        # predicted modification position differ between the genes, then it is likely that the
        # variation is genetic rather than caused by a modification.
        modification_candidates_df = modifications_df.merge(hits_df[['seed_contig_name', 'seed_alignment_start', 'gene_sequence']], how='left', on='seed_contig_name')
        modification_keys = []
        unmodified_nts = []
        snv_keys = []
        for group_key, modification_candidate_df in modification_candidates_df.groupby(['seed_contig_name', 'seed_position'], as_index=False):
            unmodified_nt = ''
            for seed_position, seed_alignment_start, gene_sequence in zip(modification_candidate_df['seed_position'], modification_candidate_df['seed_alignment_start'], modification_candidate_df['gene_sequence']):
                gene_nt = gene_sequence[int(seed_position - seed_alignment_start)]
                if unmodified_nt:
                    if gene_nt != unmodified_nt:
                        snv_keys.append(group_key)
                        break
                else:
                    unmodified_nt = gene_nt
            else:
                modification_keys.append(group_key)
                unmodified_nts.append(unmodified_nt)
        modifications_df = modifications_df.set_index(['seed_contig_name', 'seed_position'])
        modifications_df = modifications_df.loc[modification_keys]
        unmodified_nt_df = modifications_df.merge(pd.DataFrame([modification_key + (unmodified_nt, ) for modification_key, unmodified_nt in zip(modification_keys, unmodified_nts)], columns=['seed_contig_name', 'seed_position', 'unmodified_nt']), on=['seed_contig_name', 'seed_position'])
        unmodified_nt_df = unmodified_nt_df.reset_index(drop=True)

        return unmodified_nt_df


    def update_trnaseq_contigs_database(self, hits_df, unmodified_nt_df):
        trnaseq_contigs_db = self.trnaseq_contigs_db_info.load_db()
        if self.overwrite_table:
            trnaseq_contigs_db._exec(f'''DELETE FROM {tables.trna_gene_hits_table_name}''')

        # Set meta-values in the tRNA-seq contigs database to track the associated (meta)genomic databases and collection.
        genomic_contigs_db = self.genomic_contigs_db_info.load_db()
        trnaseq_contigs_db.set_meta_value('genomic_contigs_db_original_path', os.path.abspath(self.genomic_contigs_db_path))
        trnaseq_contigs_db.set_meta_value('genomic_contigs_db_project_name', genomic_contigs_db.get_meta_value('project_name'))
        trnaseq_contigs_db.set_meta_value('genomic_contigs_db_hash', genomic_contigs_db.get_meta_value('contigs_db_hash'))
        trnaseq_contigs_db.set_meta_value('genomic_contigs_db_creation_date', genomic_contigs_db.get_meta_value('creation_date'))
        genomic_contigs_db.disconnect()
        if self.genomic_profile_db_path:
            genomic_profile_db = self.genomic_profile_db_info.load_db()
            trnaseq_contigs_db.set_meta_value('genomic_profile_db_original_path', os.path.abspath(self.genomic_profile_db_path))
            trnaseq_contigs_db.set_meta_value('genomic_profile_db_creation_date', genomic_profile_db.get_meta_value('creation_date'))
            trnaseq_contigs_db.set_meta_value('genomic_profile_db_collection_name', self.collection_name)
            genomic_profile_db.disconnect()

        table_entries = []
        for hit_id, row in hits_df.iterrows():
            seed_unmodified_nt_df = unmodified_nt_df[unmodified_nt_df['seed_contig_name'] == row['seed_contig_name']]

            if len(seed_unmodified_nt_df):
                seed_unmodified_nt_series = seed_unmodified_nt_df['seed_position'].astype(str) + seed_unmodified_nt_df['unmodified_nt']
                unmodified_nt_entry = ','.join(seed_unmodified_nt_series.tolist())
            else:
                unmodified_nt_entry = ''

            # Add somewhat more than the minimum amount of information (gene callers ID of the tRNA
            # gene) to make life easier and allow inspection of the result in the absence of the
            # (meta)genomic contigs database source.
            table_entries.append([hit_id,
                                  row['seed_contig_name'],
                                  row['gene_callers_id'],
                                  row['bin_id'],
                                  row['decoded_amino_acid'],
                                  row['anticodon'],
                                  row['mismatch'],
                                  row['bitscore'],
                                  row['seed_alignment_start'],
                                  row['gene_alignment_start'],
                                  row['gene_alignment_stop'],
                                  row['trnascan_score'],
                                  unmodified_nt_entry,
                                  row['gene_sequence']])
        trnaseq_contigs_db._exec_many(f'''INSERT INTO {tables.trna_gene_hits_table_name} VALUES ({",".join(["?"] * len(tables.trna_gene_hits_table_structure))})''', table_entries)
        trnaseq_contigs_db.disconnect()



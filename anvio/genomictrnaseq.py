import os
import tempfile
import numpy as np
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
    """Using the `go` method, generates a FASTA file of permuted seed sequences from a tRNA-seq contigs database."""

    default_min_nt_frequency = 0.05
    default_max_variable_positions = 5

    def __init__(self, args={}, p=progress, r=run, do_sanity_check=True):
        self.progress = p
        self.run = r

        self.args = args

        A = lambda x: args.__dict__[x] if x in args.__dict__ else None

        # mandatory args
        self.contigs_db_path = A('trnaseq_contigs_db')
        self.modifications_txt_path = A('modifications_txt')

        # optional args for algorithm parameterization
        self.min_nt_frequency = A('min_nt_frequency')
        if self.min_nt_frequency == None:
            self.min_nt_frequency = self.default_min_nt_frequency
        self.max_variable_positions = A('max_variable_positions')
        if self.max_variable_positions == None:
            self.max_variable_positions = self.default_max_variable_positions

        # paths for output files
        self.permuted_seeds_fasta_path = A('permuted_seeds_fasta')
        self.tmp_dir = None
        if self.permuted_seeds_fasta_path == None:
            self.tmp_dir = anvio.TMP_DIR if anvio.TMP_DIR else tempfile.gettempdir()
            self.permuted_seeds_fasta_path = os.path.join(self.tmp_dir, "permuted_seeds.fa")

        if do_sanity_check:
            self.sanity_check()

        self.contigs_db_info = DBInfo(self.contigs_db_path, expecting='contigs')


    def sanity_check(self):
        """Check the feasibility of args from initialization."""

        contigs_db_info = DBInfo(self.contigs_db_path, expecting='contigs')
        if contigs_db_info.variant != 'trnaseq':
            raise ConfigError(f"The database at '{self.contigs_db_path}' was a '{contigs_db_info.variant}' variant, "
                              "not the required 'trnaseq' variant.")

        filesnpaths.is_file_exists(self.modifications_txt_path)

        if not (0 <= self.min_nt_frequency < 1):
            raise ConfigError(f"The specified minimum nucleotide frequency for permutation is '{self.min_nt_frequency}', "
                              "which is not in the required range [0, 1).")

        if self.max_variable_positions < 1:
            raise ConfigError(f"The specified maximum number of variable positions in a permuted sequence is '{self.max_variable_positions}', "
                              "but it needs to be a positive integer.")

        filesnpaths.is_output_file_writable(self.permuted_seeds_fasta_path)


    def go(self):
        """Permute tRNA seed sequences at predicted sites of modification-induced substitutions."""

        with self.contigs_db_info.load_db() as contigs_db:
            seed_contig_names_seqs = contigs_db.get_table_as_list_of_tuples('contig_sequences')

        variable_nts_gb = self.get_variable_nts_gb()

        output_fasta = fastalib.FastaOutput(self.permuted_seeds_fasta_path)
        self.run.info("FASTA file of permuted seeds", self.permuted_seeds_fasta_path)

        # Permute tRNA-seq seeds with predicted modifications.
        # Write all seed sequences to the output FASTA along with permuted sequences.
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
                # No modified positions were predicted in the seed.
                output_fasta.write_id(seed_contig_name)
                output_fasta.write_seq(seed_seq)
                continue

            variable_nts_dict = self.get_variable_nts_dict(variable_nts_df, seed_seq)
            if not variable_nts_dict:
                # No nucleotides beside the most abundant had relative frequencies above the minimum
                # threshold for permutation at any predicted modification positions in the seed.
                output_fasta.write_id(seed_contig_name)
                output_fasta.write_seq(seed_seq)
                continue

            permuted_seed_info = self.get_permuted_seed_info(variable_nts_dict, seed_seq)

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
        """Load nucleotide variability data for modification-induced substitutions.

        Load the modifications table (generated by `anvi-tabulate-trnaseq`) and drop rows for
        samples lacking reported variability for the predicted modification.

        Returns
        =======
        variable_nts_gb : pandas DataFrameGroupBy
            the modifications table with a subset of columns grouped by seed ID ("contig name")
        """

        modifications_df = pd.read_csv(self.modifications_txt_path, sep='\t', header=0, usecols=['contig_name', 'seed_position', 'A', 'C', 'G', 'T'])
        modifications_df = modifications_df.set_index('contig_name')
        modifications_df = modifications_df.dropna()
        variable_nts_gb = modifications_df.groupby('contig_name')

        return variable_nts_gb


    def get_variable_nts_dict(self, variable_nts_df, seed_seq):
        """Find nucleotides that can be substituted at each predicted modified position in a seed.

        Each row in the input table contains data for a single seed + predicted modified position +
        sample. Samples without a predicted modified position due to a lack of variability (as found
        in `anvi-merge-trnaseq`) are not present in this table.

        For each seed + modified position, calculate the average relative frequency of each nucleotide
        across samples. Ignore the nucleotide in the seed sequence.

        Record the nucleotides meeting the relative frequency threshold.

        Parameters
        ==========
        variable_nts_df : pandas DataFrame
            a table of nucleotide frequencies at predicted modified positions in a seed

        seed_seq : str
            the seed sequence to permute

        Returns
        =======
        variable_nts_dict : dict
            a dict with keys being predicted modified positions in the seed and values being tuples
            of minority nucleotides (chars)
        """

        # Loop through each predicted modification in the seed.
        variable_nts_dict = {}
        for position, gb_df in variable_nts_df.set_index('seed_position').groupby('seed_position', sort=False):
            # Calculate nucleotide relative frequencies for each sample.
            position_variable_nts_df = gb_df.div(gb_df.sum(axis=1), axis=0)
            # Drop the column of the relative frequency of the nucleotide in the seed sequence.
            position_variable_nts_df = position_variable_nts_df.drop(seed_seq[position], axis=1)
            # Calculate the mean relative frequency of each nucleotide across samples.
            position_variable_nts_series = position_variable_nts_df.mean(axis=0)
            # Filter nucleotides by threshold relative frequency.
            position_variable_nts_series = position_variable_nts_series[position_variable_nts_series > self.min_nt_frequency]

            if len(position_variable_nts_series) == 0:
                continue

            variable_nts_dict[position] = tuple(position_variable_nts_series.index)
        return variable_nts_dict


    def get_permuted_seed_info(self, variable_nts_dict, seed_seq):
        """Generate permuted sequences, returning not just the new sequences but also the permuted
        positions and substituted nucleotides in the sequences.

        Parameters
        ==========
        variable_nts_dict : dict
            the dict recording nucleotides to substitute at positions in the seed

        seed_seq : str
            the seed sequence to permute

        Returns
        =======
        permuted_seed_info : list
            a list of tuples, each tuple containing a 1) permuted seed sequence string, 2) a tuple
            of permuted position indices, and 3) a tuple of substituted nucleotide characters at
            those positions
        """

        permuted_seed_info = []
        # Loop through the different numbers of permuted positions that can be introduced in the
        # sequence, starting with 1 permuted position.
        for num_variable_positions in range(1, min(len(variable_nts_dict), self.max_variable_positions) + 1):
            # Find combinations of permuted positions, e.g., with 1 permuted position, and positions
            # 8 and 32 being variable, then one combination would simply be (8, ) and the other
            # combination (32, ).
            permutation_combinations = combinations(variable_nts_dict, num_variable_positions)
            # Loop through each combination.
            for permuted_positions in permutation_combinations:
                # Find the sets of nucleotides that will be substituted into the sequence given the
                # combination of positions. Each loop generates a new permuted sequence and entry.
                for permuted_nts in product(*[variable_nts_dict[position] for position in permuted_positions]):
                    permuted_seed_seq = seed_seq
                    # Loop through the positions to make the nucleotide substitutions.
                    for position, nt in zip(permuted_positions, permuted_nts):
                        permuted_seed_seq = permuted_seed_seq[: position] + nt + permuted_seed_seq[position + 1: ]
                    # In addition to the permuted sequence, record the permuted positions and
                    # substituted nucleotides.
                    permuted_seed_info.append((permuted_seed_seq, tuple(permuted_positions), tuple(permuted_nts)))

        return permuted_seed_info


class Integrator(object):
    """Using the `go` method, links tRNA-seq seeds to tRNA genes and adds this information to the tRNA-seq contigs database."""

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
        if self.blast_dir == None:
            self.blast_dir = anvio.TMP_DIR if anvio.TMP_DIR else tempfile.gettempdir()
        self.permuted_seeds_fasta_path = A('permuted_seeds_fasta')
        self.just_do_it = A('just_do_it')

        self.trna_genes_fasta_path = os.path.join(self.blast_dir, 'trna_genes.fa')
        self.overwrite_table = self.just_do_it

        if do_sanity_check:
            self.sanity_check()

        self.bin_conscious = True if self.genomic_profile_db_path else False

        self.trnaseq_contigs_db_info = DBInfo(self.trnaseq_contigs_db_path)
        self.genomic_contigs_db_info = DBInfo(self.genomic_contigs_db_path)
        self.genomic_profile_db_info = DBInfo(self.genomic_profile_db_path) if self.genomic_profile_db_path else None


    def sanity_check(self, check_permuted_seeds_fasta=False):
        """Check the feasibility of args from initialization."""

        trnaseq_contigs_db_info = DBInfo(self.trnaseq_contigs_db_path, expecting='contigs')
        if trnaseq_contigs_db_info.variant != 'trnaseq':
            raise ConfigError(f"The database at '{self.trnaseq_contigs_db_path}' was a '{trnaseq_contigs_db_info.variant}' variant, "
                              "not the required 'trnaseq' variant.")

        # If the tRNA-seq seeds have already been associated with tRNA genes, then `--just-do-it` is needed.
        trnaseq_contigs_db_self_table = trnaseq_contigs_db_info.get_self_table()
        associated_genomic_contigs_db_hash = trnaseq_contigs_db_self_table['genomic_contigs_db_hash']
        associated_genomic_contigs_db_project_name = trnaseq_contigs_db_self_table['genomic_contigs_db_project_name']
        associated_genomic_collection_name = trnaseq_contigs_db_self_table['genomic_collection_name']
        if not self.overwrite_table and associated_genomic_contigs_db_hash != None:
            if associated_genomic_collection_name:
                additional_message = f" A collection named {associated_genomic_collection_name} was used."
            else:
                additional_message = ""
            raise ConfigError(f"The tRNA-seq contigs database at '{self.trnaseq_contigs_db_path}' "
                              "has already been associated with tRNA genes from a (meta)genomic contigs database "
                              f"with the project name, {associated_genomic_contigs_db_project_name}, "
                              f"and hash ID, {associated_genomic_contigs_db_hash}.{additional_message} "
                              f"`anvi-integrate-trnaseq` can be run with the flag `--just-do-it` to overwrite the existing data.")

        # The tRNA-seq contigs db version must be up-to-date to update the tRNA gene hits table.
        required_version = utils.get_required_version_for_db(self.trnaseq_contigs_db_path)
        if str(trnaseq_contigs_db_info.version) != required_version:
            raise ConfigError(f"The database at '{self.trnaseq_contigs_db_path}' is outdated (this database is v{trnaseq_contigs_db_info.version} "
                              f"and your anvi'o installation wants to work with v{required_version}). "
                              "You can migrate your database without losing any data using the program `anvi-migrate` "
                              "with either of the flags `--migrate-dbs-safely` or `--migrate-dbs-quickly`.")

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
        """Link tRNA-seq seeds to tRNA genes, adding this information to the tRNA-seq contigs database."""

        self.write_trna_genes_fasta()
        search_output_path = self.blast()
        hits_df = self.filter_hits(search_output_path)
        unmodified_nt_df = self.find_unmodified_nucleotides(hits_df)
        self.update_trnaseq_contigs_database(hits_df, unmodified_nt_df)


    def write_trna_genes_fasta(self):
        """Write tRNA genes in the (meta)genomic contigs database to a FASTA file.

        This method is based on `anvi-get-sequences-for-hmm-hits`.
        """

        trna_gene_info = hmmops.SequencesForHMMHits(self.genomic_contigs_db_path, sources=set(['Transfer_RNAs']))
        splits_dict = {self.genomic_contigs_db_info.project_name: list(self.genomic_contigs_db_info.load_db().smart_get(tables.splits_info_table_name, 'split').keys())}
        hmm_seqs_dict = trna_gene_info.get_sequences_dict_for_hmm_hits_in_splits(splits_dict)

        if not len(hmm_seqs_dict):
            raise ConfigError("Unfortunately it appears that no tRNA genes were found in the (meta)genomic contigs database.")

        trna_gene_info.store_hmm_sequences_into_FASTA(hmm_seqs_dict, self.trna_genes_fasta_path, header_section_separator='|', sequence_in_header=True)


    def blast(self):
        """Align permuted tRNA-seq seeds to tRNA genes."""

        blast = BLAST(self.permuted_seeds_fasta_path, self.trna_genes_fasta_path, search_program='blastn', run=self.run, progress=self.progress, num_threads=self.num_threads)
        blast.tmp_dir = self.blast_dir
        blast.search_output_path = os.path.join(self.blast_dir, 'blast-search-results.txt')
        blast.log_file_path = os.path.join(self.blast_dir, 'blast-log.txt')
        blast.makedb(dbtype='nucl')
        blast.blast(outputfmt='6 ' + ' '.join(self.blast_search_output_cols), ungapped=True)

        return blast.search_output_path


    def filter_hits(self, search_output_path):
        """The heavy lifting to confidently associate tRNA-seq seeds with tRNA genes."""

        # Load BLAST output table.
        hits_df = pd.read_csv(search_output_path, sep='\t', header=None, names=self.blast_search_output_cols)

        # Apply some alignment filters to hits.
        hits_df = hits_df[hits_df['mismatch'] <= self.max_mismatches]
        if self.use_full_length_seeds:
            hits_df = hits_df[hits_df['sstart'] == 1]

        # Separate seed IDs and permutation info.
        hits_df[['seed_contig_name', 'seed_permutation']] = hits_df['qseqid'].str.split('|', expand=True)
        hits_df['seed_permutation'] = hits_df['seed_permutation'].fillna('')
        hits_df = hits_df.drop('qseqid', axis=1)

        # If considering a collection of bins, associate seed "contigs" with bins.
        contig_bin_id_dict = {}
        if self.bin_conscious:
            bin_contig_names_dict = ccollections.GetSplitNamesInBins(self.args).get_dict()
            for bin_id, split_names in bin_contig_names_dict.items():
                for split_name in split_names:
                    contig_bin_id_dict[split_name.split('_split_')[0]] = bin_id

        # Extract information on each hit.
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
        for index, decoded_amino_acid, seed_alignment_start, seed_length, gene_alignment_start, gene_alignment_end, gene_sequence, gene_length in zip(
            hits_df.index, hits_df['decoded_amino_acid'], hits_df['qstart'], hits_df['qlen'], hits_df['sstart'], hits_df['send'], hits_df['gene_sequence'], hits_df['slen']):
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
        # modified position has a very low frequency and so was not used in the permuted sequences.
        # The permuted sequences, none of which contain the correct nucleotide, mismatch this
        # nucleotide in the gene equally well. (2) A permutation is introduced at a predicted
        # modification position that is actually a single nucleotide variant, and different versions
        # of the SNV occur in different (meta)genomic contigs. The following procedure resolves both
        # of these possibilities, with the last step being the one that resolves the first
        # possibility. (1) Choose the permutation hitting the greatest number of genes. (2) If not
        # resolved, choose the permutation with the fewest permuted positions. (3) If not resolved,
        # break the tie by choosing the first permutation in the table, which will favor
        # permutations toward the 5' end.
        are_permutations_unresolved = True
        if hits_df.groupby('seed_contig_name').ngroups == hits_df.groupby(['seed_contig_name', 'seed_permutation']).ngroups:
            are_permutations_unresolved = False
        if are_permutations_unresolved: # (1)
            hits_df['count'] = hits_df.groupby(['seed_contig_name', 'seed_permutation'], as_index=False)['seed_contig_name'].transform(len)
            hits_df = hits_df[hits_df['count'] == hits_df.groupby('seed_contig_name')['count'].transform('max')]
            hits_df = hits_df.drop('count', axis=1)
            if hits_df.groupby('seed_contig_name').ngroups == hits_df.groupby(['seed_contig_name', 'seed_permutation']).ngroups:
                are_permutations_unresolved = False
        if are_permutations_unresolved: # (2)
            hits_df['num_permuted_positions'] = hits_df['seed_permutation'].apply(lambda p: p.count('_'))
            hits_df = hits_df[hits_df['num_permuted_positions'] == hits_df.groupby('seed_contig_name')['num_permuted_positions'].transform('min')]
            hits_df = hits_df.drop('num_permuted_positions', axis=1)
            if hits_df.groupby('seed_contig_name').ngroups == hits_df.groupby(['seed_contig_name', 'seed_permutation']).ngroups:
                are_permutations_unresolved = False
        if are_permutations_unresolved: # (3)
            hits_df = hits_df[hits_df['seed_permutation'] == hits_df.groupby('seed_contig_name')['seed_permutation'].transform('first')]

        # Seeds can be artifacts of the anvi'o de novo tRNA-seq workflow, especially in relatively
        # deeply sequenced samples with high coverages, such as tRNA-seq libraries of pure cultures.
        # `anvi-merge-trnaseq` reports up to the number of seeds set by the user. If the user asks
        # for 1,000 seeds from a bacterial isolate experiment, then ~25-50 of these seeds will be
        # true tRNA sequences and up to ~950-975 will be artifacts (containing unaccounted
        # modification-induced indels, nontemplated nucleotides, sequence errors, etc., typically at
        # low frequency), that could not be resolved as non-tRNA by the tRNA-seq workflow. To remove
        # these artifact seeds, hits to the same gene are sorted by number of mismatches in the
        # alignment and seed abundance, and only the lowest mismatch/highest seed abundance hit is
        # retained. Seed abundance is taken as the average of relative abundance in each sample
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
        """Find the unmodified nucleotides at predicted modification positions in tRNA-seq seeds using matching tRNA gene sequences."""

        # Load modification information for seeds associated with genes.
        modifications_df = pd.read_csv(self.modifications_txt_path, sep='\t', header=0, usecols=['contig_name', 'seed_position'])
        modifications_df = modifications_df.rename({'contig_name': 'seed_contig_name'}, axis=1)
        modifications_df = modifications_df.drop_duplicates()
        modifications_df = modifications_df[modifications_df['seed_contig_name'].isin(hits_df['seed_contig_name'].unique())]

        # If a seed matches multiple genes and the nucleotides at a predicted modification position
        # differ between the genes, then it is likely that the variation is genetic rather than
        # caused by a modification.
        modification_candidates_df = modifications_df.merge(hits_df[['seed_contig_name', 'seed_alignment_start', 'gene_sequence']],
                                                            how='left', on='seed_contig_name')
        modification_keys = []
        unmodified_nts = []
        snv_keys = []
        for group_key, modification_candidate_df in modification_candidates_df.groupby(['seed_contig_name', 'seed_position'], as_index=False):
            unmodified_nt = ''
            for seed_position, seed_alignment_start, gene_sequence in zip(
                modification_candidate_df['seed_position'], modification_candidate_df['seed_alignment_start'], modification_candidate_df['gene_sequence']):
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
        unmodified_nt_df = modifications_df.merge(pd.DataFrame([modification_key + (unmodified_nt, ) for modification_key, unmodified_nt in zip(modification_keys, unmodified_nts)],
                                                               columns=['seed_contig_name', 'seed_position', 'unmodified_nt']),
                                                  on=['seed_contig_name', 'seed_position'])
        unmodified_nt_df = unmodified_nt_df.reset_index(drop=True)

        return unmodified_nt_df


    def update_trnaseq_contigs_database(self, hits_df, unmodified_nt_df):
        """Add information on tRNA gene associations to the tRNA-seq contigs database."""

        trnaseq_contigs_db = self.trnaseq_contigs_db_info.load_db()

        # Erase the table contents: if this table was already filled out in a previous run, then
        # `self.sanity_check` would catch this and force the user to `--just-do-it`.
        trnaseq_contigs_db._exec(f'''DELETE FROM {tables.trna_gene_hits_table_name}''')

        # Update metadata tracking (meta)genomic associations.
        trnaseq_contigs_db.set_meta_value('genomic_contigs_db_project_name', self.genomic_contigs_db_info.hash)
        trnaseq_contigs_db.set_meta_value('genomic_contigs_db_hash', self.genomic_contigs_db_info.project_name)
        if self.collection_name:
            trnaseq_contigs_db.set_meta_value('genomic_collections_name', self.collection_name)

        # Assemble the rows of the table.
        table_entries = []
        hit_id = 0
        for row in hits_df.itertuples(index=False):
            seed_unmodified_nt_df = unmodified_nt_df[unmodified_nt_df['seed_contig_name'] == row.seed_contig_name]

            if len(seed_unmodified_nt_df):
                seed_unmodified_nt_series = seed_unmodified_nt_df['seed_position'].astype(str) + seed_unmodified_nt_df['unmodified_nt']
                unmodified_nt_entry = ','.join(seed_unmodified_nt_series.tolist())
            else:
                unmodified_nt_entry = ''

            # Add somewhat more than the minimum amount of information (gene callers ID of the tRNA
            # gene) to make life easier and allow inspection of the result in the absence of the
            # (meta)genomic contigs database source.
            table_entries.append([hit_id,
                                  row.seed_contig_name,
                                  row.gene_callers_id,
                                  row.bin_id,
                                  row.decoded_amino_acid,
                                  row.anticodon,
                                  row.mismatch,
                                  row.bitscore,
                                  row.seed_alignment_start,
                                  row.gene_alignment_start,
                                  row.gene_alignment_stop,
                                  row.trnascan_score,
                                  unmodified_nt_entry,
                                  row.gene_sequence])
            hit_id += 1
        trnaseq_contigs_db._exec_many(f'''INSERT INTO {tables.trna_gene_hits_table_name} VALUES ({",".join(["?"] * len(tables.trna_gene_hits_table_structure))})''', table_entries)
        trnaseq_contigs_db.disconnect()

        self.run.info_single(f"In total, {hit_id} tRNA-seq seeds were linked to tRNA genes.")


class Affinitizer:
    """Using the `go` method, relates changes in tRNA-seq seed abundances to the codon usage of gene functions."""

    default_min_coverage = 10
    default_min_isoacceptors = 4
    default_rarefaction_limit = 0

    def __init__(self, args={}, p=progress, r=run, do_sanity_check=True):
        self.progress = p
        self.run = r

        self.args = args
        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        self.trnaseq_contigs_db_path = A('trnaseq_contigs_db')
        self.seeds_specific_txt_path = A('seeds_specific_txt')
        self.genomic_contigs_db_path = A('contigs_db')
        self.reference_sample_name = A('reference_sample')
        self.nonreference_sample_names = A('sample_subset')
        self.min_coverage = A('min_coverage')
        if self.min_coverage == None:
            self.min_coverage = self.default_min_coverage
        self.min_isoacceptors = A('min_isoacceptors')
        if self.min_isoacceptors == None:
            self.min_isoacceptors = self.default_min_isoacceptors
        self.rarefaction_limit = A('rarefaction_limit')
        if self.rarefaction_limit == None:
            self.rarefaction_limit = self.default_rarefaction_limit

        if do_sanity_check:
            self.sanity_check()

        self.trnaseq_contigs_db_info = DBInfo(self.trnaseq_contigs_db_path)
        self.genomic_contigs_db_info = DBInfo(self.genomic_contigs_db_path)

        # Find the names of the samples to analyze in addition to the reference sample: the
        # "nonreference" samples. By default, this will be every other sample, but a subset of
        # available samples can also be used.
        if self.nonreference_sample_names == None:
            self.nonreference_sample_names = pd.read_csv(self.seeds_specific_txt_path, sep='\t', header=0, skiprows=[1, 2], usecols=['sample_name'])['sample_name'].unique().tolist()
            self.nonreference_sample_names.remove(self.reference_sample_name)
        else:
            self.nonreference_sample_names = self.nonreference_sample_names.split(',')
        self.sample_names = [self.reference_sample_name] + self.nonreference_sample_names


    def sanity_check(self):
        """Check the feasibility of args from initialization."""

        trnaseq_contigs_db_info = DBInfo(self.trnaseq_contigs_db_path, expecting='contigs')
        if trnaseq_contigs_db_info.variant != 'trnaseq':
            raise ConfigError(f"The database at '{self.trnaseq_contigs_db_path}' was a '{trnaseq_contigs_db_info.variant}' variant, "
                              "not the required 'trnaseq' variant.")
        with trnaseq_contigs_db_info.load_db() as trnaseq_contigs_db:
            if len(trnaseq_contigs_db.get_table_as_dataframe('trna_gene_hits', columns_of_interest=['seed_contig_name'])) == 0:
                raise ConfigError(f"It appears that no tRNA seeds in the tRNA-seq contigs database, '{self.trnaseq_contigs_db_path}', "
                                  f"are linked to tRNA genes in the (meta)genomic contigs database, '{self.genomic_contigs_db_path}'.")

        filesnpaths.is_file_exists(self.seeds_specific_txt_path)

        genomic_contigs_db_info = DBInfo(self.genomic_contigs_db_path, expecting='contigs')
        if genomic_contigs_db_info.variant != 'unknown':
            raise ConfigError(f"The database at '{self.genomic_contigs_db_path}' was a '{genomic_contigs_db_info.variant}' variant. "
                              "This should be a normal (meta)genomic contigs database, technically an 'unknown' variant, produced by `anvi-gen-contigs-database`.")
        if 'modules_db_hash' not in genomic_contigs_db_info.get_self_table():
            raise ConfigError(f"It appears that genes have not been annotated by KOfams in the (meta)genomic contigs database, '{self.genomic_contigs_db_path}'. "
                              "Please run `anvi-run-kegg-kofams` on the database and try again.")

        available_sample_names = pd.read_csv(self.seeds_specific_txt_path, sep='\t', header=0, skiprows=[1, 2], usecols=['sample_name'])['sample_name'].unique().tolist()
        # Check for the existence of the reference sample.
        if self.reference_sample_name not in available_sample_names:
            raise ConfigError(f"The desired reference sample name, '{self.reference_sample_name}', "
                              f"was not found in `seeds-specific-txt`, '{self.seeds_specific_txt_path}'. "
                              f"Here are the samples provided in that table: {', '.join(available_sample_names)}")
        # Check for the existence of a given subset of sample names.
        if self.nonreference_sample_names:
            nonreference_sample_names = self.nonreference_sample_names.split(',')
            bad_sample_names = set(nonreference_sample_names).difference(set(available_sample_names))
            if bad_sample_names:
                raise ConfigError("The following desired subset sample names were not found in `seeds-specific-txt`, "
                                  f"'{self.seeds_specific_txt_path}': {', '.join(bad_sample_names)}. "
                                  f"Here are the samples provided in that table: {', '.join(available_sample_names)}")
            if self.reference_sample_name in nonreference_sample_names:
                raise ConfigError(f"Please do not include the reference sample, '{self.reference_sample_name}' "
                                  f"in the subset of sample names: {', '.join(nonreference_sample_names)}. Sorry for the sclerotic idiocy.")
        else:
            nonreference_sample_names = available_sample_names
            nonreference_sample_names.remove(self.reference_sample_name)
        if len(nonreference_sample_names) == 0:
            raise ConfigError(f"There must be one or more samples beside the reference sample in `seeds-specific-txt`, '{self.seeds_specific_txt}'. "
                              f"Only the reference sample, '{self.reference_sample_name}', was found.")

        if self.min_coverage < 1:
            raise ConfigError("The minimum coverage for a tRNA isoacceptor to be detected must be an integer "
                              f"greater than or equal to 1, not the provided value of {self.min_coverage}.")

        if self.min_isoacceptors < 1:
            raise ConfigError("The minimum number of tRNA isoacceptors for translational affinity to be calculated "
                              f"must be an integer greater or equal to 1, not the provided value of {self.min_isoacceptors}.")

        if self.rarefaction_limit < 0:
            raise ConfigError("The rarefaction limit on subsampled tRNA isoacceptors must be an integer "
                              f"greater than or equal to 1, not the provided value of {self.rarefaction_limit}.")


    def go(self):
        """Relate changes in tRNA-seq seed abundances to the codon usage of gene functions."""

        isoacceptors_df = self.load_isoacceptor_data()
        if len(isoacceptors_df) == 0:
            return
        isoacceptor_abundance_dict = self.get_isoacceptor_abundances(isoacceptors_df)
        kegg_df = self.consolidate_kegg_annotations()


    def load_isoacceptor_data(self):
        """Workhorse method to load and filter data, find properties of seeds, and ultimately group isoacceptors."""

        # Load data from the tRNA-seq contigs database.
        with self.trnaseq_contigs_db_info.load_db() as trnaseq_contigs_db:
            trna_gene_hits_df = trnaseq_contigs_db.get_table_as_dataframe(
                'trna_gene_hits', columns_of_interest=['seed_contig_name', 'bin_name', 'decoded_amino_acid', 'anticodon', 'gene_sequence'])
            # Drop duplicate rows (only preserve one row) representing hits between the same
            # tRNA-seq seed and different tRNA genes with identical sequences.
            trna_gene_hits_df = trna_gene_hits_df.drop_duplicates()
            if self.genomic_collection_name == None:
                trna_gene_hits_df['bin_name'] = ''

            self.run.info("tRNA-seq seeds linked to tRNA genes", len(trna_gene_hits_df))

            # Get the gene callers IDs of the tRNA-seq seeds.
            seed_contig_names = ','.join(['"%s"' % seed_contig_name for seed_contig_name in trna_gene_hits_df['seed_contig_name'].unique()])
            contigs_where_clause = f'''contig IN ({seed_contig_names})'''
            seed_id_df = trnaseq_contigs_db.get_table_as_dataframe(
                'genes_in_contigs', columns_of_interest=['gene_callers_id', 'contig'], where_clause=contigs_where_clause)
            seed_id_df = seed_id_df.rename({'contig': 'seed_contig_name', 'gene_callers_id': 'seed_gene_callers_id'}, axis=1)
            trna_gene_hits_df = trna_gene_hits_df.merge(seed_id_df, on='seed_contig_name')

            # Find wobble position 34 in the tRNA-seq seed sequences.
            seed_gene_callers_ids = ','.join(['"%s"' % gene_callers_id for gene_callers_id in trna_gene_hits_df['seed_gene_callers_id'].unique()])
            ids_where_clause = f'''gene_callers_id IN ({seed_gene_callers_ids})'''
            wobble_position_df = trnaseq_contigs_db.get_table_as_dataframe(
                'trna_feature', columns_of_interest=['gene_callers_id', 'anticodon_loop_start'], where_clause=ids_where_clause)
            wobble_position_df = wobble_position_df.rename({'gene_callers_id': 'seed_gene_callers_id'}, axis=1)
            wobble_position_df['anticodon_start'] = wobble_position_df['anticodon_loop_start'] + 2
            wobble_position_df = wobble_position_df.drop('anticodon_loop_start', axis=1)
            trna_gene_hits_df = trna_gene_hits_df.merge(wobble_position_df, on='seed_gene_callers_id')

            # Get the tRNA-seq seed consensus sequence strings.
            seed_consensus_sequence_df = trnaseq_contigs_db.get_table_as_dataframe('contig_sequences', where_clause=contigs_where_clause)
            seed_consensus_sequence_df = seed_consensus_sequence_df.rename({'contig': 'seed_contig_name', 'sequence': 'seed_sequence'}, axis=1)
            trna_gene_hits_df = trna_gene_hits_df.merge(seed_consensus_sequence_df, on='seed_contig_name')

            # Find the nucleotides read at wobble position 34 in the tRNA-seq seeds.
            anticodon_wobble_nucleotides = []
            for anticodon_start, seed_consensus_sequence in zip(trna_gene_hits_df['anticodon_start'], trna_gene_hits_df['seed_sequence']):
                anticodon_wobble_nucleotides.append(seed_consensus_sequence[anticodon_start])
            trna_gene_hits_df['seed_anticodon_wobble_nucleotide'] = anticodon_wobble_nucleotides

            trna_gene_hits_df = trna_gene_hits_df.drop(['seed_gene_callers_id', 'anticodon_start', 'seed_sequence'], axis=1)

        # Load data from the seeds specific coverage table.
        coverage_df = pd.read_csv(self.seeds_specific_txt_path, sep='\t', header=0, skiprows=[1, 2],
                                  usecols=['contig_name', 'sample_name', 'relative_discriminator_coverage', 'discriminator_1'])
        coverage_df = coverage_df.rename({'contig_name': 'seed_contig_name'}, axis=1)
        # Select data for the samples of interest.
        coverage_df = coverage_df[coverage_df['sample_name'].isin(self.sample_names)]
        # Select data for the tRNA-seq seeds linked to tRNA genes.
        coverage_df = coverage_df[coverage_df['seed_contig_name'].isin(trna_gene_hits_df['seed_contig_name'].unique())]
        coverage_df = coverage_df[coverage_df['discriminator_1'] >= self.min_coverage]
        coverage_df = coverage_df.drop('discriminator_1', axis=1)

        # Ignore tRNA-seq seeds that do not have coverage in the reference sample.
        coverage_df = coverage_df.groupby('seed_contig_name').filter(
            lambda seed_coverage_df: self.reference_sample_name in seed_coverage_df['sample_name'].values)

        self.run.info("Filtered seeds detected in reference sample", len(coverage_df['seed_contig_name'].unique()))

        # Ignore seeds that are in only one sample.
        coverage_df = coverage_df.groupby('seed_contig_name').filter(lambda seed_coverage_df: len(seed_coverage_df) > 1)

        self.run.info("Filtered seeds detected in >1 sample", len(coverage_df['seed_contig_name'].unique()))

        if len(coverage_df) == 0:
            self.info.warning(f"No seeds remain after applying the seed detection coverage threshold of {self.min_coverage}. "
                              "This threshold must be met in both the reference sample and another sample.")

        # Evaluate the anticodon wobble nucleotide in the seed.
        effective_wobble_nucleotides = []
        for decoded_aa_type, anticodon, seed_wobble_nucleotide in zip(
            trna_gene_hits_df['decoded_amino_acid'], trna_gene_hits_df['anticodon'], trna_gene_hits_df['seed_anticodon_wobble_nucleotide']):
            if decoded_aa_type == 'Ile2':
                # tRNA-Ile2 has a wobble nucleotide of lysidine in bacteria or agmatidine in
                # archaea, which are given the same decoding weight.
                effective_wobble_nucleotides.append('L')
                continue
            elif anticodon[0] == 'A':
                # Check for modification of A34 to I, which is detected as G in tRNA-seq reads.
                # tRNA-Arg-ACG and tRNA-Leu-AAG are the only bacterial tRNAs known to contain I34.
                # No archaeal tRNAs are known to contain I34. I34 has been found in 8 eukaryotic
                # tRNAs. As far as I know, the I modification is pervasive at position 34 in the
                # tRNAs that have it, so presence of G34 in the seed consensus sequence is assumed
                # to be 100% modification.
                if seed_wobble_nucleotide == 'G':
                    effective_wobble_nucleotides.append('I')
                    continue
            effective_wobble_nucleotides.append(anticodon[0])

        trna_gene_hits_df['effective_wobble_nucleotide'] = effective_wobble_nucleotides
        trna_gene_hits_df = trna_gene_hits_df.drop(['gene_sequence', 'seed_anticodon_wobble_nucleotide'], axis=1)

        # Drop duplicate rows (preserve a single row) representing hits between the same tRNA-seq
        # seed and tRNA genes with different sequences. This should only occur if the seed is a
        # partial read of the tRNA, and the genes differ beyond the 5' end of the seed. However,
        # confirm that the hits yielded the same anticodon wobble nucleotide, just in case I'm
        # missing something.
        trna_gene_hits_df = trna_gene_hits_df.drop_duplicates()
        if trna_gene_hits_df.groupby('seed_contig_name').ngroups != trna_gene_hits_df.groupby(
            ['seed_contig_name', 'effective_wobble_nucleotide']).ngroups:
            confusing_df = trna_gene_hits_df.groupby('seed_contig_name').filter(lambda seed_df: len(seed_df) > 1)
            raise ConfigError("A strange circumstance has occurred where a tRNA-seq seed linked to tRNA genes "
                              "with different sequences was found to have different effective wobble nucleotides. "
                              f"Here are the entries for the seeds in question:\n{confusing_df.to_string()}")

        seeds_df = trna_gene_hits_df.merge(coverage_df, how='inner', on='seed_contig_name')

        # Perhaps isoacceptors in the same bin could differ in their effective wobble nucleotide:
        # say one is modified to I and the other is kept A. This oddity should be noted.
        if seeds_df.groupby(['bin_name', 'decoded_amino_acid', 'anticodon']).ngroups != seeds_df.groupby(
            ['bin_name', 'decoded_amino_acid', 'anticodon', 'effective_wobble_nucleotide']).ngroups:
            confusing_df = seeds_df.groupby(['decoded_amino_acid', 'anticodon']).filter(
                lambda isoacceptor_df: len(isoacceptor_df['effective_wobble_nucleotide'].unique()) > 1)
            self.run.warning("A very strange circumstance has been found in which apparent tRNA-seq seed isoacceptors "
                             "from the same bin differ in their anticodon wobble nucleotide. "
                             "For example, in one seed, the nucleotide could be A while in the other it is modified to I. "
                             f"Here are the entries for the seeds in question: {confusing_df.to_string()}")
            seeds_df = seeds_df.groupby(['decoded_amino_acid', 'anticodon']).filter(
                lambda isoacceptor_df: len(isoacceptor_df['effective_wobble_nucleotide'].unique()) == 1)

        seeds_df = seeds_df.drop('seed_contig_name', axis=1)
        isoacceptors_df = seeds_df.groupby(
            ['bin_name', 'decoded_amino_acid', 'anticodon', 'effective_wobble_nucleotide', 'sample_name'], as_index=False).agg('sum')

        # Remove bins lacking a diversity of isoacceptors.
        prefilter_bin_names = set(isoacceptors_df['bin_name'])
        isoacceptors_df = isoacceptors_df.groupby('bin_name').filter(
            lambda bin_df: bin_df.groupby(['decoded_amino_acid', 'anticodon']).ngroups >= self.min_isoacceptors)
        removed_bin_names = set(prefilter_bin_names).difference(set(isoacceptors_df['bin_name']))
        if removed_bin_names:
            self.info.warning(f"The following bins did not meet the bin isoacceptor threshold of {self.min_isoacceptors}: {', '.join(removed_bin_names)}")
        if len(isoacceptors_df) == 0:
            self.info.warning(f"No seeds remain after applying the bin isoacceptor threshold of {self.min_isoacceptors}.")

        return isoacceptors_df


    def get_isoacceptor_abundances(self, isoacceptors_df):
        isoacceptor_abundance_dict = {}
        for bin_name, bin_df in isoacceptors_df.groupby('bin_name'):
            reference_sample_df = bin_df[bin_df['sample_name'] == self.reference_sample_name]
            reference_sample_df = reference_sample_df.set_index(['decoded_amino_acid', 'anticodon', 'effective_wobble_nucleotide'])
            nonreference_isoacceptor_df = bin_df[bin_df['sample_name'] != self.reference_sample_name]

            decoding_keys = []
            isoacceptor_abundance_ratios = []
            for row in nonreference_isoacceptor_df.itertuples(index=False):
                decoding_key = (row.decoded_amino_acid, row.anticodon, row.effective_wobble_nucleotide)
                reference_abundance = reference_sample_df.loc[decoding_key]['relative_discriminator_coverage']
                decoding_keys.append(decoding_key)
                isoacceptor_abundance_ratios.append(row.relative_discriminator_coverage / reference_abundance)
            isoacceptor_abundance_dict[bin_name] = (decoding_keys, np.array(isoacceptor_abundance_ratios))

        return isoacceptor_abundance_dict


    def consolidate_kegg_annotations(self):
        with self.genomic_contigs_db_info.load_db() as genomic_contigs_db:
            kegg_df = genomic_contigs_db.get_table_as_dataframe('gene_functions', where_clause=f'''source IN ("KOfam", "KEGG_Module", "KEGG_Class")''')
        kofam_df = kegg_df[kegg_df['source'] == 'KOfam']
        module_df = kegg_df[kegg_df['source'] == 'KEGG_Module']
        class_df = kegg_df[kegg_df['source'] == 'KEGG_Class']

        # Iterate the annotations in each block of lines for KOfam, module, and class. Relate
        # module/class to KOfam annotations. There is the same number of module and class entries.
        # Not every KOfam is part of a module/class. There is only one edge case, presumably
        # vanishingly rare, that can prevent accurate reassignment of module/class to KOfam: a gene
        # is assigned multiple KOfams, the entries for the KOfams happen to be next to each other in
        # the table (which I think would happen randomly), and one or more of the KOfams is not part
        # of a module/class. This hypothetical edge case is resolved by assigning the module/class
        # to the first occurring KOfam.
        module_iter = iter(module_df.itertuples(index=False))
        class_iter = iter(class_df.itertuples(index=False))
        try:
            module_row = next(module_iter)
            class_row = next(class_iter)
            module_gene_callers_id = module_row.gene_callers_id
        except StopIteration:
            module_row = None
            class_row = None
            module_gene_callers_id = -1

        new_rows = []
        for kofam_row in kofam_df.itertuples(index=False):
            kofam_gene_callers_id = kofam_row.gene_callers_id
            if kofam_gene_callers_id == module_gene_callers_id:
                new_rows.append([
                    kofam_gene_callers_id,
                    kofam_row.accession,
                    kofam_row.function,
                    kofam_row.e_value,
                    module_row.accession,
                    module_row.function,
                    class_row.function
                ])

                try:
                    module_row = next(module_iter)
                    class_row = next(class_iter)
                    module_gene_callers_id = module_row.gene_callers_id
                except StopIteration:
                    module_row = None
                    class_row = None
                    module_gene_callers_id = -1
            else:
                new_rows.append([
                    kofam_gene_callers_id,
                    kofam_row.accession,
                    kofam_row.function,
                    kofam_row.e_value,
                    '',
                    '',
                    ''
                ])
        new_kegg_df = pd.DataFrame(new_rows, columns=['gene_callers_id',
                                                      'kofam_accession',
                                                      'kofam_function',
                                                      'kofam_e_value',
                                                      'kegg_module_accession',
                                                      'kegg_module_function',
                                                      'kegg_class_function'])
        return new_kegg_df


    @staticmethod
    def list_sample_names(args={}, r=run):
        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        seeds_specific_txt_path = A('seeds_specific_txt')

        if seeds_specific_txt_path == None:
            raise ConfigError("To list samples in `seeds_specific_txt`, a path to this file must be provided.")
        filesnpaths.is_file_exists(seeds_specific_txt_path)

        available_sample_names = pd.read_csv(seeds_specific_txt_path, sep='\t', header=0, skiprows=[1, 2], usecols=['sample_name'])['sample_name'].unique().tolist()
        run.info_single(f"The `seeds_specific_txt` table, '{seeds_specific_txt_path}', contains the following samples: {', '.join(available_sample_names)}")
        return available_sample_names

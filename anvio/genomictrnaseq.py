# -*- coding: utf-8
# pylint: disable=line-too-long
"""tRNA-seq integration with (meta)genomics"""

import os
import argparse
import tempfile
import numpy as np
import pandas as pd

from functools import partial
from argparse import Namespace
from itertools import combinations, product

import anvio
import anvio.utils as utils
import anvio.hmmops as hmmops
import anvio.tables as tables
import anvio.fastalib as fastalib
import anvio.terminal as terminal
import anvio.constants as constants
import anvio.codonusage as codonusage
import anvio.filesnpaths as filesnpaths
import anvio.ccollections as ccollections

from anvio.dbinfo import DBInfo
from anvio.drivers.blast import BLAST
from anvio.errors import ConfigError, FilesNPathsError
from anvio.genomedescriptions import GenomeDescriptions


__author__ = "Developers of anvi'o (see AUTHORS.txt)"
__copyright__ = "Copyleft 2015-2022, the Meren Lab (http://merenlab.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__maintainer__ = "Samuel Miller"
__email__ = "samuelmiller10@gmail.com"


run = terminal.Run()
progress = terminal.Progress()
run_quiet = terminal.Run(verbose=False)

pp = terminal.pretty_print

class SeedPermuter(object):
    """Using the `go` method, generates a FASTA file of permuted seed sequences from a tRNA-seq
    contigs database."""

    default_min_nt_frequency = 0.05
    default_max_variable_positions = 5

    def __init__(self, args={}, r=run, p=progress, do_sanity_check=True):
        self.args = args
        A = lambda x: args.__dict__[x] if x in args.__dict__ else None

        self.contigs_db_path = A('trnaseq_contigs_db')
        self.modifications_txt_path = A('modifications_txt')

        self.min_nt_frequency = A('min_nt_frequency')
        if self.min_nt_frequency is None:
            self.min_nt_frequency = self.default_min_nt_frequency
        self.max_variable_positions = A('max_variable_positions')
        if self.max_variable_positions is None:
            self.max_variable_positions = self.default_max_variable_positions

        # paths for output files
        self.permuted_seeds_fasta_path = A('permuted_seeds_fasta')
        self.tmp_dir = None
        if self.permuted_seeds_fasta_path is None:
            self.tmp_dir = anvio.TMP_DIR if anvio.TMP_DIR else tempfile.gettempdir()
            self.permuted_seeds_fasta_path = os.path.join(self.tmp_dir, "permuted_seeds.fa")

        self.run = r
        self.progress = p

        if do_sanity_check:
            self.sanity_check()

        self.contigs_db_info = DBInfo(self.contigs_db_path, expecting='contigs')


    def sanity_check(self):
        """Check the feasibility of args from initialization."""
        contigs_db_info = DBInfo(self.contigs_db_path, expecting='contigs')
        if contigs_db_info.variant != 'trnaseq':
            raise ConfigError(
                f"The database at '{self.contigs_db_path}' was a '{contigs_db_info.variant}' "
                "variant, not the required 'trnaseq' variant.")

        filesnpaths.is_file_exists(self.modifications_txt_path)

        if not (0 <= self.min_nt_frequency < 1):
            raise ConfigError(
                "The specified minimum nucleotide frequency for permutation is "
                f"'{self.min_nt_frequency}', which is not in the required range [0, 1).")

        if self.max_variable_positions < 1:
            raise ConfigError(
                "The specified maximum number of variable positions in a permuted sequence is "
                f"'{self.max_variable_positions}', but it needs to be a positive integer.")

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
                permutations = [
                    str(permuted_position) + permuted_nt
                    for permuted_position, permuted_nt in zip(permuted_positions, permuted_nts)]
                seq_id = f"{seed_contig_name}|{'_'.join(permutations)}"
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
        modifications_df = pd.read_csv(self.modifications_txt_path, sep='\t', header=0,
                                       usecols=['contig_name', 'seed_position', 'A', 'C', 'G', 'T'])
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
        for position, gb_df in variable_nts_df.set_index('seed_position').groupby(
            'seed_position', sort=False):
            # Calculate nucleotide relative frequencies for each sample.
            position_variable_nts_df = gb_df.div(gb_df.sum(axis=1), axis=0)
            # Drop the column of the relative frequency of the nucleotide in the seed sequence.
            position_variable_nts_df = position_variable_nts_df.drop(seed_seq[position], axis=1)
            # Calculate the mean relative frequency of each nucleotide across samples.
            position_variable_nts_series = position_variable_nts_df.mean(axis=0)
            # Filter nucleotides by threshold relative frequency.
            position_variable_nts_series = position_variable_nts_series[
                position_variable_nts_series > self.min_nt_frequency]

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
        for num_variable_positions in range(
            1, min(len(variable_nts_dict), self.max_variable_positions) + 1):
            # Find combinations of permuted positions, e.g., with 1 permuted position, and positions
            # 8 and 32 being variable, then one combination would simply be (8, ) and the other
            # combination (32, ).
            permutation_combinations = combinations(variable_nts_dict, num_variable_positions)
            # Loop through each combination.
            for permuted_positions in permutation_combinations:
                # Find the sets of nucleotides that will be substituted into the sequence given the
                # combination of positions. Each loop generates a new permuted sequence and entry.
                for permuted_nts in product(
                    *[variable_nts_dict[position] for position in permuted_positions]):
                    permuted_seed_seq = seed_seq
                    # Loop through the positions to make the nucleotide substitutions.
                    for position, nt in zip(permuted_positions, permuted_nts):
                        permuted_seed_seq = \
                            permuted_seed_seq[: position] + nt + permuted_seed_seq[position + 1: ]
                    # In addition to the permuted sequence, record the permuted positions and
                    # substituted nucleotides.
                    permuted_seed_info.append(
                        (permuted_seed_seq, tuple(permuted_positions), tuple(permuted_nts)))

        return permuted_seed_info


class Integrator(object):
    """Using the `go` method, links tRNA-seq seeds to tRNA genes and adds this information to the
    tRNA-seq contigs database."""
    # Here are the different possible (meta)genomic sources:
    # 1. Single contigs database without bins
    # 2. Single contigs database with collection of bins
    # 3. Single contigs database with one or more specified bins
    # 4. One or more contigs databases input as "external" genomes
    # 5. "Internal" genomes (bins) from one or more contigs databases
    # 6. A combination of "internal" and "external" genomes (4 + 5)
    # Ambiguous assignment of tRNA-seq seeds to tRNA genes can be applied to all but (1).

    default_max_mismatches = 3
    blast_search_output_cols = [
        'qseqid', 'sseqid', 'mismatch', 'qstart', 'qlen', 'sstart', 'send', 'slen', 'bitscore']

    def __init__(self, args={}, r=run, rq=run_quiet, p=progress, do_sanity_check=True):
        self.args = args
        A = lambda x: args.__dict__[x] if x in args.__dict__ else None

        self.trnaseq_contigs_db_path = A('trnaseq_contigs_db')
        self.seeds_specific_txt_path = A('seeds_specific_txt')
        self.modifications_txt_path = A('modifications_txt')

        self.genomic_contigs_db_path = A('contigs_db')
        self.genomic_profile_db_path = A('profile_db')
        self.collection_name = A('collection_name')
        self.bin_id = A('bin_id')
        self.bin_ids_file = A('bin_ids_file')

        self.internal_genomes_path = A('internal_genomes')
        self.external_genomes_path = A('external_genomes')

        self.max_mismatches = A('max_mismatches')
        if self.max_mismatches is None:
            self.max_mismatches = self.default_max_mismatches
        self.full_gene = A('full_gene')
        if self.full_gene is None:
            self.full_gene = False
        self.unambiguous_genome_assignment = A('unambiguous_genome_assignment')
        if self.unambiguous_genome_assignment is None:
            self.unambiguous_genome_assignment = False

        self.permuted_seeds_fasta_path = A('permuted_seeds_fasta')
        self.blast_dir = A('blast_dir')
        if self.blast_dir is None:
            self.blast_dir = anvio.TMP_DIR if anvio.TMP_DIR else tempfile.gettempdir()

        self.num_threads = A('num_threads') or anvio.K('num_threads')['default']
        self.remove_previous_matches = A('remove_previous_matches')

        self.trna_genes_fasta_path = os.path.join(self.blast_dir, 'trna_genes.fa')

        self.run = r
        self.run_quiet = rq
        self.progress = p

        self.trnaseq_contigs_db_info = DBInfo(self.trnaseq_contigs_db_path, expecting='contigs')

        # Store information on accessing (meta)genomes.
        self.genome_info_dict = {}

        # If `contigs_db` was provided, then `internal_genomes` and `external_genomes` should not
        # also have been provided, which is checked later in `sanity_check`.
        if self.genomic_contigs_db_path:
            contigs_db_info = DBInfo(self.genomic_contigs_db_path, expecting='contigs')
        else:
            contigs_db_info = None
        if self.genomic_profile_db_path:
            profile_db_info = DBInfo(self.genomic_profile_db_path, expecting='profile')
        else:
            profile_db_info = None
        if self.bin_ids_path:
            filesnpaths.is_file_plain_text(self.bin_ids_path)
            bin_ids = []
            with open(self.bin_ids_path) as bin_ids_file:
                for line in bin_ids_file:
                    bin_ids.append(line.rstrip())
        else:
            bin_ids = None

        if self.bin_id:
            collect = ccollections.Collections()
            collect.populate_collections_dict(self.genomic_profile_db_path)
            collect.is_bin_in_collection(bin_id)
            self.genome_info_dict[bin_id] = genome_info = {}
            genome_info['contigs_db_info'] = contigs_db_info
            genome_info['profile_db_info'] = profile_db_info
            genome_info['collection_name'] = self.collection_name
            genome_info['bin_id'] = self.bin_id
        elif self.bin_ids_path:
            collect = ccollections.Collections()
            collect.populate_collections_dict(self.genomic_profile_db_path)
            for bin_id in bin_ids:
                collect.is_bin_in_collection(bin_id)
                self.genome_info_dict[bin_id] = genome_info = {}
                genome_info['contigs_db_info'] = contigs_db_info
                genome_info['profile_db_info'] = profile_db_info
                genome_info['collection_name'] = self.collection_name
                genome_info['bin_id'] = bin_id
        elif self.collection_name:
            # There is a collection of internal genomes.
            collect = ccollections.Collections()
            collect.populate_collections_dict(self.genomic_profile_db_path)
            for bin_id in collect.get_bins_info_dict():
                self.genome_info_dict[bin_id] = genome_info = {}
                genome_info['contigs_db_info'] = contigs_db_info
                genome_info['profile_db_info'] = profile_db_info
                genome_info['collection_name'] = self.collection_name
                genome_info['bin_id'] = bin_id
        elif self.genomic_contigs_db_path:
            # The contigs database represents a genome, like an external genome.
            self.genome_info_dict[contigs_db_info.project_name] = genome_info = {}
            genome_info['contigs_db_info'] = contigs_db_info
            genome_info['profile_db_info'] = None
            genome_info['collection_name'] = None
            genome_info['bin_id'] = None

        nonunique_genome_names = []
        if self.internal_genomes_path or self.external_genomes_path:
            descriptions = GenomeDescriptions(args, run=self.run_quiet, progress=self.progress)
            descriptions.load_genomes_descriptions(init=False)

            for genome_name, genome_dict in descriptions.internal_genomes_dict.items():
                if genome_name in self.genome_info_dict:
                    nonunique_genome_names.append(genome_name)

                self.genome_info_dict[genome_name] = genome_info = {}
                genome_info['contigs_db_info'] = DBInfo(
                    genome_dict['contigs_db_path'], expecting='contigs')
                if genome_dict['profile_db_path']:
                    genome_info['profile_db_info'] = DBInfo(
                        genome_dict['profile_db_path'], expecting='profile')
                else:
                    genome_info['profile_db_info'] = None
                genome_info['collection_name'] = genome_dict['collection_id']
                genome_info['bin_id'] = genome_dict['bin_id']

            for genome_name, genome_dict in descriptions.external_genomes_dict.items():
                self.genome_info_dict[genome_name] = genome_info = {}
                genome_info['contigs_db_info'] = DBInfo(
                    genome_dict['contigs_db_path'], expecting='contigs')
                genome_info['profile_db_info'] = None
                genome_info['collection_name'] = None
                genome_info['bin_id'] = None

        if nonunique_genome_names:
            raise ConfigError(
                "Names must be unique to internal and external genomes, but the following were "
                f"not: {', '.join(nonunique_genome_names)}")

        if do_sanity_check:
            self.sanity_check()


    def sanity_check(self, check_permuted_seeds_fasta=False):
        """Check the feasibility of args from initialization."""
        if self.trnaseq_contigs_db_info.variant != 'trnaseq':
            raise ConfigError(
                f"The database at '{self.trnaseq_contigs_db_path}' was a "
                f"'{self.trnaseq_contigs_db_info.variant}' variant, not the required 'trnaseq' "
                "variant.")

        # Check that the table of seed/gene hits has the correct number of rows since the structure
        # of the table changed in the development of this program without a version upgrade of the
        # database.
        trnaseq_contigs_db = self.trnaseq_contigs_db_info.load_db()
        if (len(trnaseq_contigs_db.get_table_column_types(tables.trna_gene_hits_table_name)) !=
            len(tables.trna_gene_hits_table_types)):
            if not self.remove_previous_matches:
                raise ConfigError(
                    "The table of seed/gene matches in the tRNA-seq contigs database, "
                    f"'{self.trnaseq_contigs_db_path}', does not have the proper number of "
                    "columns, probably because it was created with an obsolete version of the "
                    "program. The table can be deleted and recreated with the proper columns by "
                    "running the program with `remove_previous_matches`.")

        # Existing seed/gene hits must be willfully overwritten or appended to.
        hit_count = trnaseq_contigs_db.get_row_counts_from_table(
            tables.trna_gene_hits_table_name)
        if hit_count:
            self.run.info(
                "Preexisting tRNA-seq seed/tRNA gene hits in the tRNA-seq contigs db", hit_count)
            if self.unambiguous_genome_assignment and not self.remove_previous_matches:
                raise ConfigError(
                    "The seeds from the tRNA-seq contigs database at "
                    f"'{self.trnaseq_contigs_db_path}' have already been associated with tRNA "
                    "genes from one or more (meta)genomic contigs databases. Either use "
                    "`remove_previous_matches` to clear existing matches, or append to existing "
                    "matches by not using `unambiguous_genome_assignments`.")
        trnaseq_contigs_db.disconnect()

        # The tRNA-seq contigs db version must be up-to-date.
        required_version = utils.get_required_version_for_db(self.trnaseq_contigs_db_path)
        if str(self.trnaseq_contigs_db_info.version) != required_version:
            raise ConfigError(
                f"The database at '{self.trnaseq_contigs_db_path}' is outdated (this database is "
                f"v{self.trnaseq_contigs_db_info.version} and your anvi'o installation wants to "
                f"work with v{required_version}). You can migrate your database without losing any "
                "data using the program `anvi-migrate` with either of the flags "
                "`--migrate-dbs-safely` or `--migrate-dbs-quickly`.")

        # Right now there are no specific checks here on the format of these tables.
        filesnpaths.is_file_exists(self.seeds_specific_txt_path)
        filesnpaths.is_file_exists(self.modifications_txt_path)

        # Do basic checks of the combinations of (meta)genomic input arguments.
        if (self.genomic_contigs_db_path and
            (self.internal_genomes_path or self.external_genomes_path)):
            raise ConfigError("`internal_genomes` and `external_genomes` cannot be used with "
                              "`contigs_db`.")

        if ((self.genomic_profile_db_path or self.collection_name) and
            not (self.genomic_contigs_db_path and
                 self.genomic_profile_db_path and
                 self.collection_name)):
            raise ConfigError("A collection must be provided using `contigs_db`, `profile_db`, and "
                              "`collection_name`.")

        if (self.bin_id and
            not (self.genomic_contigs_db_path and
                 self.genomic_profile_db_path and
                 self.collection_name)):
            raise ConfigError("A specific bin provided with `bin_id` also requires `contigs_db`, "
                              "`profile_db`, and `collection_name`.")

        # Prevent a confused user from providing a tRNA-seq contigs database in lieu of a
        # (meta)genomic contigs database.
        unrecognized = []
        for name, genome_info in self.genome_info_dict.items():
            if genome_info['contigs_db_info'].variant != 'unknown':
                unrecognized.append(name)
        if unrecognized:
            if self.genomic_contigs_db_path:
                raise ConfigError(
                    f"The purported (meta)genomic contigs database, '{unrecognized[0]}', was not "
                    "recognized as such. A proper database (technically of the 'unknown' variant) "
                    "should be generated by `anvi-gen-contigs-database`.")
            else:
                raise ConfigError(
                    "The purported (meta)genomic contigs databases for the following genomes were "
                    "not recognized as such. A proper database (technically of the 'unknown' "
                    "variant) should be generated by `anvi-gen-contigs-database`. "
                    f"{', '.join(unrecognized)}")

        # Check that there are tRNA genes annotated in the (meta)genomes.
        unannotated = []
        for name, genome_info in self.genome_info_dict.items():
            if 'Transfer_RNAs' not in hmmops.SequencesForHMMHits(
                genome_info['contigs_db_info'].path).hmm_hits_info:
                unannotated.append(name)
        if unannotated:
            if self.genomic_contigs_db_path:
                raise ConfigError(
                    "It appears that tRNA genes have not been annotated in the (meta)genomic "
                    f"contigs database, '{unannotated[0]}'. Please run `anvi-scan-trnas` on the "
                    "database and try again (this same error will arise if no tRNA genes are "
                    "found).")
            else:
                raise ConfigError(
                    "It appears that tRNA genes have not been annotated in the following "
                    "genomes. Please run `anvi-scan-trnas` on the contigs databases and try again "
                    "(this same error will arise if no tRNA genes are found). "
                    f"{', '.join(unannotated)}.")

        # Check that profile databases correspond to (meta)genomic contigs databases.
        incompatible = []
        for name, genome_info in self.genome_info_dict.items():
            if genome_info['profile_db_info']:
                try:
                    utils.is_profile_db_and_contigs_db_compatible(
                        genome_info['profile_db_info'].path, genome_info['contigs_db_info'].path)
                except ConfigError:
                    incompatible.append(name)
        if incompatible:
            if self.genomic_contigs_db_path:
                raise ConfigError(
                    f"The (meta)genomic contigs database, '{incompatible[0]}', is not compatible "
                    f"with the provided profile database, {genome_info['profile_db_info'].path}. "
                    "In fact, the profile database was not generated from the contigs database.")
            else:
                raise ConfigError(
                    "The contigs databases for the following genomes are not compatible with the "
                    "corresponding profile databases. In fact, the profile databases were not "
                    f"generated from the contigs databases. {', '.join(incompatible)}")

        # Check putative collections.
        unrecognized = []
        for name, genome_info in self.genome_info_dict.items():
            if genome_info['profile_db_info']:
                collections = ccollections.Collections()
                collections.populate_collections_dict(genome_info['profile_db_info'].path)
                if genome_info['collection_name'] not in collections.collections_dict:
                    unrecognized.append(name)
        if unrecognized:
            if self.genomic_contigs_db_path:
                raise ConfigError(
                    f"The profile database, '{self.genomic_profile_db_path}', does not contain the "
                    f"requested collection, '{self.collection_name}'.")
            else:
                raise ConfigError(
                    "The profile databases do not contain the requested collections for the "
                    f"following genomes. {', '.join(unrecognized)}")

        if self.max_mismatches < 0:
            raise ConfigError(
                "The maximum number of mismatches allowed in a seed-gene alignment must be a "
                f"non-negative integer, not the provided value of {self.max_mismatches}")

        filesnpaths.is_output_dir_writable(self.blast_dir)

        # Ignore this sanity check when using `genomictrnaseq.Permuter` FASTA output.
        if check_permuted_seeds_fasta:
            filesnpaths.is_file_fasta_formatted(self.permuted_seeds_fasta_path)

        if self.num_threads < 1:
            raise ConfigError("The number of threads (used by BLAST) must be a positive integer, "
                              f"not the provided value of {self.num_threads}")


    def go(self):
        """Link tRNA-seq seeds to tRNA genes, adding this information to the tRNA-seq contigs
        database."""
        trna_gene_seq_dict = self.write_trna_genes_fasta()
        self.blast()
        hits_df = self.filter_hits(trna_gene_seq_dict)
        unmodified_nt_df = self.find_unmodified_nucleotides(hits_df)
        self.update_trnaseq_contigs_database(hits_df, unmodified_nt_df)


    def write_trna_genes_fasta(self):
        """
        Write all tRNA gene sequences from input contigs databases to a FASTA file.

        Returns
        =======
        trna_gene_seq_dict : dict
            tRNA gene sequences keyed by tuple of contigs database name and gene callers ID.
        """
        trna_genes_fasta = open(self.trna_genes_fasta_path, 'w')

        # Unfortunately, the full subject (tRNA gene) sequence cannot be reported in the BLAST
        # output table, but these sequences are needed to filter alignments, so they are filed in a
        # dictionary.
        trna_gene_seq_dict = {}

        # Get the unique set of input contigs databases.
        contigs_db_paths = set()
        contigs_db_infos = []
        for genome_info in self.genome_info_dict.values():
            if genome_info['contigs_db_info'].path in contigs_db_paths:
                continue
            contigs_db_infos.append(genome_info['contigs_db_info'])

        for contigs_db_info in contigs_db_infos:
            contigs_db_project_name = contigs_db_info.project_name
            contigs_db_hash = contigs_db_info.hash

            trna_gene_info = hmmops.SequencesForHMMHits(
                contigs_db_info.path, sources=set(['Transfer_RNAs']))

            # Split names from the database are needed here to recover the tRNA sequence strings.
            with contigs_db_info.load_db() as contigs_db:
                splits_dict = {
                    contigs_db_hash: list(contigs_db.smart_get(
                        tables.splits_info_table_name, 'split').keys())}
            hmm_seqs_dict = trna_gene_info.get_sequences_dict_for_hmm_hits_in_splits(splits_dict)

            for gene_id, gene_entry in hmm_seqs_dict.items():
                seq_string = hmm_seqs_dict[gene_id]['sequence']

                # Record both the project name and hash of the contigs database: the project name is
                # not guaranteed to be unique.
                header = (f"{contigs_db_project_name}|"
                          f"{contigs_db_hash}|"
                          f"{gene_entry['contig']}|"
                          f"{gene_entry['gene_callers_id']}|"
                          f"{gene_entry['gene_name']}|"
                          f"{gene_entry['e_value']}|"
                          f"{gene_entry['start']}|"
                          f"{gene_entry['stop']}")
                trna_genes_fasta.write(f">{header}\n")
                trna_genes_fasta.write(f"{seq_string}\n")

                trna_gene_seq_dict[(contigs_db_hash, gene_entry['gene_callers_id'])] = seq_string

        trna_genes_fasta.close()

        return trna_gene_seq_dict


    def blast(self):
        """Align tRNA-seq seeds/permuted seeds to tRNA genes."""
        blast = BLAST(self.permuted_seeds_fasta_path,
                      self.trna_genes_fasta_path,
                      search_program='blastn',
                      run=self.run,
                      progress=self.progress,
                      num_threads=self.num_threads)
        blast.tmp_dir = self.blast_dir
        blast.search_output_path = os.path.join(self.blast_dir, 'blast-search-results.txt')
        blast.log_file_path = os.path.join(self.blast_dir, 'blast-log.txt')
        blast.additional_params_for_blast = "-ungapped"
        blast.makedb(dbtype='nucl')
        blast.blast(outputfmt='6 ' + ' '.join(self.blast_search_output_cols))


    def filter_hits(self, trna_gene_seq_dict):
        """
        Confidently associate tRNA-seq seeds with tRNA genes, filtering BLAST alignments of
        seeds/permuted seeds to genes.

        Parameters
        ==========
        trna_gene_seq_dict : dict
            tRNA gene sequences used as subjects in BLAST search keyed by tuple of contigs database
            name and gene callers ID.

        Returns
        =======
        hits_df : pandas.core.frame.DataFrame
            Each row contains a selected hit between a seed, which may be permuted, and tRNA gene.
        """
        # Load BLAST output table.
        search_output_path = os.path.join(self.blast_dir, 'blast-search-results.txt')
        hits_df = pd.read_csv(
            search_output_path, sep='\t', header=None, names=self.blast_search_output_cols)

        # Discard alignments with too many mismatches.
        hits_df = hits_df[hits_df['mismatch'] <= self.max_mismatches]
        if self.full_gene:
            # Discard alignments that do not start at the beginning of the gene.
            hits_df = hits_df[hits_df['sstart'] == 1]

        # Discard (enigmatic) indistinguishable, duplicate hits if they exist.
        hits_df = hits_df.drop_duplicates()

        # Parse seed IDs and permutation info.
        hits_df[['seed_contig_name', 'seed_permutation']] = \
            hits_df['qseqid'].str.split('|', expand=True)
        hits_df['seed_permutation'] = hits_df['seed_permutation'].fillna('')
        hits_df = hits_df.drop('qseqid', axis=1)

        # Extract information on each hit.
        contigs_db_project_names = []
        contigs_db_hashes = []
        gene_contig_names = []
        gene_callers_ids = []
        decoded_amino_acids = []
        anticodons = []
        trnascan_scores = []
        gene_starts = []
        gene_stops = []
        gene_sequences = []
        for sseqid in hits_df['sseqid']:
            split_sseqid = sseqid.split('|')
            contigs_db_project_name = split_sseqid[0]
            contigs_db_project_names.append(contigs_db_project_name)
            contigs_db_hash = split_sseqid[1]
            contigs_db_hashes.append(contigs_db_hash)
            gene_contig_names.append(split_sseqid[2])
            gene_callers_id = int(split_sseqid[3])
            gene_callers_ids.append(gene_callers_id)
            gene_name = split_sseqid[4]
            decoded_amino_acid, anticodon = gene_name.split('_')[: 2]
            decoded_amino_acids.append(decoded_amino_acid)
            anticodons.append(anticodon)
            trnascan_scores.append(float(split_sseqid[5]))
            gene_starts.append(int(split_sseqid[6]))
            gene_stops.append(int(split_sseqid[7]))
            gene_sequences.append(trna_gene_seq_dict[(contigs_db_hash, gene_callers_id)])
        hits_df['contigs_db_project_name'] = contigs_db_project_names
        hits_df['contigs_db_hash'] = contigs_db_hashes
        hits_df['gene_contig_name'] = gene_contig_names
        hits_df['decoded_amino_acid'] = decoded_amino_acids
        hits_df['anticodon'] = anticodons
        hits_df['trnascan_score'] = trnascan_scores
        hits_df['gene_callers_id'] = gene_callers_ids
        hits_df['gene_start_in_contig'] = gene_starts
        hits_df['gene_stop_in_contig'] = gene_stops
        hits_df['gene_sequence'] = gene_sequences

        # Filter individual alignments.
        retained_indices = []
        for (index,
             decoded_amino_acid,
             seed_alignment_start,
             seed_length,
             gene_alignment_start,
             gene_alignment_end,
             gene_sequence,
             gene_length) in zip(hits_df.index,
                                 hits_df['decoded_amino_acid'],
                                 hits_df['qstart'],
                                 hits_df['qlen'],
                                 hits_df['sstart'],
                                 hits_df['send'],
                                 hits_df['gene_sequence'],
                                 hits_df['slen']):
            if ((gene_length - gene_alignment_end == 0) or
                ((gene_length - gene_alignment_end == 3) and gene_sequence[-3: ] == 'CCA')):
                # The alignment ends at the end of the gene or just short of a 3'-CCA acceptor in
                # the gene (the seed/permuted seed should never contain the 3'-CCA acceptor).
                if ((seed_alignment_start == 1) and
                    (gene_alignment_end - gene_alignment_start == seed_length - 1)):
                    # The alignment starts at the beginning of the seed/permuted seed and spans the
                    # entire query.
                    retained_indices.append(index)
                elif ((decoded_amino_acid == 'His') and
                      (seed_alignment_start == 2) and
                      (gene_alignment_end - gene_alignment_start == seed_length - 2)):
                    # The alignment starts at the second position of the tRNA-His seed/permuted seed
                    # sequence, which has a post-transcriptional G at the 5' end, and spans the
                    # remaining length of the query.
                    retained_indices.append(index)
        hits_df = hits_df.loc[retained_indices]

        # Retain each seed/permuted seed's top-scoring hits.
        hits_df = hits_df[
            hits_df.groupby('seed_contig_name')['bitscore'].transform('max') == hits_df['bitscore']]

        ##################################################
        # Now add information on internal genomes of interest that contain the tRNA genes. Note that
        # external genomes derived from the same metagenome can contain the same tRNA gene: hits to
        # the same gene in different external genomes are distinguished by contigs database project
        # name/hash.
        # Create a dictionary, `contig_bin_dict`, mapping the names of contigs bearing tRNA genes to
        # bin info.
        contig_bin_dict = {}
        for genome_info in self.genome_info_dict.values():
            if not genome_info['bin_id'] and not genome_info['collection_name']:
                continue

            args = argparse.ArgumentParser()
            if genome_info['bin_id']:
                args.contigs_db = genome_info['contigs_db_info'].path
                args.profile_db = genome_info['profile_db_info'].path
                args.collection_name = collection_name = genome_info['collection_name']
                args.bin_id = genome_info['bin_id']
                search_for_bin_of_interest = True
            elif genome_info['collection_name']:
                # A single collection was supplied in the input arguments, `self.collection_name`.
                args.contigs_db = genome_info['contigs_db_info'].path
                args.profile_db = genome_info['profile_db_info'].path
                args.collection_name = collection_name = genome_info['collection_name']
                search_for_bin_of_interest = True
            else:
                search_for_bin_of_interest = False

            # Note that the same contig may be in different bins, thus the list values of the dict.
            contigs_db_hash = genome_info['contigs_db_info'].hash
            profile_db_sample_id = genome_info['profile_db_info'].sample_id
            if search_for_bin_of_interest:
                bin_contig_names_dict = ccollections.GetSplitNamesInBins(args).get_dict()
                for bin_id, split_names in bin_contig_names_dict.items():
                    for split_name in split_names:
                        contig_name = split_name.split('_split_')[0]
                        try:
                            contig_bin_dict[(contigs_db_hash, contig_name)].append(
                                (profile_db_sample_id, collection_name, bin_id))
                        except KeyError:
                            contig_bin_dict[(contigs_db_hash, contig_name)] = [
                                (profile_db_sample_id, collection_name, bin_id)]

        # Make a table of the membership of gene-bearing contigs in contigs databases/bins.
        contig_bin_rows = []
        for row in hits_df[
            ['contigs_db_hash', 'gene_contig_name']].drop_duplicates().itertuples(index=False):
            contigs_db_hash = row.contigs_db_hash
            contig_name = row.gene_contig_name
            try:
                bin_info = contig_bin_dict[(contigs_db_hash, contig_name)]
            except KeyError:
                # The contig is not binned.
                contig_bin_rows.append([contigs_db_hash, contig_name, '', '', ''])
                continue
            for profile_db_sample_id, collection_name, bin_id in bin_info:
                # Record each bin containing the contig.
                contig_bin_rows.append(
                    [contigs_db_hash, contig_name, profile_db_sample_id, collection_name, bin_id])
        contig_bin_df = pd.DataFrame(contig_bin_rows,
                                     columns=['contigs_db_hash',
                                              'gene_contig_name',
                                              'profile_db_sample_id',
                                              'collection_name',
                                              'bin_id'])
        # Merge the table of seed/gene hits with the new table of bin membership, multiplying each
        # row per hit by each bin containing the gene.
        hits_df = hits_df.merge(contig_bin_df, on=['contigs_db_hash', 'gene_contig_name'])

        # For the sake of clarity, here is what happens to each of the different possible
        # (meta)genomic sources with "unambiguous" tRNA gene assignment.
        # 1. Single contigs database without bins: The existence of genomes is not assumed, so no
        #    hits are disregarded.
        # 2. Single contigs database with collection of bins: disregard seeds with equally strong
        #    hits that are not confined to a single bin.
        # 3. Single contigs database with specified bin: disregard seeds with equally strong hits
        #    that are not confined to the bin.
        # 4. One or more contigs databases input as "external" genomes: disregard seeds with equally
        #    strong hits that are not confined to a single contigs database.
        # 5. "Internal" genomes (bins) from one or more contigs databases: disregard seeds with
        #    equally strong hits that are not confined to a single bin.
        # 6. A combination of "internal" and "external" genomes (4 + 5): disregard seeds with
        #    equally strong hits that are not confined to a single internal genome bin or external
        #    genome contigs database.
        if self.genomic_contigs_db_path and not self.collection_name and not self.bin_id: # (1)
            is_simple_contigs_db_input = True
        else:
            is_simple_contigs_db_input = False
        if self.unambiguous_genome_assignment and not is_simple_contigs_db_input:
            # Drop hits to genes in multiple bins: partly takes care of (2), (5), and (6).
            hits_df = hits_df.groupby('seed_contig_name').filter(
                lambda seed_df: len(seed_df[['contigs_db_hash',
                                             'profile_db_sample_id',
                                             'collection_name',
                                             'bin_id']].drop_duplicates()) == 1)
            # Drop hits to genes in multiple contigs databases: takes care of (4), partly takes care
            # of (5) and (6).
            hits_df = hits_df.groupby('seed_contig_name').filter(
                lambda seed_df: len(seed_df['contigs_db_hash'].drop_duplicates()) == 1)
            # Drop hits to genes in unbinned contigs: takes care of (3), finishes taking care of
            # (2), (5), and (6).
            if (self.collection_name or
                (self.internal_genomes_path and not self.external_genomes_path)):
                hits_df = hits_df[hits_df['bin_id'] != '']
            elif self.internal_genomes_path and self.external_genomes_path:
                external_genome_contigs_db_hashes = [
                    genome_info['contigs_db_hash'] for genome_info in self.genome_info_dict.values()
                    if genome_info['bin_id'] is None]
                retained_index = []
                for key, contigs_db_hash, bin_id in zip(hits_df.index, hits_df['contigs_db_hash'], hits_df['bin_id']):
                    if contigs_db_hash in external_genome_contigs_db_hashes:
                        retained_index.append(key)
                        continue
                    if bin_id is not None:
                        retained_index.append(key)
                hits_df = hits_df.loc[retained_index]
        ##################################################

        # Multiple permutations of the same seed may be retained after filtering by score. There are
        # two and possibly more ways that this can occur. (1) The unmodified nucleotide at a
        # modified position has a very low frequency and so was not used in the permuted sequences.
        # The permuted sequences, none of which contain the correct nucleotide, mismatch this
        # nucleotide in genes equally well. (2) A permutation is introduced at a predicted
        # modification position that is actually a nucleotide variant, and different versions of the
        # variant occur in different (meta)genomic contigs. The following procedure resolves both of
        # these possibilities, with the last step being the one that resolves the first possibility.
        # (1) Choose the permutation hitting the greatest number of genes. (2) If not resolved,
        # choose the permutation with the fewest permuted positions. (3) If not resolved, break the
        # tie by choosing the first permutation in the table, which will favor permutations toward
        # the 5' end.
        are_permutations_unresolved = True
        if hits_df.groupby('seed_contig_name').ngroups == hits_df.groupby(
            ['seed_contig_name', 'seed_permutation']).ngroups:
            are_permutations_unresolved = False
        if are_permutations_unresolved: # (1)
            hits_df['count'] = hits_df.groupby(
                ['seed_contig_name', 'seed_permutation'], as_index=False)[
                    'seed_contig_name'].transform(len)
            hits_df = hits_df[
                hits_df['count'] == hits_df.groupby('seed_contig_name')['count'].transform('max')]
            hits_df = hits_df.drop('count', axis=1)
            if hits_df.groupby('seed_contig_name').ngroups == hits_df.groupby(
                ['seed_contig_name', 'seed_permutation']).ngroups:
                are_permutations_unresolved = False
        if are_permutations_unresolved: # (2)
            hits_df['num_permuted_positions'] = hits_df[
                'seed_permutation'].apply(lambda p: p.count('_'))
            hits_df = hits_df[
                hits_df['num_permuted_positions'] == hits_df.groupby('seed_contig_name')[
                    'num_permuted_positions'].transform('min')]
            hits_df = hits_df.drop('num_permuted_positions', axis=1)
            if hits_df.groupby('seed_contig_name').ngroups == hits_df.groupby(
                ['seed_contig_name', 'seed_permutation']).ngroups:
                are_permutations_unresolved = False
        if are_permutations_unresolved: # (3)
            hits_df = hits_df[hits_df['seed_permutation'] == hits_df.groupby('seed_contig_name')[
                'seed_permutation'].transform('first')]

        ##################################################
        # Issues can arise in the selection of accurate seeds matching genes, especially in
        # relatively deeply sequenced samples with high coverages, such as tRNA-seq libraries of
        # pure cultures.
        # I. Seeds can be artifacts of the anvi'o de novo tRNA-seq workflow. `anvi-merge-trnaseq`
        # reports up to the number of seeds set by the user. If the user asks for 1,000 seeds from a
        # bacterial isolate experiment, then ~25-50 of these seeds will be true tRNA sequences and
        # up to ~950-975 will be artifacts (containing unaccounted modification-induced indels,
        # nontemplated nucleotides, sequence errors, etc., typically at low frequency), that could
        # not be resolved as non-tRNA by the tRNA-seq workflow. To remove these artifact seeds, hits
        # to the same gene are sorted by number of mismatches in the alignment and seed abundance,
        # and only the lowest mismatch/highest seed abundance hit is retained. Seed abundance is
        # taken as the average of relative abundance in each sample based on 3' (discriminator
        # nucleotide) coverage of the seed. For example, if there are two tRNA-seq samples in the
        # experiment, and two seeds hit the same gene each with one mismatch, but one seed has
        # relative 3' abundances of 0.02 and 0.03 in the two samples and the other seed has
        # abundances of 0.0006 and 0.00008, then the hit to the former seed will be the only one
        # retained for this gene.
        # II. The selection of the lowest mismatch seed (see the previous section) can create an
        # unintended side effect. The A -> I34 wobble position modification is typically nearly 100%
        # complete. I is detected as G in tRNA-seq reads, so the correct seed matching the gene
        # should have at least one alignment mismatch, G/A at position 34. In very deeply sequenced
        # samples, however, a seed with A34 can sometimes be detected from rare tRNA molecules
        # lacking the modification. `anvi-merge-trnaseq` would not merge the A34 and G34 seeds due
        # to the absence of a third or fourth mutated nucleotide at position 34 and the miniscule
        # frequency of A. Therefore, in this case, the lowest mismatch seed should not be selected;
        # rather, select the lowest mismatch seed with G34. The algorithm also confirms that the G34
        # seed is >10x more abundant than the A34 seed.
        coverage_df = pd.read_csv(self.seeds_specific_txt_path, sep='\t', header=0, skiprows=[1, 2])
        coverage_df = coverage_df[['contig_name', 'sample_name', 'relative_discriminator_coverage']]
        coverage_df = coverage_df.rename({'contig_name': 'seed_contig_name'}, axis=1)
        coverage_df = coverage_df[
            coverage_df['seed_contig_name'].isin(hits_df['seed_contig_name'].unique())]
        coverage_df = hits_df.merge(coverage_df, on='seed_contig_name')

        # Isolate the nucleotide at wobble position 34 in seeds that hit tRNAs with A34.
        trnaseq_contigs_db = self.trnaseq_contigs_db_info.load_db()
        # Convert seed contigs names to gene callers IDs.
        seed_contig_names_string = ','.join(
            ['"%s"' % seed_contig_name for seed_contig_name in
             hits_df[hits_df['anticodon'].str[0] == 'A']['seed_contig_name'].unique()])
        contigs_where_clause = f'''contig IN ({seed_contig_names_string})'''
        seed_id_df = trnaseq_contigs_db.get_table_as_dataframe(
            'genes_in_contigs',
            columns_of_interest=['gene_callers_id', 'contig'],
            where_clause=contigs_where_clause)
        seed_id_df = seed_id_df.rename(
            {'gene_callers_id': 'seed_gene_callers_id', 'contig': 'seed_contig_name'}, axis=1)
        # Find the index of position 34 in each of the seed sequences.
        seed_ids_string = ','.join(['"%s"' % seed_gene_callers_id for seed_gene_callers_id in
                                    seed_id_df['seed_gene_callers_id'].unique()])
        ids_where_clause = f'''gene_callers_id IN ({seed_ids_string})'''
        seed_wobble_df = trnaseq_contigs_db.get_table_as_dataframe(
            'trna_feature',
            columns_of_interest=['gene_callers_id', 'anticodon_loop_start'],
            where_clause=ids_where_clause)
        seed_wobble_df = seed_wobble_df.rename(
            {'gene_callers_id': 'seed_gene_callers_id',
             'anticodon_loop_start': 'seed_anticodon_loop_start'}, axis=1)
        seed_wobble_df['seed_anticodon_start'] = seed_wobble_df['seed_anticodon_loop_start'] + 2
        seed_wobble_df = seed_wobble_df.drop('seed_anticodon_loop_start', axis=1)
        seed_wobble_df = seed_id_df.merge(seed_wobble_df, on='seed_gene_callers_id')
        # Get the seed consensus sequence strings.
        seed_consensus_sequence_df = trnaseq_contigs_db.get_table_as_dataframe(
            'contig_sequences', where_clause=contigs_where_clause)
        trnaseq_contigs_db.disconnect()
        seed_consensus_sequence_df = seed_consensus_sequence_df.rename(
            {'contig': 'seed_contig_name', 'sequence': 'seed_sequence'}, axis=1)
        seed_wobble_df = seed_wobble_df.merge(seed_consensus_sequence_df, on='seed_contig_name')
        # Find the nucleotides read at wobble position 34 in the seeds.
        anticodon_wobble_nucleotides = []
        for anticodon_start, seed_consensus_sequence in zip(
            seed_wobble_df['seed_anticodon_start'], seed_wobble_df['seed_sequence']):
            anticodon_wobble_nucleotides.append(seed_consensus_sequence[anticodon_start])
        seed_wobble_df['seed_anticodon_wobble_nucleotide'] = anticodon_wobble_nucleotides

        def filter_multiple_hits_to_gene(gene_df): # inner function used in groupby apply
            min_mismatch_df = gene_df[gene_df['mismatch'] == gene_df['mismatch'].min()]
            if min_mismatch_df['seed_contig_name'].nunique() > 1:
                min_mismatch_df = min_mismatch_df[
                    min_mismatch_df['seed_contig_name'] == min_mismatch_df.groupby(
                        'seed_contig_name')['relative_discriminator_coverage'].mean().sort_values(
                            ascending=False).index[0]]
            if gene_df['anticodon'].iloc[0][0] != 'A':
                return min_mismatch_df
            else:
                # Case II: Address possibility that seed with I34 is neglected.
                max_coverage_seed_contig_name = gene_df.groupby('seed_contig_name')[
                        'relative_discriminator_coverage'].mean().sort_values(
                            ascending=False).index[0]
                min_mismatch_seed_contig_name = min_mismatch_df['seed_contig_name'].iloc[0]
                if min_mismatch_seed_contig_name == max_coverage_seed_contig_name:
                    # The seed with the fewest mismatches also has the highest average discriminator
                    # coverage.
                    return min_mismatch_df
                wobble_df = gene_df.merge(seed_wobble_df, how='left', on='seed_contig_name')
                min_mismatch_seed_wobble_nucleotide = wobble_df[
                    wobble_df['seed_contig_name'] == min_mismatch_seed_contig_name][
                        'anticodon'].iloc[0][0]
                if min_mismatch_seed_wobble_nucleotide != 'A':
                    # The seed with the fewest mismatches (and which does not have the highest
                    # average discriminator coverage) does not have an A at position 34. A seed is
                    # not matched to the gene.
                    return pd.DataFrame().reindex_like(min_mismatch_df)
                wobble_df = wobble_df[wobble_df['seed_anticodon_wobble_nucleotide'] == 'G']
                wobble_df = wobble_df[wobble_df['mismatch'] == wobble_df['mismatch'].min()]
                max_coverage_G34_seed_contig_name = wobble_df.groupby('seed_contig_name')[
                    'relative_discriminator_coverage'].mean().sort_values(ascending=False).index[0]
                if (wobble_df[wobble_df['seed_contig_name'] == max_coverage_G34_seed_contig_name][
                    'relative_discriminator_coverage'].mean() >= 10 * min_mismatch_df[
                        'relative_discriminator_coverage'].mean()):
                    # The seed with the fewest mismatches and a G at position 34 has an average
                    # discriminator coverage at least an order of magnitude higher than the the seed
                    # with the fewest mismatches and an A at position 34. Match the former seed to
                    # the gene instead of the latter.
                    return gene_df[gene_df['seed_contig_name'] == max_coverage_G34_seed_contig_name]
                else:
                    # The seed with the fewest mismatches and a G at position 34 does not have an
                    # average discriminator coverage at least an order of magnitude higher than the
                    # the seed with the fewest mismatches and an A at position 34. A seed is not
                    # matched to the gene.
                    return pd.DataFrame().reindex_like(min_mismatch_df)

        coverage_df = coverage_df.groupby(
            'gene_callers_id', group_keys=False).apply(filter_multiple_hits_to_gene)
        hits_df = hits_df[
            hits_df['seed_contig_name'].isin(coverage_df['seed_contig_name'].unique())]
        ##################################################

        # Add seed gene callers IDs to the table. (Both seed "contig" names and gene callers IDs are
        # unique.)
        with self.trnaseq_contigs_db_info.load_db() as trnaseq_contigs_db:
            seed_id_df = trnaseq_contigs_db.get_table_as_dataframe(
                'genes_in_contigs', columns_of_interest=['gene_callers_id', 'contig'])
        seed_id_df = seed_id_df.rename(
            {'gene_callers_id': 'seed_gene_callers_id', 'contig': 'seed_contig_name'}, axis=1)
        hits_df = hits_df.merge(seed_id_df, how='left', on='seed_contig_name')

        # Polish the columns. Order them how they will appear in the hits table in the database.
        hits_df['seed_alignment_start'] = hits_df['qstart'] - 1
        hits_df['gene_alignment_start'] = hits_df['sstart'] - 1
        hits_df = hits_df.drop(['qstart', 'qlen', 'sseqid', 'sstart'], axis=1)
        hits_df = hits_df.rename(
            {'gene_callers_id': 'gene_gene_callers_id', 'send': 'gene_alignment_stop'}, axis=1)
        hits_df = hits_df[['seed_gene_callers_id', # seed info
                           'seed_contig_name',
                           'seed_permutation',
                           'contigs_db_project_name', # gene and genome info
                           'contigs_db_hash',
                           'gene_contig_name',
                           'profile_db_sample_id',
                           'collection_name',
                           'bin_id',
                           'gene_gene_callers_id',
                           'decoded_amino_acid',
                           'anticodon',
                           'gene_start_in_contig',
                           'gene_stop_in_contig',
                           'trnascan_score',
                           'gene_sequence',
                           'mismatch', # hit info
                           'bitscore',
                           'seed_alignment_start',
                           'gene_alignment_start',
                           'gene_alignment_stop']]

        return hits_df


    def find_unmodified_nucleotides(self, hits_df):
        """
        Find the unmodified nucleotides at predicted modification positions in tRNA-seq seeds
        using matching tRNA gene sequences.

        Parameters
        ==========
        hits_df : pandas.core.frame.DataFrame
            Each row contains a selected hit between a seed and tRNA gene.

        Returns
        =======
        unmodified_nt_df : pandas.core.frame.DataFrame
            Each row contains a modification for which the underlying nucleotide could be resolved.
        """
        # Load modification information for seeds associated with genes.
        modifications_df = pd.read_csv(self.modifications_txt_path, sep='\t', header=0,
                                       usecols=['contig_name', 'seed_position'])
        modifications_df = modifications_df.rename({'contig_name': 'seed_contig_name'}, axis=1)
        modifications_df = modifications_df.drop_duplicates()
        modifications_df = modifications_df[
            modifications_df['seed_contig_name'].isin(hits_df['seed_contig_name'].unique())]

        # If a seed matches multiple genes and the nucleotides at a predicted modification position
        # differ between the genes, then it is likely that the variation is genetic rather than
        # caused by a modification.
        modification_candidates_df = modifications_df.merge(
            hits_df[['seed_contig_name', 'seed_alignment_start', 'gene_sequence']],
            how='left', on='seed_contig_name')
        modification_keys = []
        unmodified_nts = []
        variant_keys = []
        for group_key, modification_candidate_df in modification_candidates_df.groupby(
            ['seed_contig_name', 'seed_position'], as_index=False):
            unmodified_nt = ''
            for seed_position, seed_alignment_start, gene_sequence in zip(
                modification_candidate_df['seed_position'],
                modification_candidate_df['seed_alignment_start'],
                modification_candidate_df['gene_sequence']):
                gene_nt = gene_sequence[int(seed_position - seed_alignment_start)]
                if unmodified_nt:
                    if gene_nt != unmodified_nt:
                        variant_keys.append(group_key)
                        break
                else:
                    unmodified_nt = gene_nt
            else:
                modification_keys.append(group_key)
                unmodified_nts.append(unmodified_nt)
        modifications_df = modifications_df.set_index(['seed_contig_name', 'seed_position'])
        modifications_df = modifications_df.loc[modification_keys]
        unmodified_nt_df = modifications_df.merge(
            pd.DataFrame(
                [modification_key + (unmodified_nt, )
                 for modification_key, unmodified_nt in zip(modification_keys, unmodified_nts)],
                columns=['seed_contig_name', 'seed_position', 'unmodified_nt']),
            on=['seed_contig_name', 'seed_position'])
        unmodified_nt_df = unmodified_nt_df.reset_index(drop=True)

        return unmodified_nt_df


    def update_trnaseq_contigs_database(self, hits_df, unmodified_nt_df):
        """
        Add information on tRNA gene associations to the tRNA-seq contigs database.

        Parameters
        ==========
        hits_df : pandas.core.frame.DataFrame
            Each row contains a selected hit between a seed and tRNA gene.

        unmodified_nt_df : pandas.core.frame.DataFrame
            Each row contains a modification for which the underlying nucleotide could be resolved.
        """
        trnaseq_contigs_db = self.trnaseq_contigs_db_info.load_db()

        # Either clear the table of seed/gene matches with `self.remove_previous_matches` or append
        # to the table.
        row_count = trnaseq_contigs_db.get_row_counts_from_table(tables.trna_gene_hits_table_name)
        print()
        if self.remove_previous_matches:
            trnaseq_contigs_db._exec(f'''DELETE FROM {tables.trna_gene_hits_table_name}''')
            self.run.info_single(
                f"{pp(row_count)} seed/gene matches dropped from the tRNA-seq contigs database",
                cut_after=0)
        else:
            if row_count:
                self.run.info_single(
                    f"Appending to the {pp(row_count)} seed/gene matches previously stored in the "
                    "tRNA-seq contigs database",
                    cut_after=0)

        # Add the unmodified nucleotides to the rows of the table.
        table_entries = []
        hit_id = 0
        for row in hits_df.itertuples(index=False):
            seed_unmodified_nt_df = unmodified_nt_df[
                unmodified_nt_df['seed_contig_name'] == row.seed_contig_name]

            if len(seed_unmodified_nt_df):
                # Note that modification positions at which the unmodified nucleotide could not be
                # resolved will not be represented in the entry.
                seed_unmodified_nt_series = (seed_unmodified_nt_df['seed_position'].astype(str) +
                                             seed_unmodified_nt_df['unmodified_nt'])
                unmodified_nt_entry = ','.join(seed_unmodified_nt_series.tolist())
            else:
                unmodified_nt_entry = ''

            table_entries.append([hit_id,
                                  row.seed_gene_callers_id, # seed info
                                  row.seed_contig_name,
                                  row.contigs_db_project_name, # gene and genome info
                                  row.contigs_db_hash,
                                  row.gene_contig_name,
                                  row.profile_db_sample_id,
                                  row.collection_name,
                                  row.bin_id,
                                  row.gene_gene_callers_id,
                                  row.decoded_amino_acid,
                                  row.anticodon,
                                  row.gene_start_in_contig,
                                  row.gene_stop_in_contig,
                                  row.trnascan_score,
                                  row.gene_sequence,
                                  row.mismatch, # hit info
                                  row.bitscore,
                                  row.seed_alignment_start,
                                  row.gene_alignment_start,
                                  row.gene_alignment_stop,
                                  unmodified_nt_entry])
            hit_id += 1

        if (len(trnaseq_contigs_db.get_table_column_types(tables.trna_gene_hits_table_name)) !=
            len(tables.trna_gene_hits_table_types)):
            # `sanity_check` would have raised an error if the structure of the table is incorrect
            # and `remove_previous_matches` is False.
            if self.remove_previous_matches:
                trnaseq_contigs_db.drop_table(tables.trna_gene_hits_table_name)
                trnaseq_contigs_db.create_table(tables.trna_gene_hits_table_name,
                                                tables.trna_gene_hits_table_structure,
                                                tables.trna_gene_hits_table_types)
                self.run.warning(
                    "The existing table of seed/gene matches in the tRNA-seq contigs database, "
                    f"'{self.trnaseq_contigs_db_path}', did not have the proper number of columns, "
                    "probably because it was created with an obsolete version of the program, but "
                    "anvi'o was instructed to clear the table with `remove_previous_matches`, so "
                    "the table has instead been deleted entirely and replaced by one with the "
                    "proper column structure.")

        trnaseq_contigs_db._exec_many(
            f'''INSERT INTO {tables.trna_gene_hits_table_name} VALUES '''
            f'''({",".join(["?"] * len(tables.trna_gene_hits_table_structure))})''', table_entries)

        trnaseq_contigs_db.disconnect()

        self.run.info_single(
            f"{hits_df['seed_gene_callers_id'].nunique()} tRNA-seq seeds are found to match "
            f"{hits_df['gene_gene_callers_id'].nunique()} tRNA genes",
            cut_after=0, nl_after=1)


    @staticmethod
    def get_integrated_genomes(trnaseq_contigs_db_path):
        """
        Get the project names and hashes of (meta)genomic contigs databases and any bins in which
        tRNA genes linked to the 'trnaseq'-variant contigs database were found.

        Parameters
        ==========
        trnaseq_contigs_db_path : str
            Path to 'trnaseq'-variant contigs database.

        Returns
        =======
        integrated_genome_dict : dict
            Nested dictionary with levels for contigs database, profile database, collection, and
            bins.
        """
        trnaseq_contigs_db_info = DBInfo(trnaseq_contigs_db_path, expecting='contigs')

        if trnaseq_contigs_db_info.variant != 'trnaseq':
            raise ConfigError(
                f"The database at '{trnaseq_contigs_db_path}' was a "
                f"'{trnaseq_contigs_db_info.variant}' variant, not the required 'trnaseq' variant.")

        with trnaseq_contigs_db_info.load_db() as trnaseq_contigs_db:
            trna_gene_hits_df = trnaseq_contigs_db.get_table_as_dataframe(
                tables.trna_gene_hits_table_name,
                columns_of_interest=['gene_contigs_db_project_name',
                                     'gene_contigs_db_hash',
                                     'profile_db_sample_id',
                                     'collection_name',
                                     'bin_id'])
        integrated_genome_dict = {}
        for row in trna_gene_hits_df.itertuples(index=False):
            if not row.gene_contigs_db_project_name or not row.gene_contigs_db_hash:
                raise ConfigError(
                    "For some reason a row of the tRNA gene hits table in the tRNA-seq contigs "
                    f"database at '{trnaseq_contigs_db_path}' does not have a proper (meta)genomic "
                    "contigs database identifier indicating the source of the gene. The contigs "
                    "database should be identified by both a project name and a hash. Here is all "
                    "of the information from the erroneous row. Contigs database project name: "
                    f"{row.gene_contigs_db_project_name} ; Contigs database hash: "
                    f"{row.gene_contigs_db_hash} ; Profile database sample ID: "
                    f"{row.profile_db_sample_id} ; Collection name: {row.collection_name} ; Bin "
                    f"ID: {row.bin_id}")

            contigs_db_key = (row.gene_contigs_db_project_name, row.gene_contigs_db_hash)
            try:
                profile_db_dict = integrated_genome_dict[contigs_db_key]
            except KeyError:
                integrated_genome_dict[contigs_db_key] = profile_db_dict = {}

            if not row.profile_db_sample_id and not row.collection_name and not row.bin_id:
                continue
            elif row.profile_db_sample_id and row.collection_name and row.bin_id:
                pass
            else:
                raise ConfigError(
                    "For some reason a row of tRNA gene hits only contains some but not all of the "
                    "information needed to identify a bin. A profile database sample ID, "
                    "collection name, and bin ID should all be provided. Here is all of the "
                    "information from the erroneous row. Contigs database project name: "
                    f"{row.gene_contigs_db_project_name} ; Contigs database hash: "
                    f"{row.gene_contigs_db_hash} ; Profile database sample ID: "
                    f"{row.profile_db_sample_id} ; Collection name: {row.collection_name} ; Bin "
                    f"ID: {row.bin_id}")

            try:
                collections_dict = profile_db_dict[row.profile_db_sample_id]
            except KeyError:
                profile_db_dict[row.profile_db_sample_id] = collections_dict = {}

            try:
                bin_ids = collections_dict[row.collection_name]
            except KeyError:
                collections_dict[row.collection_name] = bin_ids = set()
            bin_ids.add(row.bin_id)

        return integrated_genome_dict


class Affinitizer:
    """Using the `go` method, relates changes in tRNA-seq seed abundances to the codon usage of gene
    functions."""

    default_min_coverage_ratio = 5
    default_function_sources = ['KEGG_BRITE']
    default_min_coverage = 10
    default_min_isoacceptors = 4

    recognized_anticodon_wobble_modifications = ['I', 'L']
    # Default decoding weights are from Table 2 of dos Reis, Savva, and Wernisch (2004), used in
    # their tRNA Adaptation Index (tAI) metric: https://doi.org/10.1093/nar/gkh834
    default_decoding_weights_df = pd.DataFrame([
        [1, 1, 1, 0],
        [1, 1, 0, 1],
        [1, 0, 1, 0.41],
        [0, 1, 0.68, 1],
        [0.9999, 0.28, 1, 0],
        [0.89, 1, 1, 1]],
        index=['A', 'C', 'G', 'T', 'I', 'L'],
        columns=['A', 'C', 'G', 'T'])

    builtin_function_blacklists = {
        # Exclude terms associated with eukaryotes, archaea, and photoautotrophs. (The name of the
        # blacklist is misleading in that no terms have yet been included to exclude other types of
        # bacterial autotrophs.)
        'bacterial_heterotrophs': [
            '[eE]ukaryot',
            '[mM]itochondria',
            '[eE]xosom',
            '[pP]roteasom',
            'CD molecule',
            '[eE]ndocytosis',
            '[eE]xocytosis'
            '[aA]rchaea',
            '[pP]hotosynthe'
        ]
    }

    affinity_normalization_methods = ['min_max', 'min_max_mean', 'magnitude_min_max']


    def __init__(self, args={}, r=run, rq=run_quiet, p=progress, do_sanity_check=True):
        self.args = args
        A = lambda x: args.__dict__[x] if x in args.__dict__ else None

        self.trnaseq_contigs_db_path = A('trnaseq_contigs_db')
        self.seeds_specific_txt_path = A('seeds_specific_txt')

        self.reference_sample_name = A('reference_sample')
        self.nonreference_sample_names = A('nonreference_samples')

        self.genomic_contigs_db_path = A('contigs_db')
        self.genomic_profile_db_path = A('profile_db')
        self.collection_name = A('collection_name')
        self.bin_id = A('bin_id')
        self.bin_ids_path = A('bin_ids_file')

        self.internal_genomes_path = A('internal_genomes')
        self.external_genomes_path = A('external_genomes')

        self.seed_assignment = A('seed_assignment')
        self.min_coverage_ratio = A('min_coverage_ratio')

        self.function_sources = A('function_sources')
        self.function_accessions = A('function_accessions')
        if self.function_accessions is None:
            self.function_accessions = []
        self.function_accessions_dict = A('function_accessions_dict')
        if self.function_accessions_dict is None:
            self.function_accessions_dict = {}
        self.function_names = A('function_names')
        if self.function_names is None:
            self.function_names = []
        self.function_names_dict = A('function_names_dict')
        if self.function_names_dict is None:
            self.function_names_dict = {}
        self.select_functions_txt = A('select_functions_txt')
        self.lax_function_sources = A('lax_function_sources')
        if self.lax_function_sources is None:
            self.lax_function_sources = False
        self.function_blacklist_txt = A('function_blacklist_txt')
        self.gene_affinity = A('gene_affinity')
        if self.gene_affinity is None:
            self.gene_affinity = False
        self.gene_caller_ids = A('gene_caller_ids')
        if self.function_sources is None:
            if self.gene_affinity or self.function_accessions_dict or self.function_names_dict:
                self.function_sources = []
            else:
                self.function_sources = self.default_function_sources

        self.min_coverage = A('min_coverage')
        if self.min_coverage is None:
            self.min_coverage = self.default_min_coverage
        self.min_isoacceptors = A('min_isoacceptors')
        if self.min_isoacceptors is None:
            self.min_isoacceptors = self.default_min_isoacceptors
        self.exclude_unmodified_anticodons = A('exclude_unmodified_anticodons')
        self.exclude_modified_anticodons = A('exclude_modified_anticodons')
        self.decoding_weights_df = A('decoding_weights')
        if self.decoding_weights_df is None:
            self.decoding_weights_df = self.default_decoding_weights_df

        self.min_analyzed_codons = A('min_analyzed_codons')
        if self.min_analyzed_codons is None:
            self.min_analyzed_codons = 0
        self.function_min_total_codons = A('function_min_total_codons')
        if self.function_min_total_codons is None:
            self.function_min_total_codons = 0
        self.gene_min_total_codons = A('gene_min_total_codons')
        if self.gene_min_total_codons is None:
            self.gene_min_total_codons = 0
        self.exclude_codons = A('exclude_codons')
        self.exclude_amino_acids = A('exclude_amino_acids')

        self.rarefaction_limit = A('rarefaction_limit')
        if self.rarefaction_limit is None:
            self.rarefaction_limit = 0

        self.run = r
        self.run_quiet = rq
        self.progress = p

        self.trnaseq_contigs_db_info = DBInfo(self.trnaseq_contigs_db_path, expecting='contigs')

        # Store information on accessing (meta)genomes.
        self.genome_info_dict = {}

        # If `contigs_db` was provided, then `internal_genomes` and `external_genomes` should not
        # also have been provided, which is checked later in `sanity_check`.
        if self.genomic_contigs_db_path:
            contigs_db_info = DBInfo(self.genomic_contigs_db_path, expecting='contigs')
        else:
            contigs_db_info = None
        if self.genomic_profile_db_path:
            profile_db_info = DBInfo(self.genomic_profile_db_path, expecting='profile')
            profile_db_sample_id = profile_db_info.get_self_table()['sample_id']
        else:
            profile_db_info = None
            profile_db_sample_id = None
        if self.bin_ids_path:
            filesnpaths.is_file_plain_text(self.bin_ids_path)
            bin_ids = []
            with open(self.bin_ids_path) as bin_ids_file:
                for line in bin_ids_file:
                    bin_ids.append(line.rstrip())
        else:
            bin_ids = None

        if self.bin_id:
            collect = ccollections.Collections()
            collect.populate_collections_dict(self.genomic_profile_db_path)
            collect.is_bin_in_collection(bin_id)
            self.genome_info_dict[bin_id] = genome_info = {}
            genome_info['contigs_db_info'] = contigs_db_info
            genome_info['profile_db_info'] = profile_db_info
            genome_info['profile_db_sample_id'] = profile_db_sample_id
            genome_info['collection_name'] = self.collection_name
            genome_info['bin_id'] = self.bin_id
        elif self.bin_ids_path:
            collect = ccollections.Collections()
            collect.populate_collections_dict(self.genomic_profile_db_path)
            for bin_id in bin_ids:
                collect.is_bin_in_collection(bin_id)
                self.genome_info_dict[bin_id] = genome_info = {}
                genome_info['contigs_db_info'] = contigs_db_info
                genome_info['profile_db_info'] = profile_db_info
                genome_info['profile_db_sample_id'] = profile_db_sample_id
                genome_info['collection_name'] = self.collection_name
                genome_info['bin_id'] = bin_id
        elif self.collection_name:
            # There is a collection of internal genomes.
            collect = ccollections.Collections()
            collect.populate_collections_dict(self.genomic_profile_db_path)
            for bin_id in collect.get_bins_info_dict():
                self.genome_info_dict[bin_id] = genome_info = {}
                genome_info['contigs_db_info'] = contigs_db_info
                genome_info['profile_db_info'] = profile_db_info
                genome_info['profile_db_sample_id'] = profile_db_sample_id
                genome_info['collection_name'] = self.collection_name
                genome_info['bin_id'] = bin_id
        elif self.genomic_contigs_db_path:
            # The contigs database represents a genome, like an external genome.
            self.genome_info_dict[contigs_db_info.project_name] = genome_info = {}
            genome_info['contigs_db_info'] = contigs_db_info
            genome_info['profile_db_info'] = None
            genome_info['profile_db_sample_id'] = None
            genome_info['collection_name'] = None
            genome_info['bin_id'] = None

        nonunique_genome_names = []
        if self.internal_genomes_path or self.external_genomes_path:
            descriptions = GenomeDescriptions(args, run=self.run_quiet, progress=self.progress)
            descriptions.load_genomes_descriptions(init=False)

            for genome_name, genome_dict in descriptions.internal_genomes_dict.items():
                if genome_name in self.genome_info_dict:
                    nonunique_genome_names.append(genome_name)

                self.genome_info_dict[genome_name] = genome_info = {}
                genome_info['contigs_db_info'] = DBInfo(
                    genome_dict['contigs_db_path'], expecting='contigs')
                if genome_dict['profile_db_path']:
                    genome_info['profile_db_info'] = DBInfo(
                        genome_dict['profile_db_path'], expecting='profile')
                    genome_info['profile_db_sample_id'] = genome_info[
                        'profile_db_info'].get_self_table()['sample_id']
                else:
                    genome_info['profile_db_info'] = None
                    genome_info['profile_db_sample_id'] = None
                genome_info['collection_name'] = genome_dict['collection_id']
                genome_info['bin_id'] = genome_dict['bin_id']

            for genome_name, genome_dict in descriptions.external_genomes_dict.items():
                if genome_name in self.genome_info_dict:
                    nonunique_genome_names.append(genome_name)

                self.genome_info_dict[genome_name] = genome_info = {}
                genome_info['contigs_db_info'] = DBInfo(
                    genome_dict['contigs_db_path'], expecting='contigs')
                genome_info['profile_db_info'] = None
                genome_info['profile_db_sample_id'] = None
                genome_info['collection_name'] = None
                genome_info['bin_id'] = None

        if nonunique_genome_names:
            raise ConfigError(
                "Names must be unique to internal and external genomes, but the following were "
                f"not: {', '.join(nonunique_genome_names)}")

        if do_sanity_check:
            # No object attributes are assigned or modified in `sanity_check`.
            self.sanity_check()

        # Find the names of the samples to analyze in addition to the reference sample. By default,
        # every other sample will be a non-reference sample, but a subset of available samples can
        # be provided instead.
        if self.nonreference_sample_names is None:
            self.nonreference_sample_names = pd.read_csv(
                self.seeds_specific_txt_path,
                sep='\t',
                header=0,
                skiprows=[1, 2],
                usecols=['sample_name'])['sample_name'].unique().tolist()
            self.nonreference_sample_names.remove(self.reference_sample_name)
        self.sample_names = [self.reference_sample_name] + self.nonreference_sample_names

        self.function_sources += \
            list(set(list(self.function_accessions_dict) + list(self.function_names_dict)))

        if self.function_accessions:
            self.function_accessions_dict[self.function_sources[0]] = self.function_accessions
        if self.function_names:
            self.function_names_dict[self.function_sources[0]] = self.function_names

        if self.select_functions_txt:
            select_functions_df = pd.read_csv(
                self.select_functions_txt, sep='\t', header=None,
                names=['function_source', 'function_accession', 'function_name'])
            select_functions_df = select_functions_df.fillna('')
            for row in select_functions_df.itertuples():
                if row.accession:
                    try:
                        self.function_accessions_dict[row.source].append(row.accession)
                    except KeyError:
                        self.function_accessions_dict[row.source] = [row.accession]
                elif row.name:
                    try:
                        self.function_names_dict[row.source].append(row.name)
                    except KeyError:
                        self.function_names_dict[row.source] = [row.name]
            for function_source in self.function_accessions_dict:
                if function_source not in self.function_sources:
                    self.function_sources.append(function_source)
            for function_source in self.function_names_dict:
                if function_source not in self.function_sources:
                    self.function_sources.append(function_source)

        self.function_blacklist_patterns = []
        if self.function_blacklist_txt:
            for line in open(self.function_blacklist_txt):
                self.function_blacklist_patterns.append(line.rstrip())

        # Find which codons are decoded by which anticodons. A decoding efficiency of 1 in
        # `decoding_weights_df` means that the anticodon does not decode the codon. The following
        # dictionary maps each anticodon to a list of decoded codons. Note that anticodons for tRNAs
        # that do not exist are included in the dictionary, such as TTA decoding TAA (STP): only
        # entries with observed tRNAs in the tRNA-seq data and corresponding codons in genes
        # contribute to the result.
        self.nucleotide_decoding_dict = {}
        for anticodon in constants.anticodon_to_AA:
            anticodon_wobble_nucleotide = anticodon[0]
            # Two anticodon wobble modifications are treated separately, as in tAI: ANN -> INN and
            # CAT -> LAT. These have separate entries in `decoding_weights_df`.
            derived_anticodons = [anticodon]
            if anticodon_wobble_nucleotide == 'A':
                derived_anticodons.append('I' + anticodon[1: ])
            elif anticodon == 'CAT':
                derived_anticodons.append('LAT')

            for derived_anticodon in derived_anticodons:
                derived_anticodon_wobble_nucleotide = derived_anticodon[0]
                decoding_weights_series = self.decoding_weights_df.loc[
                    derived_anticodon_wobble_nucleotide]
                decoded_codons = []
                for codon in constants.codons:
                    codon_RC = constants.codon_to_codon_RC[codon]
                    if anticodon[1: ] != codon_RC[1: ]:
                        # The anticodon must form Watson-Crick base pairs with the codon at the
                        # non-wobble positions.
                        continue

                    codon_wobble_nucleotide = codon[2]
                    decoding_weight = decoding_weights_series.loc[codon_wobble_nucleotide]
                    if decoding_weight < 1:
                        decoded_codons.append(codon)
                self.nucleotide_decoding_dict[derived_anticodon] = decoded_codons


    def sanity_check(self):
        """Check the feasibility of args from initialization."""
        if self.trnaseq_contigs_db_info.variant != 'trnaseq':
            raise ConfigError(
                f"The database at '{self.trnaseq_contigs_db_path}' was a "
                f"'{self.trnaseq_contigs_db_info.variant}' variant, not the required 'trnaseq' "
                "variant.")

        with self.trnaseq_contigs_db_info.load_db() as trnaseq_contigs_db:
            if trnaseq_contigs_db.get_row_counts_from_table(tables.trna_gene_hits_table_name) == 0:
                raise ConfigError(
                    "No tRNA seeds in the tRNA-seq contigs database, "
                    f"'{self.trnaseq_contigs_db_path}', are linked to tRNA genes, a task performed "
                    "by `anvi-integrate-trnaseq`.")

        ##################################################
        # Check for tRNA-seq reference and non-reference samples.
        filesnpaths.is_file_exists(self.seeds_specific_txt_path)

        if self.reference_sample_name is None:
            raise ConfigError("A reference sample must be provided.")

        available_sample_names = pd.read_csv(
            self.seeds_specific_txt_path,
            sep='\t',
            header=0,
            skiprows=[1, 2],
            usecols=['sample_name'])['sample_name'].unique().tolist()

        if self.reference_sample_name not in available_sample_names:
            raise ConfigError(
                f"The provided reference sample name, '{self.reference_sample_name}', is not found "
                f"in `seeds_specific_txt`, '{self.seeds_specific_txt_path}', the table of specific "
                "coverages of tRNA-seq seeds in different samples. Here are the samples provided "
                f"in that table: {', '.join(available_sample_names)}")

        if self.nonreference_sample_names:
            nonreference_sample_names = self.nonreference_sample_names
            bad_sample_names = set(self.nonreference_sample_names).difference(
                set(available_sample_names))
            if bad_sample_names:
                raise ConfigError(
                    "The following provided non-reference sample names were not found in "
                    f"`seeds_specific_txt`, '{self.seeds_specific_txt_path}', the table of "
                    "specific coverages of tRNA-seq seeds in different samples: "
                    f"{', '.join(bad_sample_names)}. Here are the samples provided in that table: "
                    f"{', '.join(available_sample_names)}")

            if self.reference_sample_name in self.nonreference_sample_names:
                raise ConfigError(
                    f"Please do not include the reference sample, '{self.reference_sample_name}' "
                    f"in the subset of sample names: {', '.join(nonreference_sample_names)}. Sorry "
                    "for the sclerotic idiocy.")
        else:
            nonreference_sample_names = available_sample_names
            nonreference_sample_names.remove(self.reference_sample_name)
        if len(nonreference_sample_names) == 0:
            raise ConfigError(
                "There must be one or more samples beside the reference sample in "
                f"`seeds-specific-txt`, '{self.seeds_specific_txt}'. Only the reference sample, "
                f"'{self.reference_sample_name}', was found.")
        ##################################################

        # Do basic checks of the combinations of (meta)genomic input arguments.
        if (self.genomic_contigs_db_path and
            (self.internal_genomes_path or self.external_genomes_path)):
            raise ConfigError(
                "`internal_genomes` and `external_genomes` cannot be used with `contigs_db`.")

        if ((self.genomic_profile_db_path or self.collection_name) and
            not (
                self.genomic_contigs_db_path and
                self.genomic_profile_db_path and
                self.collection_name)):
            raise ConfigError(
                "A collection must be provided using `contigs_db`, `profile_db`, and "
                "`collection_name`.")

        if (self.bin_id and
            not (
                self.genomic_contigs_db_path and
                self.genomic_profile_db_path and
                self.collection_name)):
            raise ConfigError(
                "A specific bin provided with `bin_id` also requires `contigs_db`, `profile_db`, "
                "and `collection_name`.")

        ##################################################
        # Check that the (meta)genomic inputs correspond to tRNA genes linked to seeds in the
        # tRNA-seq contigs database. At least one of the (meta)genomic inputs must contain linked
        # genes: not all need contain linked genes, as the user may have integrated a number of
        # internal and external genomes, not all of which yielded seed/gene matches. Report
        # (meta)genomes containing and lacking linked genes.
        with self.trnaseq_contigs_db_info.load_db() as trnaseq_contigs_db:
            genome_id_df = trnaseq_contigs_db.get_table_as_dataframe(
                tables.trna_gene_hits_table_name,
                columns_of_interest=[
                    'gene_contigs_db_hash',
                    'profile_db_sample_id',
                    'collection_name',
                    'bin_id'])
            genome_id_df = genome_id_df.drop_duplicates()
            genome_id_df = genome_id_df.fillna('')
            genome_ids = []
            for row in genome_id_df.itertuples(index=False):
                genome_ids.append((
                    row.gene_contigs_db_hash,
                    row.profile_db_sample_id,
                    row.collection_name,
                    row.bin_id))

        linked_bin_info_dict = {}
        unlinked_bin_info_dict = {}
        linked_nonbin_info_dict = {}
        unlinked_nonbin_info_dict = {}
        for genome_name, genome_info in self.genome_info_dict.items():
            contigs_db_hash = genome_info['contigs_db_info'].hash
            profile_db_sample_id = genome_info['profile_db_sample_id']
            collection_name = genome_info['collection_name']
            bin_id = genome_info['bin_id']
            if bin_id:
                # Considering a bin.
                if (contigs_db_hash,
                    profile_db_sample_id,
                    collection_name,
                    bin_id) in genome_ids:
                    linked_bin_info_dict[genome_name] = genome_info
                else:
                    unlinked_bin_info_dict[genome_name] = genome_info
            else:
                # Considering an external genome or metagenome.
                if (contigs_db_hash, '', '', '') in genome_ids:
                    linked_nonbin_info_dict[genome_name] = genome_info
                else:
                    unlinked_nonbin_info_dict[genome_name] = genome_info

        if not linked_bin_info_dict and not linked_nonbin_info_dict:
            if len(self.genome_info_dict) > 1:
                first_variable_message_string = "None of the provided (meta)genomes contained"
                second_variable_message_string = "(meta)genomes were"
            else:
                first_variable_message_string = "The provided (meta)genome did not contain"
                second_variable_message_string = "(meta)genome was"
            raise ConfigError(
                f"{first_variable_message_string} tRNA genes that matched seeds in the tRNA-seq "
                f"contigs database at '{self.trnaseq_contigs_db_path}'. Either the "
                f"{second_variable_message_string} not run with `anvi-integrate-trnaseq` or did "
                "not contain any tRNA genes matching any tRNA-seq seeds.")

        linked_genome_names = []
        unlinked_genome_names = []
        for genome_name in self.genome_info_dict:
            if genome_name in linked_bin_info_dict or genome_name in linked_nonbin_info_dict:
                linked_genome_names.append(genome_name)
            else:
                unlinked_genome_names.append(genome_name)
        self.run.info_single(
            "(Meta)genome(s) containing tRNA genes linked to seeds: "
            f"{', '.join(linked_genome_names)}")
        if unlinked_genome_names:
            self.run.info_single(
                "(Meta)genomes lacking tRNA genes linked to seeds: "
                f"{', '.join(unlinked_genome_names)}")
        ##################################################

        if self.min_coverage_ratio <= 1:
            raise ConfigError(
                "`min_coverage_ratio` must be >1, as in the context of the `seed_assignment` "
                "option, `ambiguous_choose`, it only makes sense for an ambiguous seed to be "
                "assigned to the genome with the maximum abundance of unambiguous seeds.")

        if len(self.function_sources) > 1 and (self.function_accessions or self.function_names):
            raise ConfigError(
                "`function_accessions` and `function_names` require a single function source, but "
                f"multiple `function_sources` were provided: {', '.join(self.function_sources)}.")

        if self.function_sources and (self.function_accessions_dict or self.function_names_dict):
            raise ConfigError(
                "`function_sources` cannot be provided alongside `function_accesions_dict` and "
                "`function_names_dict`.")

        if self.select_functions_txt:
            filesnpaths.is_file_tab_delimited(self.select_functions_txt)

        if self.select_functions_txt and (
            self.function_sources or self.function_accessions or self.function_names):
            raise ConfigError(
                "`select_functions_txt` cannot be provided alongside `function_sources`, "
                "`function_accessions`, and `function_names`.")

        if self.function_blacklist_txt:
            filesnpaths.is_file_tab_delimited(self.function_blacklist_txt)

        if self.gene_affinity and (
            self.function_sources or self.function_accessions or self.function_accessions_dict or
            self.function_names or self.function_names_dict or self.select_functions_txt or
            self.function_blacklist_txt):
            raise ConfigError(
                "`gene_affinity`, used to compute affinity for genes rather than functions, cannot "
                "be used alongside function options.")

        # Report input (meta)genomes that did not have requested function sources run on them. By
        # default, with `lax_function_sources` being True, requested sources must have been run on
        # every (meta)genome.
        if self.function_sources:
            function_sources = self.function_sources
        elif self.function_accessions_dict or self.function_names_dict or self.select_functions_txt:
            function_sources = list(self.function_accessions_dict)
            function_sources += list(self.function_names_dict)
            if self.select_functions_txt:
                select_functions_df = pd.read_csv(
                    self.select_functions_txt, sep='\t', header=None,
                    names=['function_source', 'function_accession', 'function_name'])
                function_sources += select_functions_df['function_source'].unique().tolist()
        present_function_source_genome_dict = {
            function_source: [] for function_source in function_sources}
        missing_function_source_genome_dict = {
            function_source: [] for function_source in function_sources}
        for genome_name, genome_info in self.genome_info_dict.items():
            genome_function_sources = genome_info['contigs_db_info'].get_self_table()[
                'gene_function_sources'].split(',')
            for function_source in function_sources:
                if function_source in genome_function_sources:
                    present_function_source_genome_dict[function_source].append(genome_name)
                else:
                    missing_function_source_genome_dict[function_source].append(genome_name)
        message = ""
        for function_source, genome_names in missing_function_source_genome_dict.items():
            if len(genome_names):
                message += f"{function_source}: {', '.join(genome_names)} ; "
        if message:
            message = "(Meta)genomes lack requested function source(s): " + message
            message = message[: -3]
            if self.lax_function_sources:
                self.run.info_single(message)
            else:
                raise ConfigError(message)

        if self.gene_caller_ids and not self.gene_affinity:
            raise ConfigError("`gene_caller_ids` requires the `gene_affinity` option.")

        if self.min_coverage < 1:
            raise ConfigError(
                "The minimum coverage for a tRNA isoacceptor to be detected must be a positive "
                f"integer, not the provided `min_coverage`, {self.min_coverage}.")

        if self.min_isoacceptors < 1:
            raise ConfigError(
                "The minimum number of tRNA isoacceptors for translational affinity to be "
                "calculated must be a positive integer, not the provided `min_isoacceptors` of "
                f"{self.min_isoacceptors}.")

        if self.function_min_total_codons and self.gene_affinity:
            raise ConfigError(
                "`function_min_total_codons` cannot be provided in conjunction with "
                "`gene_affinity`. Filtering functions by codon count does not affect the selection "
                "of genes for gene affinity calculations.")
        if self.gene_min_total_codons and not self.gene_affinity:
            raise ConfigError(
                "`gene_min_total_codons` cannot be provided without `gene_affinity`. Genes cannot "
                "be filtered by codon count prior to calculation of functional affinity in this "
                "implementation.")

        unrecognized_unmodified_anticodons = set(self.exclude_unmodified_anticodons).difference(
            constants.anticodon_to_AA) if self.exclude_unmodified_anticodons else []
        if unrecognized_unmodified_anticodons:
            raise ConfigError(
                "The following unmodified anticodons are unrecognized and cannot be excluded from "
                f"affinity calculations: {', '.join(unrecognized_unmodified_anticodons)}")

        unrecognized_modified_anticodons = []
        if self.exclude_modified_anticodons:
            for modified_anticodon in self.exclude_modified_anticodons:
                if modified_anticodon[0] not in self.recognized_anticodon_wobble_modifications:
                    unrecognized_modified_anticodons.append(modified_anticodon)
                    continue
                if (modified_anticodon[1] not in constants.unambiguous_nucleotides or
                    modified_anticodon[2] not in constants.unambiguous_nucleotides):
                    unrecognized_modified_anticodons.append(modified_anticodon)
        if unrecognized_modified_anticodons:
            raise ConfigError(
                "The following modified anticodons are unrecognized and cannot be excluded from "
                f"affinity calculations: {', '.join(unrecognized_modified_anticodons)}")

        unrecognized_codons = set(self.exclude_codons).difference(
            constants.codons) if self.exclude_codons else []
        if unrecognized_codons:
            raise ConfigError(
                "Unrecognized codons were provided to `exclude_codons`: "
                f"{', '.join(unrecognized_codons)}")

        unrecognized_amino_acids = set(self.exclude_amino_acids).difference(
            constants.amino_acids) if self.exclude_amino_acids else []
        if unrecognized_amino_acids:
            raise ConfigError(
                "Unrecognized amino acids were provided to `exclude_amino_acids`: "
                f"{', '.join(unrecognized_amino_acids)}")


    def go(self):
        """Relate changes in tRNA-seq seed abundances to the codon usage of gene functions."""
        isoacceptor_abund_ratios_df = self.get_isoacceptors()
        if len(isoacceptor_abund_ratios_df) == 0:
            run.warning("Affinity could not be calculated given the lack of tRNA-seq data passing "
                        "the filters.")
            return

        weighted_frequencies_df = self.get_weighted_codon_frequencies(isoacceptor_abund_ratios_df)
        if len(weighted_frequencies_df) == 0:
            run.warning(
                "Affinity could not be calculated given the lack of "
                f"{'genes' if self.gene_affinity else 'functions'} passing the codon filters.")
            return

        affinities_df = self.get_affinities(isoacceptor_abund_ratios_df, weighted_frequencies_df)
        if len(affinities_df) == 0:
            run.warning(
                "Affinity could not be calculated despite the presence of the prerequisite "
                f"filtered tRNA-seq and {'gene' if self.gene_affinity else 'function'} codon "
                "frequency data.")

        return affinities_df


    def get_isoacceptors(self):
        """
        Get a table of per-genome isoacceptor non-reference/reference abundance ratios to use in
        affinity calculations. Isoacceptors are groups of tRNA-seq seeds with the same anticodon
        that are assigned to a genome.

        Returns
        =======
        isoacceptor_abund_ratios_df : pandas.core.frame.DataFrame
            Each row contains genome and isoacceptor information identifying non-reference/reference
            abundance ratios. This table is "long," with one column of abundance ratio data.
        """
        # Load data from the tRNA-seq contigs database.
        trnaseq_contigs_db = self.trnaseq_contigs_db_info.load_db()

        trna_gene_hits_df = trnaseq_contigs_db.get_table_as_dataframe(
            tables.trna_gene_hits_table_name,
            columns_of_interest=['seed_gene_callers_id',
                                 'seed_contig_name',
                                 'gene_contigs_db_hash',
                                 'profile_db_sample_id',
                                 'collection_name',
                                 'bin_id',
                                 'decoded_amino_acid',
                                 'anticodon'])

        # Dereplicate duplicate rows representing hits between the same tRNA-seq seed and
        # different tRNA genes with identical sequences. (These rows would be distinguished by
        # the column `gene_gene_callers_id` if that were loaded.)
        trna_gene_hits_df = trna_gene_hits_df.drop_duplicates()

        for column_name in ('profile_db_sample_id', 'collection_name', 'bin_id'):
            trna_gene_hits_df[column_name] = trna_gene_hits_df[column_name].fillna('')

        if self.seed_assignment == 'unambiguous_db':
            # Disregard seeds with ambiguous matches in the database, not just seeds with
            # ambiguous matches among the possible subset of genomes input to this program.
            trna_gene_hits_df = trna_gene_hits_df.groupby('seed_gene_callers_id').filter(
                lambda seed_df: len(seed_df) == 1)

        if self.seed_assignment == 'ambiguous_choose':
            # The identities of the ambiguous seeds are needed later when choosing genomes.
            ambiguous_seed_gene_callers_ids = trna_gene_hits_df.groupby(
                'seed_gene_callers_id').filter(lambda seed_df: len(seed_df) > 1)[
                    'seed_gene_callers_id'].unique()

        # Select seeds matching input genomes. Replace the four columns needed to uniquely identify
        # a genome with a single column of the unique genome name as given in the input.
        select_genome_ids = []
        genome_id_name_dict = {}
        for genome_name, genome_info in self.genome_info_dict.items():
            contigs_db_hash = genome_info['contigs_db_info'].hash
            if genome_info['profile_db_info'] is None:
                profile_db_sample_id = ''
            else:
                profile_db_sample_id = genome_info['profile_db_sample_id']
            if genome_info['collection_name'] is None:
                collection_name = ''
            else:
                collection_name = genome_info['collection_name']
            if genome_info['bin_id'] is None:
                bin_id = ''
            else:
                bin_id = genome_info['bin_id']
            genome_id = (contigs_db_hash, profile_db_sample_id, collection_name, bin_id)
            select_genome_ids.append(genome_id)
            genome_id_name_dict[genome_id] = genome_name
        trna_gene_hits_df = trna_gene_hits_df.set_index(
            ['gene_contigs_db_hash', 'profile_db_sample_id', 'collection_name', 'bin_id'])
        trna_gene_hits_df = trna_gene_hits_df.loc[select_genome_ids]
        trna_gene_hits_df = trna_gene_hits_df.reset_index()
        trna_gene_hits_df['genome_name'] = [
            genome_id_name_dict[genome_id] for genome_id in zip(
                trna_gene_hits_df['gene_contigs_db_hash'],
                trna_gene_hits_df['profile_db_sample_id'],
                trna_gene_hits_df['collection_name'],
                trna_gene_hits_df['bin_id'])]
        trna_gene_hits_df = trna_gene_hits_df.drop(
            ['gene_contigs_db_hash', 'profile_db_sample_id', 'collection_name', 'bin_id'], axis=1)

        if self.seed_assignment == 'unambiguous_genome':
            # Disregard seeds with ambiguous matches among genomes input to this program.
            trna_gene_hits_df = trna_gene_hits_df.groupby('seed_gene_callers_id').filter(
                lambda seed_df: len(seed_df) == 1)

        ##################################################
        # Find the nucleotide at anticodon wobble position 34 in each seed consensus sequence.
        seed_ids_string = ','.join(['"%s"' % gene_callers_id for gene_callers_id in
                                    trna_gene_hits_df['seed_gene_callers_id'].unique()])
        ids_where_clause = f'''gene_callers_id IN ({seed_ids_string})'''
        wobble_position_df = trnaseq_contigs_db.get_table_as_dataframe(
            'trna_feature',
            columns_of_interest=['gene_callers_id', 'anticodon_loop_start'],
            where_clause=ids_where_clause)
        wobble_position_df = wobble_position_df.rename(
            {'gene_callers_id': 'seed_gene_callers_id'}, axis=1)
        wobble_position_df['anticodon_start'] = wobble_position_df['anticodon_loop_start'] + 2
        wobble_position_df = wobble_position_df.drop('anticodon_loop_start', axis=1)
        trna_gene_hits_df = trna_gene_hits_df.merge(
            wobble_position_df, on='seed_gene_callers_id')

        seed_contig_names_string = ','.join(['"%s"' % seed_contig_name for seed_contig_name in
                                                trna_gene_hits_df['seed_contig_name'].unique()])
        contigs_where_clause = f'''contig IN ({seed_contig_names_string})'''
        seed_consensus_sequence_df = trnaseq_contigs_db.get_table_as_dataframe(
            'contig_sequences', where_clause=contigs_where_clause)
        seed_consensus_sequence_df = seed_consensus_sequence_df.rename(
            {'contig': 'seed_contig_name', 'sequence': 'seed_sequence'}, axis=1)
        trna_gene_hits_df = trna_gene_hits_df.merge(
            seed_consensus_sequence_df, on='seed_contig_name')

        anticodon_wobble_nucleotides = []
        for anticodon_start, seed_consensus_sequence in zip(
            trna_gene_hits_df['anticodon_start'], trna_gene_hits_df['seed_sequence']):
            anticodon_wobble_nucleotides.append(seed_consensus_sequence[anticodon_start])
        trna_gene_hits_df['seed_anticodon_wobble_nucleotide'] = anticodon_wobble_nucleotides

        trna_gene_hits_df = trna_gene_hits_df.drop(['anticodon_start', 'seed_sequence'], axis=1)

        # Evaluate the anticodon wobble nucleotide in the seed.
        effective_wobble_nucleotides = []
        for decoded_aa_type, anticodon, seed_wobble_nucleotide in zip(
            trna_gene_hits_df['decoded_amino_acid'],
            trna_gene_hits_df['anticodon'],
            trna_gene_hits_df['seed_anticodon_wobble_nucleotide']):
            if decoded_aa_type == 'Ile2':
                # tRNA-Ile2 has a wobble nucleotide of lysidine in bacteria or agmatidine in
                # archaea, which are given the same decoding weight and both symbolized by 'L'.
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
        trna_gene_hits_df = trna_gene_hits_df.drop('seed_anticodon_wobble_nucleotide', axis=1)
        ##################################################

        # No more data is loaded from the tRNA-seq contigs database.
        trnaseq_contigs_db.disconnect()

        # Dereplicate duplicate rows representing hits between the same seed and genes with
        # different sequences. This should only occur if the seed is a partial read of the tRNA, and
        # the genes differ beyond the 5' end of the seed. However, confirm that the hits yielded the
        # same anticodon wobble nucleotide, in case I'm missing something.
        trna_gene_hits_df = trna_gene_hits_df.drop_duplicates()
        if trna_gene_hits_df.groupby('seed_gene_callers_id').ngroups != trna_gene_hits_df.groupby(
            ['seed_gene_callers_id', 'effective_wobble_nucleotide']).ngroups:
            confusing_df = trna_gene_hits_df.groupby('seed_gene_callers_id').filter(
                lambda seed_df: len(seed_df) > 1)
            print(confusing_df.to_string())
            raise ConfigError(
                "A strange circumstance has occurred where a tRNA-seq seed linked to tRNA genes "
                "with different sequences was found to have different effective wobble "
                "nucleotides. The tabular entries for the seeds in question were printed above "
                "this error message.")

        # Load seed specific coverages in each sample from the tabular file.
        coverage_df = pd.read_csv(
            self.seeds_specific_txt_path,
            sep='\t',
            header=0,
            skiprows=[1, 2],
            usecols=['gene_callers_id',
                     'sample_name',
                     'discriminator_1',
                     'relative_discriminator_coverage'])
        coverage_df = coverage_df.rename({'gene_callers_id': 'seed_gene_callers_id',
                                          'sample_name': 'trnaseq_sample_name'}, axis=1)
        # Select data for the samples of interest.
        coverage_df = coverage_df[coverage_df['trnaseq_sample_name'].isin(self.sample_names)]
        # Select data for the seeds linked to tRNA genes.
        coverage_df = coverage_df[coverage_df['seed_gene_callers_id'].isin(
            trna_gene_hits_df['seed_gene_callers_id'].unique())]
        # Add the tRNA-seq sample coverage data, multiplying the rows per seed by the number of
        # samples in which the seed in measured.
        trna_gene_hits_df = trna_gene_hits_df.merge(
            coverage_df, how='inner', on='seed_gene_callers_id')

        ##################################################
        # At this point, `trna_gene_hits_df` only contains unambiguous seeds if `seed_assignment` is
        # `unambiguous_genome` or `unambiguous_db`, and still contains ambiguous seeds if
        # `ambiguous_all` or `ambiguous_choose`.
        if self.seed_assignment == 'ambiguous_choose':
            # Include ambiguous seeds if they can be assigned to a single genome due to the
            # preponderance of unambiguous seed abundance in that genome compared to the other
            # genomes in which the ambiguous seed is found.
            unambiguous_hits_df = trna_gene_hits_df.groupby(
                'seed_gene_callers_id').filter(lambda seed_df: len(seed_df) == 1)
            unambiguous_coverage_df = unambiguous_hits_df[
                ['genome_name', 'trnaseq_sample_name', 'discriminator_1']].groupby(
                    ['genome_name', 'trnaseq_sample_name'], as_index=False).aggregate('sum')
            # Store unambiguous coverages in a nested dictionary: tRNA-seq sample -> genome ->
            # summed unambiguous coverage. Note that not every genome need be represented in every
            # tRNA-seq sample.
            unambiguous_coverage_dict = {}
            for row in unambiguous_coverage_df.itertuples(index=False):
                try:
                    genome_coverage_dict = unambiguous_coverage_dict[row.trnaseq_sample_name]
                except KeyError:
                    unambiguous_coverage_dict[row.trnaseq_sample_name] = genome_coverage_dict = {}
                genome_coverage_dict[genome_name] = row.discriminator_1

            def choose_genome(ambiguous_seed_df):
                # First, for each genome/tRNA-seq sample in which the ambiguous seed is found
                # ("select" genomes), retrieve the summed unambiguous coverage from
                # `unambiguous_coverage_dict`. Note that the ambiguous seed may be found in genomes
                # for which there is zero unambiguous coverage in the same samples with coverage of
                # the ambiguous seed -- indeed, it is possible for a genome to have zero unambiguous
                # coverage in all samples.
                select_unambiguous_coverage_dict = {}
                for trnaseq_sample_name, sample_df in ambiguous_seed_df.groupby(
                    'trnaseq_sample_name'):
                    try:
                        select_genome_coverage_dict = select_unambiguous_coverage_dict[
                            trnaseq_sample_name]
                    except KeyError:
                        select_unambiguous_coverage_dict[trnaseq_sample_name] = \
                            select_genome_coverage_dict = {}
                    try:
                        genome_coverage_dict = unambiguous_coverage_dict[trnaseq_sample_name]
                    except:
                        # No input genome whatsoever has coverage of unambiguous seeds in the
                        # tRNA-seq sample -- the ambiguous seed has coverage in the sample and is
                        # linked to genomes.
                        continue
                    for genome_id in sample_df.groupby('genome_name').groups.keys():
                        try:
                            unambiguous_coverage = genome_coverage_dict[genome_id]
                        except KeyError:
                            # The genome does not have coverage of unambiguous seeds in the tRNA-seq
                            # sample -- there is coverage of the ambiguous seed linked to the
                            # genome.
                            unambiguous_coverage = 0
                        select_genome_coverage_dict[genome_id] = unambiguous_coverage
                # Second, for each select tRNA-seq sample, find the genome with the highest
                # unambiguous coverage. If this exceeds the second-highest unambiguous coverage
                # by the threshold ratio for the same genome in every sample, then the genome is
                # chosen as the source of the ambiguous seed.
                samples_lacking_unambiguous_coverage = []
                chosen_genome_id = ''
                for trnaseq_sample_name, select_genome_coverage_dict in \
                    select_unambiguous_coverage_dict.items():
                    if len(select_genome_coverage_dict) == 0:
                        samples_lacking_unambiguous_coverage.append(trnaseq_sample_name)
                        continue
                    select_genome_unambiguous_coverages = sorted(
                        select_genome_coverage_dict.items(), key=lambda item: -item[1])
                    first_unambiguous_coverage = select_genome_unambiguous_coverages[0][1]
                    second_unambiguous_coverage = select_genome_unambiguous_coverages[0][2]
                    if second_unambiguous_coverage == 0:
                        sample_chosen_genome_id = select_genome_unambiguous_coverages[0][0]
                    elif first_unambiguous_coverage / second_unambiguous_coverage >= \
                        self.min_coverage_ratio:
                        sample_chosen_genome_id = select_genome_unambiguous_coverages[0][0]
                    else:
                        # The threshold ratio is not met in the sample, so no genome can be chosen
                        # for the ambiguous seed.
                        sample_chosen_genome_id = None
                    if sample_chosen_genome_id != chosen_genome_id:
                        return ambiguous_seed_df[: 0]
                    if not chosen_genome_id:
                        chosen_genome_id = sample_chosen_genome_id
                return ambiguous_seed_df.set_index('genome_name').loc[
                    chosen_genome_id].reset_index()

            ambiguous_hits_df = trna_gene_hits_df[
                trna_gene_hits_df['seed_gene_callers_id'].isin(ambiguous_seed_gene_callers_ids)]
            # Retain entries for ambiguous seeds for which a genome could be chosen.
            ambiguous_hits_df = ambiguous_hits_df.groupby('seed_gene_callers_id').filter(
                choose_genome)
            trna_gene_hits_df = pd.concat(
                [unambiguous_hits_df, ambiguous_hits_df], ignore_index=True)
        ##################################################

        # Group seeds in each genome by isoacceptor. This involves summing isoacceptor seed
        # coverages in each tRNA-seq sample.
        trna_gene_hits_df = trna_gene_hits_df.set_index(
            ['seed_gene_callers_id', 'seed_contig_name'])
        isoacceptors_df = trna_gene_hits_df.groupby([
            'genome_name',
            'decoded_amino_acid',
            'anticodon',
            'effective_wobble_nucleotide',
            'trnaseq_sample_name'],
            as_index=False).aggregate('sum')

        ##################################################
        # Filter isoacceptors.

        # Drop isoacceptors with excluded anticodons.
        if self.exclude_unmodified_anticodons:
            isoacceptors_df = isoacceptors_df[
                ~isoacceptors_df['anticodon'].isin(self.exclude_unmodified_anticodons)]
        if self.exclude_modified_anticodons:
            isoacceptors_df['effective_anticodon'] = (
                isoacceptors_df['effective_wobble_nucleotide'] +
                isoacceptors_df['anticodon'].str[1: ])
            isoacceptors_df = isoacceptors_df[
                ~isoacceptors_df['effective_anticodon'].isin(self.exclude_modified_anticodons)]
            isoacceptors_df = isoacceptors_df.drop('effective_anticodon', axis=1)

        # Drop isoacceptors without data for the reference sample. Development note: It now looks
        # like the reference sample will always have an entry for the isoacceptor in the genome,
        # even with zero coverage, and therefore this filter is redundant with the next minimum
        # coverage threshold filter. I am conservatively retaining it in case I'm missing something,
        # since I now don't know why I decided to include this filter in the first place.
        isoacceptors_df = isoacceptors_df.groupby(
            ['genome_name', 'decoded_amino_acid', 'anticodon']).filter(
                lambda genome_isoacceptor_df: self.reference_sample_name in genome_isoacceptor_df[
                    'trnaseq_sample_name'].values)

        # Drop isoacceptors that do not meet the minimum coverage threshold in the reference sample.
        isoacceptors_df = isoacceptors_df.groupby(
            ['genome_name', 'decoded_amino_acid', 'anticodon']).filter(
                lambda genome_isoacceptor_df: genome_isoacceptor_df[
                    genome_isoacceptor_df['trnaseq_sample_name'] == self.reference_sample_name][
                        'discriminator_1'] >= self.min_coverage)

        # Drop isoacceptors that have coverage in only one sample, preventing the isoacceptors from
        # contributing to affinity.
        isoacceptors_df = isoacceptors_df.groupby(
            ['genome_name', 'decoded_amino_acid', 'anticodon']).filter(
                 lambda genome_isoacceptor_df: genome_isoacceptor_df[
                     'trnaseq_sample_name'].nunique() > 1)

        # Drop isoacceptor data for tRNA-seq samples in which the isoacceptor does not meet the
        # minimum coverage threshold.
        isoacceptors_df = isoacceptors_df[isoacceptors_df['discriminator_1'] >= self.min_coverage]

        # Drop data for genomes in tRNA-seq samples lacking a minimum isoacceptor diversity.
        isoacceptors_df = isoacceptors_df.groupby(['genome_name', 'trnaseq_sample_name']).filter(
            lambda genome_sample_df:
                genome_sample_df['anticodon'].nunique() > self.min_isoacceptors)

        # Drop genomes altogether if the reference sample lacked a minimum isoacceptor diversity.
        isoacceptors_df = isoacceptors_df.groupby('genome_name').filter(
            lambda genome_df: self.reference_sample_name in genome_df['trnaseq_sample_name'].values)
        ##################################################

        # Create a table of isoacceptor non-reference/reference abundance ratios.
        isoacceptor_abund_ratios_rows = []
        for genome_id, genome_df in isoacceptors_df.groupby('genome_name'):
            reference_sample_df = genome_df[
                genome_df['trnaseq_sample_name'] == self.reference_sample_name]
            if len(reference_sample_df) == 0:
                # The genome is not detected in the reference sample, precluding affinity
                # calculations.
                continue
            reference_sample_df = reference_sample_df.set_index(
                ['decoded_amino_acid', 'anticodon', 'effective_wobble_nucleotide'])

            nonreference_samples_df = genome_df[
                genome_df['trnaseq_sample_name'] != self.reference_sample_name]
            if len(nonreference_samples_df) == 0:
                # The genome is not detected in non-reference samples, precluding affinity
                # calculations.
                continue

            for trnaseq_sample_name, nonreference_sample_df in nonreference_samples_df.groupby(
                'trnaseq_sample_name'):
                for row in nonreference_sample_df.itertuples(index=False):
                    decoding_key = (
                        row.decoded_amino_acid, row.anticodon, row.effective_wobble_nucleotide)

                    try:
                        reference_isoacceptor_series = reference_sample_df.loc[decoding_key]
                    except KeyError:
                        # The isoacceptor is not detected in the reference sample, so the
                        # abundance ratio is infinite.
                        isoacceptor_abund_ratios_rows.append(
                            (genome_id, ) + decoding_key + (trnaseq_sample_name, np.nan))
                        continue
                    reference_abundance = reference_isoacceptor_series[
                        'relative_discriminator_coverage']

                    isoacceptor_abundance_ratio = \
                        row.relative_discriminator_coverage / reference_abundance
                    isoacceptor_abund_ratios_rows.append(
                        (genome_id, ) +
                        decoding_key +
                        (trnaseq_sample_name, isoacceptor_abundance_ratio))
        isoacceptor_abund_ratios_df = pd.DataFrame(
            isoacceptor_abund_ratios_rows,
            columns=['genome_name',
                     'decoded_amino_acid',
                     'anticodon',
                     'effective_wobble_nucleotide',
                     'nonreference_trnaseq_sample_name',
                     'abundance_ratio'])

        return isoacceptor_abund_ratios_df


    def get_weighted_codon_frequencies(self, isoacceptor_abund_ratios_df):
        """
        Get a table of weighted codon frequencies per isoacceptor for the functions or genes used in
        the affinity calculations. A weighted codon frequency for an isoacceptor is the sum of codon
        frequencies weighted by decoding efficiencies given the anticodon and wobble modification.

        Parameters
        ==========
        isoacceptor_abund_ratios_df : pandas.core.frame.DataFrame
            Each row of this table, returned by `get_isoacceptors`, contains genome and isoacceptor
            information identifying non-reference/reference abundance ratios. This table is only
            used to retrieve the isoacceptors identified in the tRNA-seq data: abundance data is not
            used in any way by this method.

        Returns
        =======
        weighted_frequencies_df : pandas.core.frame.DataFrame
            This table of weighted codon frequencies has rows representing functions (or genes) in
            input genomes and columns representing each isoacceptor identified in the tRNA-seq data
            regardless of genome source.
        """
        args = argparse.Namespace()
        args.contigs_db = self.genomic_contigs_db_path
        args.profile_db = self.genomic_profile_db_path
        args.collection_name = self.collection_name
        args.bin_id = self.bin_id
        args.internal_genomes = self.internal_genomes_path
        args.external_genomes = self.external_genomes_path
        args.function_sources = self.function_sources
        args.shared_function_sources = not self.lax_function_sources
        args.ignore_start_codons = True

        if self.genomic_contigs_db_path:
            codon_usage = codonusage.SingleGenomeCodonUsage(
                args, r=self.run, rq=self.run_quiet, p=self.progress)
            codon_frequency_df = codon_usage.get_frequencies(
                from_function_sources=not self.gene_affinity,
                return_functions=not self.gene_affinity,
                gene_caller_ids=self.gene_caller_ids,
                function_accessions=self.function_accessions_dict,
                function_names=self.function_names_dict,
                gene_min_codons=self.gene_min_total_codons,
                function_min_codons=self.function_min_total_codons,
                drop_amino_acids=self.exclude_amino_acids)
            # Match single- and multi-genome table formats by adding genome name to the index.
            new_index_cols = codon_frequency_df.index.names
            codon_frequency_df = codon_frequency_df.reset_index()
            codon_frequency_df['genome_name'] = list(self.genome_info_dict)[0]
            codon_frequency_df = codon_frequency_df.set_index(new_index_cols)
        else:
            codon_usage = codonusage.MultiGenomeCodonUsage(
                args, r=self.run, rq=self.run_quiet, p=self.progress)
            codon_frequency_df = codon_usage.get_frequencies(
                from_function_sources=not self.gene_affinity,
                return_functions=not self.gene_affinity,
                function_accessions=self.function_accessions_dict,
                function_names=self.function_names_dict,
                gene_min_codons=self.gene_min_total_codons,
                function_min_codons=self.function_min_total_codons,
                drop_amino_acids=self.exclude_amino_acids)

        for pattern in self.function_blacklist_patterns:
            codon_frequency_df = codon_frequency_df[
                ~codon_frequency_df['function_name'].str.contains(pattern)]

        if self.exclude_codons:
            # Exclusion of select individual codons (`exclude_codons`) is not an option in the
            # `codonusage` `get_frequencies` methods: apply this filter here.
            codon_frequency_df = codon_frequency_df.drop(self.exclude_codons, axis=1)
            if self.gene_min_total_codons:
                codon_frequency_df = codon_frequency_df[
                    codon_frequency_df.sum(axis=1) >= self.gene_min_total_codons]
            elif self.function_min_total_codons:
                codon_frequency_df = codon_frequency_df[
                    codon_frequency_df.sum(axis=1) >= self.function_min_total_codons]

        # Weight codon frequencies by decoding efficiency. For each isoacceptor found from the
        # tRNA-seq data (regardless of genome) sum the weighted frequencies of decoded codons in the
        # functions (or genes). Produce a table of function (or gene) x anticodon, with the values
        # being summed weighted codon counts.
        isoacceptors_df = isoacceptor_abund_ratios_df[
            ['decoded_amino_acid', 'anticodon', 'effective_wobble_nucleotide']].drop_duplicates()
        # Disregard tRNA-iMet and tRNA-fMet. Start codons were excluded from the codon frequencies
        # in `codon_frequency_df`.
        isoacceptors_df = isoacceptors_df[
            ~isoacceptors_df['decoded_amino_acid'].isin(['iMet', 'fMet'])]
        col_dict = {}
        for row in isoacceptors_df.itertuples(index=False):
            anticodon_wobble_nucleotide = row.effective_wobble_nucleotide
            decoding_weights_series = self.decoding_weights_df.loc[anticodon_wobble_nucleotide]
            effective_anticodon = anticodon_wobble_nucleotide + row.anticodon[1: ]
            # Get summed weighted codon counts for the isoacceptor across all functions or genes.
            summed_weighted_codon_counts = pd.Series(0, index=codon_frequency_df.index)
            for codon in self.nucleotide_decoding_dict[effective_anticodon]:
                decoding_weight = decoding_weights_series.loc[codon[2]]
                summed_weighted_codon_counts = \
                    summed_weighted_codon_counts + (1 - decoding_weight) * codon_frequency_df[codon]
                col_dict[effective_anticodon] = summed_weighted_codon_counts
        weighted_frequencies_df = pd.DataFrame.from_dict(col_dict)
        weighted_frequencies_df.index = codon_frequency_df.index

        return weighted_frequencies_df


    def get_affinities(self, isoacceptor_abund_ratios_df, weighted_frequencies_df):
        """
        Calculate affinities of tRNA isoacceptors to functions or genes for each genome.

        Parameters
        ==========
        isoacceptor_abund_ratios_df : pandas.core.frame.DataFrame
            Each row of this table, returned by `get_isoacceptors`, contains genome and isoacceptor
            information identifying non-reference/reference abundance ratios.
        weighted_frequencies_df : pandas.core.frame.DataFrame
            This table of weighted codon frequencies has rows representing functions (or genes) in
            input genomes and columns representing each isoacceptor identified in the tRNA-seq data
            regardless of genome source.

        Returns
        =======
        affinities_df : pandas.core.frame.DataFrame
            The row index columns of the affinity table are, in order, genome name, function source,
            function accession, function name -- or, if genes are analyzed, there are no function
            indices, but a gene callers ID index column. Columns are isoacceptor anticodons, with
            modified wobble nucleotide if applicable (e.g., ICG, LAT).
        """
        isoacceptor_abund_ratios_df['effective_anticodon'] = (
            isoacceptor_abund_ratios_df['effective_wobble_nucleotide'] +
            isoacceptor_abund_ratios_df['anticodon'].str[1: ])

        isoacceptor_abund_ratios_gb = isoacceptor_abund_ratios_df.groupby('genome_name')
        weighted_frequencies_gb = weighted_frequencies_df.groupby('genome_name')

        filtered_genome_names = []
        genome_affinities_dfs = []
        for genome_name in self.genome_info_dict:
            try:
                genome_isoacceptor_abund_ratios_df = isoacceptor_abund_ratios_gb.get_group(
                    genome_name)
            except KeyError:
                filtered_genome_names.append(genome_name)
                continue
            try:
                genome_weighted_frequencies_df = weighted_frequencies_gb.get_group(genome_name)
            except KeyError:
                filtered_genome_names.append(genome_name)
                continue

            sample_affinities_dict = {}
            for trnaseq_sample_name, sample_isoacceptor_abund_ratios_df in \
                genome_isoacceptor_abund_ratios_df.groupby('nonreference_trnaseq_sample_name'):
                # Take the dot product of the m function or gene x n isoacceptor matrix of
                # isoacceptor weighted codon frequencies and the n x 1 matrix (Series) of
                # non-reference/reference sample isoacceptor abundance ratios, yielding the affinity
                # of the tRNA pool to each function or gene in the genome.
                abund_ratios = sample_isoacceptor_abund_ratios_df[
                    ['effective_anticodon', 'abundance_ratio']].set_index('effective_anticodon')[
                        'abundance_ratio']
                # Every "effective" anticodon (including modified wobble nucleotide, if applicable)
                # from the tRNA-seq data should be represented by a column in the weighted codon
                # frequencies table output by `get_weighted_codon_frequencies`.
                sample_affinities_dict[trnaseq_sample_name] = \
                    genome_weighted_frequencies_df[abund_ratios.index].dot(abund_ratios)

            genome_affinities_df = pd.DataFrame.from_dict(sample_affinities_dict)
            genome_affinities_df['genome_name'] = genome_name
            genome_affinities_dfs.append(genome_affinities_df)

        if filtered_genome_names:
            self.run.info_single("The following genomes did not pass the filters for affinity "
                                 f"calculation: {', '.join(filtered_genome_names)}")

        if affinities_df:
            affinities_df = pd.concat(genome_affinities_dfs, axis=0)
        else:
            affinities_df = pd.DataFrame()

        return affinities_df


    @staticmethod
    def list_sample_names(seeds_specific_txt_path):
        """List samples in the tRNA-seq input files."""
        if seeds_specific_txt_path is None:
            raise ConfigError(
                "To list samples in `seeds_specific_txt`, a path to this file must be provided.")
        filesnpaths.is_file_exists(seeds_specific_txt_path)

        available_sample_names = pd.read_csv(
            seeds_specific_txt_path,
            sep='\t',
            header=0,
            skiprows=[1, 2],
            usecols=['sample_name'])['sample_name'].unique().tolist()

        return available_sample_names

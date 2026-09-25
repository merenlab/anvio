"""Convert genomic codon composition into tRNA isoacceptor (anticodon-level) decoding demand.

Unlike `anvio.genomictrnaseq.Affinitizer`, which infers tRNA wobble modification state from
tRNA-seq sequencing evidence, this module infers it from the presence or absence of the genes
that encode the modifying enzymes (e.g., a tadA-like adenosine deaminase converting A34 to
inosine, or a TilS-like enzyme converting C34 to lysidine at tRNA-Ile2). This lets it work from a
contigs database alone, with no tRNA-seq data required.

The decoding-weights table and the anticodon-to-codon translation logic below are intentionally
duplicated from `anvio.genomictrnaseq.Affinitizer` rather than imported, so that this module has
no dependency on the tRNA-seq-specific `genomictrnaseq` module.
"""

import argparse
import pandas as pd

import anvio
import anvio.db as db
import anvio.tables as t
import anvio.terminal as terminal
import anvio.constants as constants
import anvio.codonusage as codonusage
import anvio.ccollections as ccollections
import anvio.filesnpaths as filesnpaths
import anvio.hmmops as hmmops

from anvio.dbinfo import DBInfo
from anvio.errors import ConfigError


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__


run = terminal.Run()
progress = terminal.Progress()
run_quiet = terminal.Run(verbose=False)


# Bacterial mean wobble s(i,j) values from Table 4 of Sabi and Tuller (2014,
# https://doi.org/10.1093/dnares/dsu017). Row index is the anticodon wobble nucleotide (position
# 34; I=inosine derived from A, L=lysidine derived from C); column headers are the codon wobble
# nucleotide (position 3). A cell value of 1 means the anticodon cannot decode that codon.
default_decoding_weights_df = pd.DataFrame([
    [1, 1, 1, 0],
    [1, 1, 0, 1],
    [1, 0, 1, 0.6294],
    [0, 1, 0.698, 1],
    [0.8773, 0.4211, 1, 0],
    [0.7309, 1, 1, 1]],
    index=['A', 'C', 'G', 'T', 'I', 'L'],
    columns=['A', 'C', 'G', 'T'])

recognized_anticodon_wobble_modifications = ['I', 'L']

TRNA_MODIFICATION_ENZYME_LIST_COLUMNS = [
    'modifying_enzyme_name', 'modification', 'canonical_position', 'expected_reference',
    'isoacceptor_specificity', 'aa', 'anticodon', 'function_accession', 'function_source']

# The only two wobble modifications this program's decoding-weight model understands. Other
# enzyme-list rows (e.g. miaA/i6A at position 37) are valid entries but not actionable here.
MODIFICATION_TO_WOBBLE_LETTER = {
    'inosine': 'I', 'I': 'I',
    'lysidine': 'L', 'L': 'L', 'k2C': 'L'}


def build_nucleotide_decoding_dict(decoding_weights_df):
    """Map each anticodon (plain and wobble-modified forms) to the codons it can decode.

    Parameters
    ==========
    decoding_weights_df : pandas.core.frame.DataFrame
        Shaped like `default_decoding_weights_df`.

    Returns
    =======
    nucleotide_decoding_dict : dict
        Maps an anticodon string (e.g. 'AAT', or a derived modified form like 'IAT'/'LAT') to a
        list of codon strings it can decode. Covers the full anticodon universe in
        `constants.anticodon_to_AA`, not just anticodons known to exist in any particular genome.
    """
    nucleotide_decoding_dict = {}
    for anticodon in constants.anticodon_to_AA:
        anticodon_wobble_nucleotide = anticodon[0]
        # Two anticodon wobble modifications are recognized: ANN -> INN and CAT -> LAT.
        derived_anticodons = [anticodon]
        if anticodon_wobble_nucleotide == 'A':
            derived_anticodons.append('I' + anticodon[1: ])
        elif anticodon == 'CAT':
            derived_anticodons.append('LAT')

        for derived_anticodon in derived_anticodons:
            derived_anticodon_wobble_nucleotide = derived_anticodon[0]
            decoding_weights_series = decoding_weights_df.loc[derived_anticodon_wobble_nucleotide]
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
            nucleotide_decoding_dict[derived_anticodon] = decoded_codons

    return nucleotide_decoding_dict


def parse_decoding_weights_table(decoding_weights_txt):
    """Load a user-supplied decoding weights table, validating it against the default shape."""
    if not decoding_weights_txt:
        return None

    filesnpaths.is_file_tab_delimited(decoding_weights_txt)

    decoding_weights_df = pd.read_csv(decoding_weights_txt, sep='\t', index_col=0, header=0)

    is_formatted = True
    if decoding_weights_df.shape != default_decoding_weights_df.shape:
        is_formatted = False
    else:
        if not (decoding_weights_df.index == default_decoding_weights_df.index).all():
            is_formatted = False
        if not (decoding_weights_df.columns == default_decoding_weights_df.columns).all():
            is_formatted = False
    if not is_formatted:
        raise ConfigError(
            f"The input decoding weights table at '{decoding_weights_txt}' is not formatted "
            f"properly. The tab-delimited file generated by running this program with the sole "
            f"argument, `--get-default-decoding-weights`, should be used as a template, with the "
            f"row index of anticodon wobble nucleotides and column header of codon wobble "
            f"nucleotides left unchanged.")

    return decoding_weights_df


def parse_trna_modification_enzyme_list(enzyme_list_path):
    """Load and validate a %(trna-modification-enzyme-list)s file.

    Returns
    =======
    enzyme_list_df : pandas.core.frame.DataFrame or None
        None if `enzyme_list_path` is falsy, signaling that no enzyme-presence-driven wobble
        modification should be inferred at all.
    """
    if not enzyme_list_path:
        return None

    filesnpaths.is_file_exists(enzyme_list_path)
    filesnpaths.is_file_tab_delimited(enzyme_list_path)

    enzyme_list_df = pd.read_csv(enzyme_list_path, sep='\t')

    missing_columns = [
        c for c in TRNA_MODIFICATION_ENZYME_LIST_COLUMNS if c not in enzyme_list_df.columns]
    if missing_columns:
        raise ConfigError(
            f"The tRNA modification enzyme list at '{enzyme_list_path}' is missing the "
            f"following required column(s): {', '.join(missing_columns)}.")

    try:
        enzyme_list_df['canonical_position'] = enzyme_list_df['canonical_position'].astype(int)
    except (ValueError, TypeError) as e:
        raise ConfigError(
            f"The `canonical_position` column of the tRNA modification enzyme list at "
            f"'{enzyme_list_path}' must contain only integers. Pandas complained with the "
            f"following error: {e}")

    for row in enzyme_list_df.itertuples(index=False):
        wobble_letter = MODIFICATION_TO_WOBBLE_LETTER.get(row.modification)
        if wobble_letter is None:
            # Not a modification this program's decoding-weight model can act on. Nothing to
            # validate here -- see `resolve_effective_anticodons`, which will skip and warn about
            # this row at runtime.
            continue

        anticodon = str(row.anticodon).upper()
        if wobble_letter == 'L' and anticodon != 'CAT':
            raise ConfigError(
                f"In the tRNA modification enzyme list at '{enzyme_list_path}', the row for "
                f"`{row.modifying_enzyme_name}` has a `modification` value ('{row.modification}') "
                f"that implies a C34-to-lysidine wobble modification, but its `anticodon` column "
                f"is '{row.anticodon}', not 'CAT'. Only tRNA-Ile2 (anticodon CAT) is known to "
                f"carry this modification -- please check this row.")
        if wobble_letter == 'I' and (not anticodon or anticodon[0] != 'A'):
            raise ConfigError(
                f"In the tRNA modification enzyme list at '{enzyme_list_path}', the row for "
                f"`{row.modifying_enzyme_name}` has a `modification` value ('{row.modification}') "
                f"that implies an A34-to-inosine wobble modification, but its `anticodon` column "
                f"is '{row.anticodon}', which does not start with 'A'. Please check this row.")

    return enzyme_list_df


def get_genome_trna_gene_calls(contigs_db_path, split_names_of_interest=set(), run=run_quiet):
    """Get the tRNA gene calls in a genome (or a bin within it) from `hmm_hits`.

    Requires `anvi-scan-trnas` (or `anvi-run-hmms --also-scan-trnas`) to have been run on the
    contigs database. Raises a `ConfigError` (via `hmmops.SequencesForHMMHits`) if it hasn't.

    Returns
    =======
    trna_gene_calls_df : pandas.core.frame.DataFrame
        Columns: `gene_callers_id`, `aa_type` (the amino acid or isotype label tRNAscan-SE
        assigned, e.g. 'Arg', 'Ile2'), `anticodon` (the raw, unmodified anticodon sequence).
    """
    sequences_for_hmm_hits = hmmops.SequencesForHMMHits(
        contigs_db_path,
        sources={'Transfer_RNAs'},
        split_names_of_interest=split_names_of_interest,
        load_sequences=False,
        run=run)

    # `hmm_hits` is a dict of dicts keyed by entry_id. The amino-acid/isotype label and raw
    # anticodon are recovered from the composite `gene_name`, formatted as f'{aa_type}_{anticodon}'
    # by `anvio.tables.trnahits.TablesForTransferRNAs`.
    hits_df = pd.DataFrame.from_dict(sequences_for_hmm_hits.hmm_hits, orient='index')
    if hits_df.empty:
        return pd.DataFrame(columns=['gene_callers_id', 'aa_type', 'anticodon'])

    aa_types, anticodons = zip(*(name.rsplit('_', 1) for name in hits_df['gene_name']))
    hits_df = hits_df.assign(aa_type=aa_types, anticodon=anticodons)

    return (hits_df[['gene_callers_id', 'aa_type', 'anticodon']]
            .drop_duplicates()
            .reset_index(drop=True))


def get_genome_gene_functions(contigs_db_path, gene_caller_ids_of_interest=None):
    """Get the `gene_functions` table of a contigs database, optionally restricted to certain genes.

    Parameters
    ==========
    gene_caller_ids_of_interest : set or None
        When given, restricts the returned rows to these gene caller IDs -- used to scope the
        enzyme-presence check to a single bin rather than an entire (meta)genome assembly.

    Returns
    =======
    gene_functions_df : pandas.core.frame.DataFrame
        Columns: `gene_callers_id`, `source`, `accession`.
    """
    contigs_db = db.DB(contigs_db_path, None, ignore_version=True)
    gene_functions_df = contigs_db.get_table_as_dataframe(
        t.gene_function_calls_table_name,
        columns_of_interest=['gene_callers_id', 'source', 'accession'],
        error_if_no_data=False)

    if gene_caller_ids_of_interest is not None and len(gene_functions_df):
        gene_functions_df = gene_functions_df[
            gene_functions_df['gene_callers_id'].isin(gene_caller_ids_of_interest)]

    contigs_db.disconnect()

    return gene_functions_df


def get_gene_caller_ids_in_splits(contigs_db_path, split_names_of_interest):
    """Resolve a set of split names to the gene caller IDs they contain."""
    if not split_names_of_interest:
        return None

    contigs_db = db.DB(contigs_db_path, None, ignore_version=True)
    genes_in_splits_df = contigs_db.get_table_as_dataframe(
        t.genes_in_splits_table_name, columns_of_interest=['split', 'gene_callers_id'])
    contigs_db.disconnect()

    return set(genes_in_splits_df[
        genes_in_splits_df['split'].isin(split_names_of_interest)]['gene_callers_id'])


def resolve_effective_anticodons(trna_gene_calls_df, gene_functions_df, enzyme_list_df, run=run_quiet):
    """Map each (amino acid/isotype, raw anticodon) pair present in a genome to its effective anticodon.

    An isoacceptor absent from `trna_gene_calls_df` never enters the returned mapping at all --
    that is how an isoacceptor's tRNA gene being entirely absent from a genome is distinguished
    from it being present but unmodified.

    Returns
    =======
    effective_anticodons : dict
        Maps `(aa_type, raw_anticodon)` to the effective anticodon string, which is the raw
        anticodon unchanged unless an enzyme list identifies it as wobble-modified in this genome
        (in which case the wobble nucleotide is replaced with 'I' or 'L').
    """
    effective_anticodons = {
        (row.aa_type, row.anticodon): row.anticodon
        for row in trna_gene_calls_df.itertuples(index=False)}

    if enzyme_list_df is None or not len(effective_anticodons):
        return effective_anticodons

    for row in enzyme_list_df.itertuples(index=False):
        wobble_letter = MODIFICATION_TO_WOBBLE_LETTER.get(row.modification)
        if row.canonical_position != 34 or wobble_letter is None:
            run.warning(
                f"The tRNA modification enzyme list has an entry for `{row.modifying_enzyme_name}` "
                f"(modification: '{row.modification}', canonical position: {row.canonical_position}) "
                f"that this program does not know how to translate into a wobble-nucleotide change, "
                f"so it will be ignored. Only modifications resolving to inosine ('I') or lysidine "
                f"('L') at canonical position 34 affect isoacceptor decoding in this program.",
                header="Enzyme list entry not actionable", lc='yellow')
            continue

        if wobble_letter == 'L':
            # Lysidine only ever occurs at tRNA-Ile2 (anticodon CAT); no other CAT-anticodon
            # isoacceptor (Met, fMet, iMet) is a valid substrate, regardless of what the enzyme
            # list's `isoacceptor_specificity` column claims.
            candidates = [
                key for key in effective_anticodons
                if key[0] == 'Ile2' and key[1].upper() == 'CAT']
        elif row.isoacceptor_specificity == 'Specific':
            candidates = [
                key for key in effective_anticodons
                if key[0].upper() == str(row.aa).upper()
                and key[1].upper() == str(row.anticodon).upper()]
        else:
            # Non-specific A34 -> inosine row: apply to every A-wobble isoacceptor present.
            candidates = [key for key in effective_anticodons if key[1][0] == 'A']

        if not candidates:
            continue

        enzyme_present = len(gene_functions_df) and (
            (gene_functions_df['source'] == row.function_source) &
            (gene_functions_df['accession'] == row.function_accession)).any()
        if not enzyme_present:
            continue

        for aa_type, anticodon in candidates:
            effective_anticodons[(aa_type, anticodon)] = wobble_letter + anticodon[1:]

    return effective_anticodons


class TRNAIsoacceptorFrequencyCalculator:
    """Converts a genome's codon composition into tRNA isoacceptor decoding demand."""

    def __init__(self, args=argparse.Namespace(), r=run, rq=run_quiet, p=progress):
        A = lambda x: args.__dict__[x] if x in args.__dict__ else None

        self.args = args
        self.genomic_contigs_db_path = A('contigs_db')
        self.genomic_profile_db_path = A('profile_db')
        self.collection_name = A('collection_name')
        self.bin_id = A('bin_id')
        self.internal_genomes_path = A('internal_genomes')
        self.external_genomes_path = A('external_genomes')
        self.gene_caller_ids = A('gene_caller_ids')
        self.function_sources = A('function_sources')
        self.all_brite_categories = A('all_brite_categories')
        self.shared_function_sources = A('shared_function_sources')
        self.enzyme_list_path = A('trna_modification_enzyme_list')

        self.decoding_weights_df = A('decoding_weights')
        if self.decoding_weights_df is None:
            self.decoding_weights_df = default_decoding_weights_df
        self.nucleotide_decoding_dict = build_nucleotide_decoding_dict(self.decoding_weights_df)
        self.enzyme_list_df = parse_trna_modification_enzyme_list(self.enzyme_list_path)

        self.run = r
        self.run_quiet = rq
        self.progress = p

        if self.enzyme_list_df is None:
            self.run.info_single(
                "No tRNA modification enzyme list provided -- every isoacceptor will be treated "
                "as canonically unmodified.")
        else:
            self.run.info("tRNA modification enzyme list", self.enzyme_list_path)


    def _get_single_genome_info_dict(self):
        """Build a one-entry genome_info_dict for the single-contigs-db case.

        Mirrors the genome-naming convention of `anvio.genomictrnaseq.Affinitizer`: a bin is keyed
        by its bin ID, and a plain contigs database is keyed by its project name.
        """
        if self.bin_id:
            genome_name = self.bin_id
        else:
            genome_name = DBInfo(self.genomic_contigs_db_path, expecting='contigs').project_name

        return {genome_name: {
            'contigs_db': self.genomic_contigs_db_path,
            'profile_db': self.genomic_profile_db_path,
            'collection_name': self.collection_name,
            'bin_id': self.bin_id}}


    def _get_split_names_of_interest(self, genome_info):
        """Resolve the splits belonging to a bin, or an empty set for a whole contigs database."""
        if not genome_info.get('bin_id'):
            return set()

        bin_args = argparse.Namespace(
            contigs_db=genome_info['contigs_db'],
            profile_db=genome_info['profile_db'],
            collection_name=genome_info['collection_name'],
            bin_id=genome_info['bin_id'])
        return ccollections.GetSplitNamesInBins(bin_args).get_split_names_only()


    def get_frequencies(self,
                        from_function_sources=False, return_functions=False,
                        gene_caller_ids=None, function_accessions=None, function_names=None,
                        expect_functions=False, gene_min_codons=0, function_min_codons=0,
                        min_codon_filter='both', drop_amino_acids=None,
                        sequence_min_amino_acids=0, pansequence_min_amino_acids=(0, 1.0),
                        relative=False, sum_genes=False, average_genes=False,
                        deviation_from_genome_average=False, infinity_to_zero=False):
        """Get a table of tRNA isoacceptor decoding demand, at gene or function level.

        Parameters
        ==========
        See `anvio.codonusage.SingleGenomeCodonUsage.get_frequencies`/`MultiGenomeCodonUsage.get_frequencies`
        for the codon-composition-related parameters, which are passed through unchanged.

        relative : bool
            Normalize each row to proportions of that row's total weighted decoding demand.
        sum_genes, average_genes : bool
            Collapse to a single row per genome (sum or mean across genes/functions).
        deviation_from_genome_average : bool
            Report each row as a fold-change against that genome's average isoacceptor profile.
            Requires `relative=True` -- otherwise a longer gene or function would appear to
            deviate from the genome average on every isoacceptor purely because it has more
            codons, not because its isoacceptor usage is actually biased. Not combinable with
            `sum_genes`/`average_genes` (that combination is degenerate -- every value would be
            1.0).
        infinity_to_zero : bool
            Under `deviation_from_genome_average`, replace any resulting +/-inf (a nonzero value
            divided by a zero genome-average) with 0.0, and also replace a 0/0 quotient (an
            isoacceptor whose value and genome average were both exactly zero) with 0.0. An
            isoacceptor that is genuinely absent from a genome (`NaN` prior to the division)
            is left untouched either way -- this only ever resolves a real division artifact.

        Returns
        =======
        isoacceptor_frequency_df : pandas.core.frame.DataFrame
            Index: `['genome_name', 'gene_caller_id']` (gene-level) or
            `['genome_name', 'function_source', 'function_accession', 'function_name']`
            (function-level), collapsed to just `['genome_name']` if `sum_genes`/`average_genes`
            is set without `deviation_from_genome_average`.
            Columns: sorted `f'{aa_type}_{effective_anticodon}'` strings.
            Values: absolute weighted decoding demand by default; proportions if `relative`;
            fold-change vs. genome average if `deviation_from_genome_average`. `NaN` marks an
            isoacceptor whose tRNA gene does not exist in a given genome (multi-genome only).
        """
        if deviation_from_genome_average and (sum_genes or average_genes):
            raise ConfigError(
                "`--deviation-from-genome-average` cannot be combined with `--sum`/`--average`: "
                "once genes or functions are collapsed to a single row per genome, that row *is* "
                "the genome average, so every deviation value would trivially be 1.0. Please use "
                "one or the other.")

        if deviation_from_genome_average and not relative:
            raise ConfigError(
                "`--deviation-from-genome-average` requires `--relative`. Absolute weighted "
                "decoding demand scales with how many codons a gene or function has, so without "
                "`--relative`, a longer gene would appear to deviate from the genome average on "
                "every isoacceptor just by virtue of its length, not because its isoacceptor "
                "usage is actually biased. Please add `--relative`.")

        codon_usage_args = argparse.Namespace()
        codon_usage_args.contigs_db = self.genomic_contigs_db_path
        codon_usage_args.profile_db = self.genomic_profile_db_path
        codon_usage_args.collection_name = self.collection_name
        codon_usage_args.bin_id = self.bin_id
        codon_usage_args.gene_caller_ids = gene_caller_ids
        codon_usage_args.internal_genomes = self.internal_genomes_path
        codon_usage_args.external_genomes = self.external_genomes_path
        codon_usage_args.function_sources = self.function_sources
        codon_usage_args.all_brite_categories = self.all_brite_categories
        codon_usage_args.shared_function_sources = self.shared_function_sources
        # Start-codon-driven demand doesn't belong in elongator isoacceptor demand.
        codon_usage_args.ignore_start_codons = True
        # For now, this program only supports the standard genetic code.
        codon_usage_args.codon_to_amino_acid = None

        if self.internal_genomes_path or self.external_genomes_path:
            codon_usage = codonusage.MultiGenomeCodonUsage(codon_usage_args, r=self.run, rq=self.run_quiet, p=self.progress)
            codon_frequency_df = codon_usage.get_frequencies(
                from_function_sources=from_function_sources,
                return_functions=return_functions,
                function_accessions=function_accessions,
                function_names=function_names,
                expect_functions=expect_functions,
                gene_min_codons=gene_min_codons,
                function_min_codons=function_min_codons,
                min_codon_filter=min_codon_filter,
                drop_amino_acids=drop_amino_acids,
                sequence_min_amino_acids=sequence_min_amino_acids,
                pansequence_min_amino_acids=pansequence_min_amino_acids)
            genome_info_dict = codon_usage.genome_info_dict
        else:
            codon_usage = codonusage.SingleGenomeCodonUsage(codon_usage_args, r=self.run, rq=self.run_quiet, p=self.progress)
            codon_frequency_df = codon_usage.get_frequencies(
                from_function_sources=from_function_sources,
                return_functions=return_functions,
                gene_caller_ids=gene_caller_ids,
                function_accessions=function_accessions,
                function_names=function_names,
                expect_functions=expect_functions,
                gene_min_codons=gene_min_codons,
                function_min_codons=function_min_codons,
                min_codon_filter=min_codon_filter,
                drop_amino_acids=drop_amino_acids,
                sequence_min_amino_acids=sequence_min_amino_acids,
                pansequence_min_amino_acids=pansequence_min_amino_acids)
            genome_info_dict = self._get_single_genome_info_dict()
            # Match single- and multi-genome table formats by adding genome name to the index.
            genome_name = list(genome_info_dict)[0]
            new_index_cols = ['genome_name'] + codon_frequency_df.index.names
            codon_frequency_df = codon_frequency_df.reset_index()
            codon_frequency_df['genome_name'] = genome_name
            codon_frequency_df = codon_frequency_df.set_index(new_index_cols)

        genome_isoacceptor_dfs = []
        for genome_name, genome_info in genome_info_dict.items():
            if genome_name not in codon_frequency_df.index.get_level_values('genome_name'):
                continue

            genome_codon_frequency_df = codon_frequency_df.xs(genome_name, level='genome_name')
            split_names_of_interest = self._get_split_names_of_interest(genome_info)
            gene_caller_ids_of_interest = get_gene_caller_ids_in_splits(
                genome_info['contigs_db'], split_names_of_interest)

            trna_gene_calls_df = get_genome_trna_gene_calls(
                genome_info['contigs_db'], split_names_of_interest, run=self.run_quiet)
            gene_functions_df = get_genome_gene_functions(
                genome_info['contigs_db'], gene_caller_ids_of_interest)
            effective_anticodons = resolve_effective_anticodons(
                trna_gene_calls_df, gene_functions_df, self.enzyme_list_df, run=self.run)

            num_modified = sum(
                1 for raw, effective in
                ((anticodon, effective) for (_, anticodon), effective in effective_anticodons.items())
                if raw != effective)
            self.run.info(
                f"{genome_name}: isoacceptors found / modified",
                f"{len(effective_anticodons)} / {num_modified}")

            col_dict = {}
            for (aa_type, raw_anticodon), effective_anticodon in effective_anticodons.items():
                if effective_anticodon not in self.nucleotide_decoding_dict:
                    continue

                decoding_weights_series = self.decoding_weights_df.loc[effective_anticodon[0]]
                summed_weighted_codon_counts = pd.Series(0.0, index=genome_codon_frequency_df.index)
                for codon in self.nucleotide_decoding_dict[effective_anticodon]:
                    if codon not in genome_codon_frequency_df.columns:
                        continue
                    decoding_weight = decoding_weights_series.loc[codon[2]]
                    summed_weighted_codon_counts = (
                        summed_weighted_codon_counts +
                        (1 - decoding_weight) * genome_codon_frequency_df[codon])
                col_dict[f'{aa_type}_{effective_anticodon}'] = summed_weighted_codon_counts

            genome_isoacceptor_df = pd.DataFrame.from_dict(col_dict)
            genome_isoacceptor_df.index = genome_codon_frequency_df.index
            new_index_cols = ['genome_name'] + genome_isoacceptor_df.index.names
            genome_isoacceptor_df = genome_isoacceptor_df.reset_index()
            genome_isoacceptor_df['genome_name'] = genome_name
            genome_isoacceptor_df = genome_isoacceptor_df.set_index(new_index_cols)
            genome_isoacceptor_dfs.append(genome_isoacceptor_df)

        if not genome_isoacceptor_dfs:
            raise ConfigError("No genomes yielded any tRNA isoacceptor data.")

        isoacceptor_frequency_df = pd.concat(genome_isoacceptor_dfs, axis=0, sort=True)
        isoacceptor_frequency_df = isoacceptor_frequency_df[sorted(isoacceptor_frequency_df.columns)]

        if relative:
            isoacceptor_frequency_df = isoacceptor_frequency_df.div(
                isoacceptor_frequency_df.sum(axis=1, skipna=True), axis=0)

        # Tracks which cells held a real (non-NaN) value before the division below, so that
        # `infinity_to_zero` can later distinguish a 0/0 quotient (a genuine isoacceptor whose
        # genome-wide average happens to be exactly zero) from an isoacceptor that was already
        # `NaN` because its tRNA gene doesn't exist in that genome at all -- only the former
        # should ever be converted to 0.0.
        was_present_before_deviation = None

        if deviation_from_genome_average:
            was_present_before_deviation = isoacceptor_frequency_df.notna()
            genome_average_df = isoacceptor_frequency_df.groupby(level='genome_name').mean()
            isoacceptor_frequency_df = isoacceptor_frequency_df.div(
                genome_average_df, level='genome_name')
        elif sum_genes:
            isoacceptor_frequency_df = isoacceptor_frequency_df.groupby(level='genome_name').sum()
        elif average_genes:
            isoacceptor_frequency_df = isoacceptor_frequency_df.groupby(level='genome_name').mean()

        if infinity_to_zero:
            isoacceptor_frequency_df = isoacceptor_frequency_df.replace([float('inf'), float('-inf')], 0.0)
            if was_present_before_deviation is not None:
                # A `NaN` here can only come from a 0/0 quotient (the isoacceptor was present, but
                # both its value and its genome's average were exactly zero). Leave a genuinely
                # absent isoacceptor's `NaN` (present before the division) untouched.
                zero_over_zero = isoacceptor_frequency_df.isna() & was_present_before_deviation
                isoacceptor_frequency_df = isoacceptor_frequency_df.where(~zero_over_zero, 0.0)

        return isoacceptor_frequency_df

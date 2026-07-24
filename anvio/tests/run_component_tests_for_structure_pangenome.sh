#!/bin/bash
# Component test for the PANGENOME input type of anvi-gen-structure-database and
# anvi-update-structure-database.
#
# A structure database built from a pangenome keys every protein on a surrogate integer `protein_id`.
# This is unlike the contigs-db input, where `protein_id == gene_callers_id` (identity), because a
# pangenome gene id is only unique WITHIN a genome. The source therefore mints fresh surrogate ids and
# records the real provenance -- genome_name, gene_callers_id, gene_cluster_id -- in the `proteins`
# table. The assertions below focus on that bookkeeping and, above all, on protein_id STABILITY across
# update/rerun, which is the one behavior where the pangenome path genuinely diverges from contigs-db.
#
# The default MODELLER engine is slow and may find no template for a given single-copy-core gene, so
# the `proteins` table (rows, surrogate ids, provenance, idempotency) is asserted STRICTLY while
# whether a structure was actually produced is treated softly.
#
# Run it with no arguments (output goes to a fresh sandbox/test-output), e.g.:
#     bash run_component_tests_for_structure_pangenome.sh
# Unlike the SCVs/SAAVs structure test, this one takes no 'new'/'make'/'display' mode argument; its
# single optional first argument is an output directory (see SETUP_WITH_OUTPUT_DIR in 00.sh). It runs
# the default MODELLER engine and needs internet, because MODELLER downloads the template structures it
# finds from the RCSB PDB (there is no local --pdb-db).

source 00.sh
set -e

# ---------------------------------------------------------------------------------------------------
# assertion helpers
# ---------------------------------------------------------------------------------------------------
assert_eq() {
    # assert_eq <actual> <expected> <message>
    if [ "$1" != "$2" ]; then
        echo ""
        echo "ASSERTION FAILED: $3"
        echo "  expected: '$2'"
        echo "  actual:   '$1'"
        exit 1
    fi
    echo "  OK: $3 (== '$2')"
}

# interrogate a structure db's proteins table / self table
proteins_count() { sqlite3 "$1" "SELECT COUNT(*) FROM proteins;"; }
input_type_of()  { sqlite3 "$1" "SELECT value FROM self WHERE key='input_type';"; }
structures_count() { sqlite3 "$1" "SELECT COUNT(*) FROM structures;"; }
# number of proteins rows whose provenance fields are fully populated (genome + gene call + cluster)
proteins_fully_provenanced() {
    sqlite3 "$1" "SELECT COUNT(*) FROM proteins \
                  WHERE genome_name IS NOT NULL \
                    AND gene_callers_id IS NOT NULL \
                    AND gene_cluster_id IS NOT NULL;"
}
# stable protein_id <-> natural-key mapping, sorted, for before/after diffs
proteins_map() {
    sqlite3 -separator $'\t' "$1" \
        "SELECT protein_id, genome_name, gene_callers_id, gene_cluster_id \
         FROM proteins ORDER BY protein_id;"
}

# ---------------------------------------------------------------------------------------------------
# setup: build a small pangenome (genomes storage + pan db WITH alignments) from the mock genomes
# ---------------------------------------------------------------------------------------------------
SETUP_WITH_OUTPUT_DIR $1 $2 $3

INFO "Setting up the pangenome for the structure test"
mkdir -p $output_dir/
cp $files/mock_data_for_pangenomics/*.db                 $output_dir/
cp $files/mock_data_for_pangenomics/external-genomes.txt $output_dir/
cd $output_dir/

INFO "Migrating the mock genome databases"
anvi-migrate *.db --migrate-quickly

INFO "Generating an anvi'o genomes storage"
anvi-gen-genomes-storage -e external-genomes.txt -o TEST-GENOMES.db --no-progress

# Alignments are ON (the default) because representative-picking relies on the aligned sequences.
INFO "Running the pangenome analysis (with gene-cluster alignments)"
anvi-pan-genome -g TEST-GENOMES.db -n TEST --use-ncbi-blast --no-progress $thread_controller

# Two deterministic single-copy-core gene clusters of this pangenome. Each occurs exactly once in all
# five genomes (so all-members mode yields five proteins per cluster and representative mode yields one),
# and both have a high-identity, high-coverage template in MODELLER's pdb_95 search database (~98-99%
# identity over ~99% of the gene), so MODELLER reliably models their structures and the create step is
# not left empty. The ids are stable for these mock genomes and the default clustering, exactly as the
# pangenomics component test relies on hard-coded gene cluster ids for the same data.
GC_A=GC_00000068
GC_B=GC_00000088
NUM_GENOMES=5

# common, engine-agnostic knobs -- keep MODELLER fast; ColabFold is a separate opt-in path elsewhere.
GEN_COMMON="--pan-db TEST-PAN.db -g TEST-GENOMES.db --very-fast --num-models 1 --num-threads 1 --debug"

# ===================================================================================================
# TEST 1 -- create, all-members mode: every gene of the selected cluster becomes a protein
# ===================================================================================================
INFO "anvi-gen-structure-database from a pangenome (all members of one gene cluster)"
anvi-gen-structure-database $GEN_COMMON \
                            --gene-cluster-ids $GC_A \
                            -o STRUCTURE_PAN_ALL.db

assert_eq "$(input_type_of STRUCTURE_PAN_ALL.db)" "pangenome" \
          "structure db records input_type=pangenome"
assert_eq "$(proteins_count STRUCTURE_PAN_ALL.db)" "$NUM_GENOMES" \
          "one protein per member gene of $GC_A"
assert_eq "$(proteins_fully_provenanced STRUCTURE_PAN_ALL.db)" "$NUM_GENOMES" \
          "every protein carries genome_name + gene_callers_id + gene_cluster_id"
assert_eq "$(sqlite3 STRUCTURE_PAN_ALL.db "SELECT COUNT(DISTINCT protein_id) FROM proteins;")" "$NUM_GENOMES" \
          "surrogate protein_ids are distinct"
assert_eq "$(sqlite3 STRUCTURE_PAN_ALL.db "SELECT COUNT(DISTINCT gene_cluster_id) FROM proteins;")" "1" \
          "all members trace back to the single selected gene cluster"
# soft: structures may be fewer than proteins if MODELLER found no template for some SCG
INFO "  (soft) structures produced for all-members db: $(structures_count STRUCTURE_PAN_ALL.db)/$NUM_GENOMES"

# ===================================================================================================
# TEST 2 -- create, representative mode: one picked protein per gene cluster
# ===================================================================================================
INFO "anvi-gen-structure-database from a pangenome (representative of two gene clusters)"
anvi-gen-structure-database $GEN_COMMON \
                            --gene-cluster-ids $GC_A,$GC_B \
                            --select-representative \
                            -o STRUCTURE_PAN_REP.db

assert_eq "$(proteins_count STRUCTURE_PAN_REP.db)" "2" \
          "exactly one representative protein per selected gene cluster"
assert_eq "$(proteins_fully_provenanced STRUCTURE_PAN_REP.db)" "2" \
          "each representative resolves to a real (genome, gene_callers_id) and its cluster"
assert_eq "$(sqlite3 STRUCTURE_PAN_REP.db "SELECT COUNT(DISTINCT gene_cluster_id) FROM proteins;")" "2" \
          "the two representatives come from the two distinct clusters"

# ===================================================================================================
# TEST 3 -- update idempotency: the crux. Surrogate ids must be STABLE across add and rerun.
# ===================================================================================================
INFO "Snapshotting the all-members proteins table before update"
proteins_map STRUCTURE_PAN_ALL.db > proteins_before_update.tsv
cat proteins_before_update.tsv

INFO "anvi-update-structure-database: ADD a second gene cluster to the pangenome structure db"
anvi-update-structure-database --pan-db TEST-PAN.db -g TEST-GENOMES.db \
                               -s STRUCTURE_PAN_ALL.db \
                               --gene-cluster-ids $GC_B \
                               --num-threads 1 --debug

# after adding GC_B we expect the original members PLUS the new cluster's members
assert_eq "$(proteins_count STRUCTURE_PAN_ALL.db)" "$((NUM_GENOMES * 2))" \
          "adding a second cluster appends its members"

INFO "Verifying the pre-existing protein_id <-> natural-key mapping was NOT disturbed"
proteins_map STRUCTURE_PAN_ALL.db | head -n $NUM_GENOMES > proteins_after_update_head.tsv
# the first NUM_GENOMES rows (lowest protein_ids) must be byte-identical to the pre-update snapshot
if ! diff -u proteins_before_update.tsv proteins_after_update_head.tsv; then
    echo ""
    echo "ASSERTION FAILED: update re-minted or reshuffled existing protein_ids"
    exit 1
fi
echo "  OK: existing protein_id <-> (genome, gene_callers_id, gene_cluster_id) mapping unchanged after add"

INFO "anvi-update-structure-database: RERUN the original cluster -- must NOT duplicate proteins"
proteins_map STRUCTURE_PAN_ALL.db > proteins_before_rerun.tsv
anvi-update-structure-database --pan-db TEST-PAN.db -g TEST-GENOMES.db \
                               -s STRUCTURE_PAN_ALL.db \
                               --gene-cluster-ids $GC_A \
                               --rerun \
                               --num-threads 1 --debug
proteins_map STRUCTURE_PAN_ALL.db > proteins_after_rerun.tsv

assert_eq "$(proteins_count STRUCTURE_PAN_ALL.db)" "$((NUM_GENOMES * 2))" \
          "rerun adds no new protein rows"
if ! diff -u proteins_before_rerun.tsv proteins_after_rerun.tsv; then
    echo ""
    echo "ASSERTION FAILED: rerun changed the proteins table"
    exit 1
fi
echo "  OK: rerun left the proteins table identical (idempotent)"

INFO "All pangenome structure-database assertions passed."

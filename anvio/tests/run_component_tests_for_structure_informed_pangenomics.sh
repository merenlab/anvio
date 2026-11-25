#!/bin/bash
# Component test for anvi-pan-genome-structure, which merges the gene clusters of a conventional pangenome
# into protein structure-informed gene clusters (PSGCs) based on the similarity of their structures.
#
# Real structure prediction needs MODELLER (and internet) or ColabFold (and a GPU). Instead, this test builds
# SYNTHETIC structures with `gen-synthetic-structures.py` (see its docstring): protein-like backbones for the
# exact amino acid sequences of the genes, whose geometry is set by a 'fold seed'. Two gene clusters with the
# same fold seed have near-identical structures and must end up in the same PSGC; gene clusters with different
# fold seeds must not. The synthetic structures enter a structure-db through the real pangenome path of
# anvi-gen-structure-database (--external-structures), so the whole chain is exercised:
#
#     anvi-pan-genome -> anvi-gen-structure-database --pan-db -> anvi-pan-genome-structure
#
# It needs foldseek and MCL in the PATH, but no internet or GPU.
#
# Run it from anvio/tests with no arguments (output goes to a fresh sandbox/test-output):
#     bash run_component_tests_for_structure_informed_pangenomics.sh

source 00.sh
set -e
set -o pipefail

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

expect_failure() {
    # expect_failure <expected text in the output> <message> <command ...>
    # runs the command, and asserts that it fails with an output that contains the expected text
    local expected="$1"
    local message="$2"
    shift 2
    if "$@" > expect_failure.log 2>&1; then
        echo ""
        echo "ASSERTION FAILED: $message"
        echo "  the command was expected to fail, but it succeeded: $*"
        exit 1
    fi
    if ! tr -s ' \n' '  ' < expect_failure.log | grep -q -- "$expected"; then
        echo ""
        echo "ASSERTION FAILED: $message"
        echo "  the command failed as expected, but its output does not contain '$expected':"
        cat expect_failure.log
        exit 1
    fi
    echo "  OK: $message (failed with '$expected')"
}

expect_in_log() {
    # expect_in_log <log file> <expected text> <message>
    if ! tr -s ' \n' '  ' < "$1" | grep -q -- "$2"; then
        echo ""
        echo "ASSERTION FAILED: $3"
        echo "  '$2' is not in '$1'"
        exit 1
    fi
    echo "  OK: $3"
}

file_md5()        { if command -v md5sum > /dev/null; then md5sum "$1" | cut -d' ' -f1; else md5 -q "$1"; fi; }
structure_of_gc() { awk -F'\t' -v gc="$2" '$1 == gc {print $2}' "$1"; }
num_psgcs()       { sqlite3 "$1" "SELECT COUNT(DISTINCT gene_cluster_id) FROM gene_clusters;"; }
psgc_of_gc()      { sqlite3 "$1" "SELECT protein_structure_informed_gene_cluster_id FROM gc_psgc_associations WHERE gene_cluster_id='$2';"; }
meta_value()      { sqlite3 "$1" "SELECT value FROM self WHERE key='$2';"; }
# the PSGC partition of gene clusters, for comparing two structure-informed pan-dbs
psgc_partition()  { sqlite3 -separator $'\t' "$1" "SELECT gene_cluster_id, protein_structure_informed_gene_cluster_id FROM gc_psgc_associations ORDER BY gene_cluster_id;"; }
# every gene of every PSGC, in the input gene cluster it comes from: if PSGCs are made of whole gene clusters,
# each PSGC gene must sit in a gene cluster that is associated with that very PSGC.
genes_outside_their_gcs() {
    sqlite3 "$1" "SELECT COUNT(*) FROM gene_clusters g \
                  JOIN gc_tracker t ON t.genome_name = g.genome_name AND t.gene_caller_id = g.gene_caller_id \
                  LEFT JOIN gc_psgc_associations a ON a.gene_cluster_id = t.gene_cluster_id \
                  WHERE a.protein_structure_informed_gene_cluster_id IS NULL \
                     OR a.protein_structure_informed_gene_cluster_id != g.gene_cluster_id;"
}
# number of genes in the input gene clusters that are associated with a PSGC
genes_in_associated_gcs() {
    sqlite3 "$1" "SELECT COUNT(*) FROM gc_tracker t JOIN gc_psgc_associations a ON a.gene_cluster_id = t.gene_cluster_id;"
}

# ---------------------------------------------------------------------------------------------------
# setup: a conventional pangenome with alignments
# ---------------------------------------------------------------------------------------------------
SETUP_WITH_OUTPUT_DIR $1 $2 $3

INFO "Setting up the structure-informed pangenome analysis directory"
mkdir -p $output_dir/
cp $files/mock_data_for_pangenomics/*.db                      $output_dir/
cp $files/mock_data_for_pangenomics/external-genomes.txt      $output_dir/
cp $files/mock_data_for_pangenomics/default-state.json        $output_dir/
cp $files/mock_data_for_pangenomics/group-information.txt     $output_dir/
cp $files/mock_data_for_pangenomics/gen-synthetic-structures.py $output_dir/
cp $files/example_description.md                              $output_dir/
cd $output_dir/

# output directories for the runs below
mkdir -p KEEP TXT FOLDSEEK STRICT UNKNOWN ALL RENAMED FOUR BAD

INFO "Migrating all databases"
anvi-migrate *db --migrate-quickly

INFO "Generating an anvi'o genomes storage"
anvi-gen-genomes-storage -e external-genomes.txt \
                         -o TEST-GENOMES.db \
                         --no-progress

# alignments are ON (the default), because representatives are picked from the aligned sequences
INFO "Running a conventional pangenome analysis"
anvi-pan-genome -g TEST-GENOMES.db \
                -n TEST \
                --no-progress \
                $thread_controller

NUM_GCS=$(sqlite3 TEST-PAN.db "SELECT COUNT(DISTINCT gene_cluster_id) FROM gene_clusters;")
INPUT_PAN_MD5=$(file_md5 TEST-PAN.db)

# five single-copy core gene clusters, longest first, among those with complete genes of at least 150 aa. the
# two longest share a fold seed (their lengths are close enough for a TM-score well above 0.5), the three
# others each get their own.
INFO "Picking gene clusters for synthetic structures"
sqlite3 -separator $'\t' TEST-PAN.db \
    "ATTACH 'TEST-GENOMES.db' AS g; \
     SELECT p.gene_cluster_id FROM gene_clusters p \
     JOIN g.gene_info i ON i.genome_name = p.genome_name AND i.gene_caller_id = p.gene_caller_id \
     GROUP BY p.gene_cluster_id \
     HAVING COUNT(*) = 5 AND COUNT(DISTINCT p.genome_name) = 5 AND SUM(i.partial) = 0 AND MIN(LENGTH(i.aa_sequence)) >= 150 \
     ORDER BY MIN(LENGTH(i.aa_sequence)) DESC, p.gene_cluster_id LIMIT 5;" > picked-gcs.txt
cat picked-gcs.txt

GC_PAIR_A=$(sed -n 1p picked-gcs.txt)
GC_PAIR_B=$(sed -n 2p picked-gcs.txt)
GC_SOLO_1=$(sed -n 3p picked-gcs.txt)
GC_SOLO_2=$(sed -n 4p picked-gcs.txt)
GC_SOLO_3=$(sed -n 5p picked-gcs.txt)
NUM_GCS_WITH_STRUCTURE=5

printf "$GC_PAIR_A\t1\n$GC_PAIR_B\t1\n$GC_SOLO_1\t2\n$GC_SOLO_2\t3\n$GC_SOLO_3\t4\n" > fold-seeds.txt

INFO "Building synthetic structures for the representatives of these gene clusters"
python gen-synthetic-structures.py TEST-PAN.db TEST-GENOMES.db fold-seeds.txt SYNTHETIC_REPS

INFO "Generating a structure-db from the pangenome, with the synthetic structures as external structures"
anvi-gen-structure-database --pan-db TEST-PAN.db \
                            -g TEST-GENOMES.db \
                            --external-structures SYNTHETIC_REPS/external-structures.txt \
                            --skip-DSSP \
                            -o TEST-STRUCTURE.db \
                            --no-progress

# ===================================================================================================
# the default run, from a structure-db
# ===================================================================================================
INFO "Running the structure-informed pangenome analysis from the structure-db"
anvi-pan-genome-structure -p TEST-PAN.db \
                          -g TEST-GENOMES.db \
                          --structure-db TEST-STRUCTURE.db \
                          --description example_description.md \
                          --no-progress \
                          $thread_controller 2>&1 | tee default-run.log

assert_eq "$(meta_value TEST-STRUCTURE-PAN.db db_variant)" "structure-informed" \
          "the output is a structure-informed pan-db, next to the input pan-db"
assert_eq "$(meta_value TEST-STRUCTURE-PAN.db structure_source)" "structure-db" \
          "the pan-db records where the structures came from"
assert_eq "$(sqlite3 TEST-STRUCTURE-PAN.db 'SELECT COUNT(DISTINCT gene_cluster_id) FROM gc_tracker;')" "$NUM_GCS" \
          "every gene cluster of the input pan-db is tracked"
assert_eq "$(sqlite3 TEST-STRUCTURE-PAN.db 'SELECT COUNT(*) FROM gc_psgc_associations;')" "$NUM_GCS_WITH_STRUCTURE" \
          "every gene cluster with a structure (and only those) is associated with a PSGC"
assert_eq "$(num_psgcs TEST-STRUCTURE-PAN.db)" "4" \
          "the two gene clusters with the same fold merge, the three others stay apart"
assert_eq "$(psgc_of_gc TEST-STRUCTURE-PAN.db $GC_PAIR_A)" "$(psgc_of_gc TEST-STRUCTURE-PAN.db $GC_PAIR_B)" \
          "$GC_PAIR_A and $GC_PAIR_B (same fold) are in the same PSGC"
assert_eq "$(sqlite3 TEST-STRUCTURE-PAN.db "SELECT COUNT(DISTINCT protein_structure_informed_gene_cluster_id) FROM gc_psgc_associations WHERE gene_cluster_id IN ('$GC_SOLO_1', '$GC_SOLO_2', '$GC_SOLO_3', '$GC_PAIR_A');")" "4" \
          "gene clusters with different folds are in different PSGCs"
assert_eq "$(genes_outside_their_gcs TEST-STRUCTURE-PAN.db)" "0" \
          "every gene of a PSGC comes from a gene cluster associated with that PSGC"
assert_eq "$(sqlite3 TEST-STRUCTURE-PAN.db 'SELECT COUNT(*) FROM gene_clusters;')" "$(genes_in_associated_gcs TEST-STRUCTURE-PAN.db)" \
          "PSGCs hold all the genes of their gene clusters, and nothing else"
expect_in_log default-run.log "do not have a structure" \
              "the gene clusters without a structure are reported"

# ===================================================================================================
# gene clusters without a structure can be kept
# ===================================================================================================
INFO "Keeping the gene clusters without a structure"
anvi-pan-genome-structure -p TEST-PAN.db \
                          -g TEST-GENOMES.db \
                          --structure-db TEST-STRUCTURE.db \
                          -o KEEP/TEST-KEEP-STRUCTURE-PAN.db \
                          --keep-gcs-without-structure \
                          --no-progress \
                          $thread_controller

assert_eq "$(sqlite3 KEEP/TEST-KEEP-STRUCTURE-PAN.db 'SELECT COUNT(*) FROM gc_psgc_associations;')" "$NUM_GCS" \
          "every gene cluster of the input pan-db is associated with a PSGC"
assert_eq "$(num_psgcs KEEP/TEST-KEEP-STRUCTURE-PAN.db)" "$((NUM_GCS - 1))" \
          "each gene cluster without a structure is its own PSGC"
assert_eq "$(sqlite3 KEEP/TEST-KEEP-STRUCTURE-PAN.db 'SELECT COUNT(*) FROM gene_clusters;')" "$(sqlite3 TEST-PAN.db 'SELECT COUNT(*) FROM gene_clusters;')" \
          "the structure-informed pangenome covers every gene of the input pangenome"
assert_eq "$(genes_outside_their_gcs KEEP/TEST-KEEP-STRUCTURE-PAN.db)" "0" \
          "every gene of a PSGC comes from a gene cluster associated with that PSGC"
assert_eq "$(sqlite3 KEEP/TEST-KEEP-STRUCTURE-PAN.db "SELECT COUNT(*) FROM item_additional_data WHERE data_key='has_structure' AND data_value='1';")" "4" \
          "only the PSGCs made from structures have has_structure=1"

# ===================================================================================================
# a structures-txt gives the same PSGCs as the structure-db it mirrors
# ===================================================================================================
INFO "Running the structure-informed pangenome analysis from a structures-txt"
anvi-pan-genome-structure -p TEST-PAN.db \
                          -g TEST-GENOMES.db \
                          --structures-txt SYNTHETIC_REPS/structures-txt.txt \
                          -o TXT/TEST-TXT-STRUCTURE-PAN.db \
                          --intermediate-data-dir TXT/INTERMEDIATE \
                          --no-progress \
                          $thread_controller

assert_eq "$(psgc_partition TXT/TEST-TXT-STRUCTURE-PAN.db)" "$(psgc_partition TEST-STRUCTURE-PAN.db)" \
          "the structures-txt and the structure-db give the same PSGCs"

INFO "Re-using the foldseek search results of the previous run"
anvi-pan-genome-structure -p TEST-PAN.db \
                          -g TEST-GENOMES.db \
                          --foldseek-search-results TXT/INTERMEDIATE/result \
                          -o FOLDSEEK/TEST-FOLDSEEK-STRUCTURE-PAN.db \
                          --keep-gcs-without-structure \
                          --no-progress \
                          $thread_controller

assert_eq "$(meta_value FOLDSEEK/TEST-FOLDSEEK-STRUCTURE-PAN.db structure_source)" "foldseek-search-results" \
          "the pan-db records where the structures came from"
assert_eq "$(psgc_partition FOLDSEEK/TEST-FOLDSEEK-STRUCTURE-PAN.db)" "$(psgc_partition KEEP/TEST-KEEP-STRUCTURE-PAN.db)" \
          "foldseek search results give the same PSGCs as the structures they came from"

INFO "A stricter TM-score threshold keeps the two similar gene clusters apart"
anvi-pan-genome-structure -p TEST-PAN.db \
                          -g TEST-GENOMES.db \
                          --structures-txt SYNTHETIC_REPS/structures-txt.txt \
                          -o STRICT/TEST-STRICT-STRUCTURE-PAN.db \
                          --min-tm-score 0.99 \
                          --no-progress \
                          $thread_controller

assert_eq "$(num_psgcs STRICT/TEST-STRICT-STRUCTURE-PAN.db)" "$NUM_GCS_WITH_STRUCTURE" \
          "with --min-tm-score 0.99, no two gene clusters merge"

INFO "A structures-txt entry that is not a gene cluster of the pan-db is ignored with a warning"
{ cat SYNTHETIC_REPS/structures-txt.txt; echo -e "GC_NOT_IN_THE_PAN\t$(sed -n 2p SYNTHETIC_REPS/structures-txt.txt | cut -f2)"; } > structures-txt-with-unknown-gc.txt
anvi-pan-genome-structure -p TEST-PAN.db \
                          -g TEST-GENOMES.db \
                          --structures-txt structures-txt-with-unknown-gc.txt \
                          -o UNKNOWN/TEST-UNKNOWN-STRUCTURE-PAN.db \
                          --no-progress \
                          $thread_controller > unknown-gc-run.log 2>&1
expect_in_log unknown-gc-run.log "do not correspond to any gene cluster in the input pan-db" \
              "the unknown gene cluster is reported"
assert_eq "$(psgc_partition UNKNOWN/TEST-UNKNOWN-STRUCTURE-PAN.db)" "$(psgc_partition TEST-STRUCTURE-PAN.db)" \
          "the unknown gene cluster does not change the PSGCs"

INFO "Malformed structures-txt files are rejected"
PDB_OF_PAIR_A=$(structure_of_gc SYNTHETIC_REPS/structures-txt.txt $GC_PAIR_A)
echo -e "gene_cluster\tpath\n$GC_PAIR_A\t$PDB_OF_PAIR_A" > bad-header.txt
expect_failure "The first column of" "a structures-txt with an unknown header" \
    anvi-pan-genome-structure -p TEST-PAN.db -g TEST-GENOMES.db --structures-txt bad-header.txt -o BAD/X-STRUCTURE-PAN.db
echo -e "gene_id\tpath\n$GC_PAIR_A\t$PDB_OF_PAIR_A\n$GC_PAIR_A\t$PDB_OF_PAIR_A" > duplicate-gene-id.txt
expect_failure "more than once" "a structures-txt with a duplicate gene_id" \
    anvi-pan-genome-structure -p TEST-PAN.db -g TEST-GENOMES.db --structures-txt duplicate-gene-id.txt -o BAD/X-STRUCTURE-PAN.db
echo -e "gene_id\tpath\n$GC_PAIR_A\tdoes/not/exist.pdb" > missing-file.txt
expect_failure "exist" "a structures-txt pointing at a missing file" \
    anvi-pan-genome-structure -p TEST-PAN.db -g TEST-GENOMES.db --structures-txt missing-file.txt -o BAD/X-STRUCTURE-PAN.db
cp $PDB_OF_PAIR_A structure-in-another-format.cif
echo -e "gene_id\tpath\n$GC_PAIR_A\tstructure-in-another-format.cif" > not-a-pdb.txt
expect_failure "accepts only .pdb files" "a structures-txt pointing at a non-.pdb structure" \
    anvi-pan-genome-structure -p TEST-PAN.db -g TEST-GENOMES.db --structures-txt not-a-pdb.txt -o BAD/X-STRUCTURE-PAN.db

# ===================================================================================================
# more than one protein with a structure for a gene cluster
# ===================================================================================================
INFO "Building synthetic structures for all members of $GC_PAIR_A"
printf "$GC_PAIR_A\t1\n" > fold-seeds-all-members.txt
python gen-synthetic-structures.py TEST-PAN.db TEST-GENOMES.db fold-seeds-all-members.txt SYNTHETIC_ALL --all-members

INFO "A structure-db with all members of a gene cluster: the representative of the gene cluster is used"
anvi-gen-structure-database --pan-db TEST-PAN.db -g TEST-GENOMES.db \
                            --external-structures SYNTHETIC_ALL/external-structures.txt \
                            --skip-DSSP -o TEST-STRUCTURE-ALL-MEMBERS.db --no-progress
anvi-pan-genome-structure -p TEST-PAN.db \
                          -g TEST-GENOMES.db \
                          --structure-db TEST-STRUCTURE-ALL-MEMBERS.db \
                          -o ALL/TEST-ALL-STRUCTURE-PAN.db \
                          --intermediate-data-dir ALL/INTERMEDIATE \
                          --no-progress \
                          $thread_controller > all-members-run.log 2>&1
expect_in_log all-members-run.log "used the structure of the representative" \
              "anvi'o tells which structure it used for the gene cluster"
REP_GENE_OF_PAIR_A=$(basename $(structure_of_gc SYNTHETIC_REPS/structures-txt.txt $GC_PAIR_A))
cmp ALL/INTERMEDIATE/structures-from-structure-db/$GC_PAIR_A.pdb SYNTHETIC_REPS/$REP_GENE_OF_PAIR_A
echo "  OK: the exported structure of $GC_PAIR_A is the one of its representative ($REP_GENE_OF_PAIR_A)"

INFO "A structure-db with several proteins for a gene cluster, none of which is its representative, is rejected"
grep -v "/$REP_GENE_OF_PAIR_A$" SYNTHETIC_ALL/external-structures.txt > external-structures-without-representative.txt
anvi-gen-structure-database --pan-db TEST-PAN.db -g TEST-GENOMES.db \
                            --external-structures external-structures-without-representative.txt \
                            --skip-DSSP -o TEST-STRUCTURE-NO-REP.db --no-progress
expect_failure "will not guess" "several proteins for a gene cluster, but not its representative" \
    anvi-pan-genome-structure -p TEST-PAN.db -g TEST-GENOMES.db --structure-db TEST-STRUCTURE-NO-REP.db -o BAD/X-STRUCTURE-PAN.db

# ===================================================================================================
# a structure-db and a pan-db that don't go together
# ===================================================================================================
INFO "A pan-db with the same genes in differently named gene clusters (as if it was computed again)"
{ echo -e "genome_name\tgene_caller_id\tgene_cluster_name"
  sqlite3 -separator $'\t' TEST-PAN.db "SELECT genome_name, gene_caller_id, 'RENAMED_' || gene_cluster_id FROM gene_clusters;"
} > renamed-gene-clusters.txt
anvi-pan-genome -g TEST-GENOMES.db \
                -n RENAMED \
                --gene-clusters-txt renamed-gene-clusters.txt \
                -o RENAMED/RENAMED-PAN.db \
                --no-progress \
                $thread_controller
expect_failure "do not belong to the gene cluster they were predicted for" "a structure-db made from another pan-db" \
    anvi-pan-genome-structure -p RENAMED/RENAMED-PAN.db -g TEST-GENOMES.db --structure-db TEST-STRUCTURE.db -o BAD/X-STRUCTURE-PAN.db

INFO "A pan-db made from another genomes storage"
head -n 5 external-genomes.txt > external-genomes-four.txt
anvi-gen-genomes-storage -e external-genomes-four.txt -o FOUR-GENOMES.db --no-progress
anvi-pan-genome -g FOUR-GENOMES.db -n FOUR -o FOUR/FOUR-PAN.db --no-progress $thread_controller
expect_failure "different genomes storage" "a structure-db made from a pangenome of another genomes storage" \
    anvi-pan-genome-structure -p FOUR/FOUR-PAN.db -g FOUR-GENOMES.db --structure-db TEST-STRUCTURE.db -o BAD/X-STRUCTURE-PAN.db

INFO "A structure-db that was not made from a pangenome"
cp TEST-STRUCTURE.db TEST-STRUCTURE-FROM-CONTIGS.db
sqlite3 TEST-STRUCTURE-FROM-CONTIGS.db "UPDATE self SET value='contigs_db' WHERE key='input_type';"
expect_failure "generated from a 'contigs_db' input" "a structure-db made from a contigs-db" \
    anvi-pan-genome-structure -p TEST-PAN.db -g TEST-GENOMES.db --structure-db TEST-STRUCTURE-FROM-CONTIGS.db -o BAD/X-STRUCTURE-PAN.db

# ===================================================================================================
# command line sanity
# ===================================================================================================
INFO "Command line sanity checks"
expect_failure "exactly one source of structures" "two structure sources at once" \
    anvi-pan-genome-structure -p TEST-PAN.db -g TEST-GENOMES.db --structure-db TEST-STRUCTURE.db --structures-txt SYNTHETIC_REPS/structures-txt.txt -o BAD/X-STRUCTURE-PAN.db
expect_failure "exactly one source of structures" "no structure source" \
    anvi-pan-genome-structure -p TEST-PAN.db -g TEST-GENOMES.db -o BAD/X-STRUCTURE-PAN.db
expect_failure "conventional, sequence-based pan-db" "a structure-informed pan-db as input" \
    anvi-pan-genome-structure -p TEST-STRUCTURE-PAN.db -g TEST-GENOMES.db --structure-db TEST-STRUCTURE.db -o BAD/X-STRUCTURE-PAN.db
expect_failure "must end with" "an output file name without the proper suffix" \
    anvi-pan-genome-structure -p TEST-PAN.db -g TEST-GENOMES.db --structure-db TEST-STRUCTURE.db -o BAD/X-PAN.db
expect_failure "already exists" "an existing output file" \
    anvi-pan-genome-structure -p TEST-PAN.db -g TEST-GENOMES.db --structure-db TEST-STRUCTURE.db

assert_eq "$(file_md5 TEST-PAN.db)" "$INPUT_PAN_MD5" \
          "none of the runs above modified the input pan-db"

# ===================================================================================================
# comparison with the input pan-db (opt-in)
# ===================================================================================================
INFO "Comparing the structure-informed pan-db with the input pan-db"
mkdir -p COMPARE
cp TEST-PAN.db COMPARE/TEST-PAN.db
anvi-pan-genome-structure -p COMPARE/TEST-PAN.db \
                          -g TEST-GENOMES.db \
                          --structure-db TEST-STRUCTURE.db \
                          --compare-with-input-pan \
                          --no-progress \
                          $thread_controller

assert_eq "$(sqlite3 COMPARE/TEST-STRUCTURE-PAN.db "SELECT COUNT(DISTINCT data_key) > 0 FROM item_additional_data WHERE data_key LIKE 'sequence_based_%';")" "1" \
          "the structure-informed pan-db has comparison layers"
assert_eq "$(file_md5 COMPARE/TEST-PAN.db)" "$INPUT_PAN_MD5" \
          "the comparison does not modify the input pan-db"
SHOW_FILE COMPARE/TEST-STRUCTURE-VS-SEQUENCE.txt

# ===================================================================================================
# downstream programs work with a structure-informed pan-db
# ===================================================================================================
INFO "Importing group information as misc data for layers"
anvi-import-misc-data -p KEEP/TEST-KEEP-STRUCTURE-PAN.db \
                      -t layers \
                      group-information.txt \
                      --no-progress

INFO "Estimating enriched functions per pan group"
anvi-compute-functional-enrichment-in-pan -p KEEP/TEST-KEEP-STRUCTURE-PAN.db \
                                          -g TEST-GENOMES.db \
                                          --category group \
                                          --annotation-source COG20_FUNCTION \
                                          -o functions-enrichment.txt \
                                          -F functional-occurence.txt \
                                          --no-progress
SHOW_FILE functions-enrichment.txt
SHOW_FILE functional-occurence.txt

INFO "Importing a collection of PSGCs (the ones made from structures, and the others)"
sqlite3 -separator $'\t' KEEP/TEST-KEEP-STRUCTURE-PAN.db \
    "SELECT item_name, CASE WHEN data_value = '1' THEN 'FROM_STRUCTURES' ELSE 'WITHOUT_STRUCTURE' END \
     FROM item_additional_data WHERE data_key = 'has_structure';" > psgcs-collection.txt
anvi-import-collection -p KEEP/TEST-KEEP-STRUCTURE-PAN.db \
                       -C test_collection psgcs-collection.txt \
                       --no-progress

INFO "Summarizing the pan, using the test collection (in quick mode)"
anvi-summarize -p KEEP/TEST-KEEP-STRUCTURE-PAN.db \
               -g TEST-GENOMES.db \
               -C test_collection \
               -o TEST_SUMMARY_QUICK \
               --init-pan-mode structure-informed \
               --quick \
               --no-progress

INFO "Summarizing the pan, using the test collection"
anvi-summarize -p KEEP/TEST-KEEP-STRUCTURE-PAN.db \
               -g TEST-GENOMES.db \
               -C test_collection \
               -o TEST_SUMMARY \
               --init-pan-mode structure-informed \
               --no-progress

INFO "Splitting bins in the pan genome into smaller, self-contained pan databases"
anvi-split -p KEEP/TEST-KEEP-STRUCTURE-PAN.db \
           -g TEST-GENOMES.db \
           -C test_collection \
           -o TEST_SPLIT_PAN

INFO "Resulting split pans"
ls -l TEST_SPLIT_PAN/*-PAN.db

INFO "Taking a look at the make up one of the split pans"
anvi-db-info TEST_SPLIT_PAN/FROM_STRUCTURES-PAN.db
assert_eq "$(meta_value TEST_SPLIT_PAN/FROM_STRUCTURES-PAN.db db_variant)" "structure-informed" \
          "a split structure-informed pan-db is structure-informed, too"
assert_eq "$(sqlite3 TEST_SPLIT_PAN/FROM_STRUCTURES-PAN.db 'SELECT COUNT(*) FROM gc_psgc_associations;')" "$NUM_GCS_WITH_STRUCTURE" \
          "the split pan-db keeps the gene clusters of its PSGCs"

INFO "Listing collections available"
anvi-show-collections-and-bins -p KEEP/TEST-KEEP-STRUCTURE-PAN.db \
                               --no-progress

INFO "Importing the default state for pretty outputs"
anvi-import-state -p KEEP/TEST-KEEP-STRUCTURE-PAN.db -s default-state.json -n default

INFO "Displaying the structure-informed pangenome"
anvi-display-pan -p KEEP/TEST-KEEP-STRUCTURE-PAN.db \
                 -g TEST-GENOMES.db \
                 --title "A mock structure-informed pangenome analysis" \
                 --no-progress \
                 $dry_run_controller

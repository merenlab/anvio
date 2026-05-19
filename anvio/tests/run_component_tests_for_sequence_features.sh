#!/bin/bash
source 00.sh

# Setup #############################
SETUP_WITH_OUTPUT_DIR $1 $2 $3
#####################################

# Helper: assert a sqlite3 query returns an expected value. Aborts the test
# with set -e if the actual output does not match.
assert_query() {
    local db="$1"
    local query="$2"
    local expected="$3"
    local label="$4"
    local actual
    actual=$(sqlite3 "$db" "$query")
    if [ "$actual" != "$expected" ]; then
        echo "ASSERT FAILED [$label]: expected '$expected', got '$actual' for query: $query"
        exit 1
    else
        echo "  OK [$label]: $actual"
    fi
}

NEW_TABLES="contigs_sequence_features feature_types feature_relationships feature_qualifiers CDS_features"

####################################################################################################
# Test 1: v24 -> v25 migration
####################################################################################################
INFO "Test 1: migrating a v24 contigs database forward to v25"

cp $files/sequence_features/CONTIGS-v24.db $output_dir/MIGRATED.db

# capture pre-migration row counts for invariant tables we expect to remain untouched
PRE_CONTIGS=$(sqlite3 $output_dir/MIGRATED.db "SELECT COUNT(*) FROM contigs_basic_info")
PRE_SPLITS=$(sqlite3 $output_dir/MIGRATED.db "SELECT COUNT(*) FROM splits_basic_info")
PRE_GENES=$(sqlite3 $output_dir/MIGRATED.db "SELECT COUNT(*) FROM genes_in_contigs")

assert_query $output_dir/MIGRATED.db "SELECT value FROM self WHERE key='version'" "24" "pre-migration version is 24"

anvi-migrate $output_dir/MIGRATED.db --migrate-quickly

assert_query $output_dir/MIGRATED.db "SELECT value FROM self WHERE key='version'" "25" "post-migration version is 25"

for tbl in $NEW_TABLES; do
    assert_query $output_dir/MIGRATED.db "SELECT COUNT(*) FROM sqlite_master WHERE type='table' AND name='$tbl'" "1" "$tbl table exists after migration"
done

assert_query $output_dir/MIGRATED.db "SELECT COUNT(*) FROM feature_types" "5" "feature_types has 5 builtin rows"
assert_query $output_dir/MIGRATED.db "SELECT COUNT(*) FROM feature_types WHERE is_builtin=1" "5" "all 5 feature_types are builtin"
assert_query $output_dir/MIGRATED.db "SELECT GROUP_CONCAT(feature_type, ',') FROM (SELECT feature_type FROM feature_types ORDER BY feature_type)" "CDS,exon,gene,intron,mRNA" "feature_types contains the five expected names"
assert_query $output_dir/MIGRATED.db "SELECT has_dedicated_table FROM feature_types WHERE feature_type='CDS'" "1" "CDS feature type is flagged has_dedicated_table"
assert_query $output_dir/MIGRATED.db "SELECT has_dedicated_table FROM feature_types WHERE feature_type='gene'" "0" "gene feature type is not flagged has_dedicated_table"

for tbl in contigs_sequence_features feature_relationships feature_qualifiers CDS_features; do
    assert_query $output_dir/MIGRATED.db "SELECT COUNT(*) FROM $tbl" "0" "$tbl is empty after migration"
done

# invariant: untouched tables retain their pre-migration row counts
assert_query $output_dir/MIGRATED.db "SELECT COUNT(*) FROM contigs_basic_info" "$PRE_CONTIGS" "contigs_basic_info row count is preserved"
assert_query $output_dir/MIGRATED.db "SELECT COUNT(*) FROM splits_basic_info"  "$PRE_SPLITS"  "splits_basic_info row count is preserved"
assert_query $output_dir/MIGRATED.db "SELECT COUNT(*) FROM genes_in_contigs"   "$PRE_GENES"   "genes_in_contigs row count is preserved"

# verify the indexes anvi'o expects on the new tables actually exist
for idx in idx_pk_contigs_sequence_features idx_pk_feature_types idx_pk_feature_relationships idx_pk_feature_qualifiers idx_pk_CDS_features idx_contigs_sequence_features_contig_start idx_contigs_sequence_features_contig_stop idx_contigs_sequence_features_feature_type idx_contigs_sequence_features_group_id idx_contigs_sequence_features_gcid idx_feature_relationships_parent; do
    assert_query $output_dir/MIGRATED.db "SELECT COUNT(*) FROM sqlite_master WHERE type='index' AND name='$idx'" "1" "$idx exists after migration"
done

####################################################################################################
# Test 2: fresh v25 contigs database initialization
####################################################################################################
INFO "Test 2: anvi-gen-contigs-database produces a v25 db with the new tables empty"

anvi-gen-contigs-database -f $files/contigs.fa \
                          -o $output_dir/FRESH.db \
                          --project-name "SequenceFeaturesFreshInit" \
                          --skip-gene-calling \
                          -L 1000 \
                          --no-progress \
                          $thread_controller

assert_query $output_dir/FRESH.db "SELECT value FROM self WHERE key='version'" "25" "fresh db is at version 25"

for tbl in $NEW_TABLES; do
    assert_query $output_dir/FRESH.db "SELECT COUNT(*) FROM sqlite_master WHERE type='table' AND name='$tbl'" "1" "$tbl table exists in fresh db"
done

assert_query $output_dir/FRESH.db "SELECT COUNT(*) FROM feature_types" "5" "feature_types has 5 builtin rows in fresh db"
for tbl in contigs_sequence_features feature_relationships feature_qualifiers CDS_features; do
    assert_query $output_dir/FRESH.db "SELECT COUNT(*) FROM $tbl" "0" "$tbl is empty in fresh db"
done

# pre-existing tables behave normally
assert_query $output_dir/FRESH.db "SELECT COUNT(*) FROM contigs_basic_info" "6" "fresh db has 6 contigs (matching contigs.fa)"

# fresh-init end state and migration end state must agree on schema-defining
# objects (CREATE TABLE / CREATE INDEX SQL). Compare via sqlite_master.
INFO "Comparing the schema of the freshly-initialised db against the migrated db"
FRESH_SCHEMA=$(sqlite3 $output_dir/FRESH.db    "SELECT sql FROM sqlite_master WHERE tbl_name IN ('contigs_sequence_features','feature_types','feature_relationships','feature_qualifiers','CDS_features') AND sql IS NOT NULL ORDER BY name")
MIG_SCHEMA=$(sqlite3   $output_dir/MIGRATED.db "SELECT sql FROM sqlite_master WHERE tbl_name IN ('contigs_sequence_features','feature_types','feature_relationships','feature_qualifiers','CDS_features') AND sql IS NOT NULL ORDER BY name")
if [ "$FRESH_SCHEMA" != "$MIG_SCHEMA" ]; then
    echo "ASSERT FAILED: schemas of freshly-initialised and migrated databases differ:"
    diff <(echo "$FRESH_SCHEMA") <(echo "$MIG_SCHEMA")
    exit 1
fi
echo "  OK [schema parity]: fresh-init and migration produce identical sqlite_master entries"

####################################################################################################
echo
echo "Sequence-features component tests 1 & 2 PASSED."

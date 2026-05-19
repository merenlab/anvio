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
# Test 2.5: ContigsSuperclass read accessor probe (manually-populated db)
####################################################################################################
INFO "Test 2.5: probing ContigsSuperclass sequence-features read accessors"

cp $output_dir/FRESH.db $output_dir/PROBE.db
python $files/sequence_features/probe_read_accessors.py $output_dir/PROBE.db

####################################################################################################
# Test 2.6: TablesForSequenceFeatures + db.transaction() write-side probe
####################################################################################################
INFO "Test 2.6: probing TablesForSequenceFeatures.populate_features and the transaction context manager"

cp $output_dir/FRESH.db $output_dir/POPULATE.db
python $files/sequence_features/probe_populate_features.py $output_dir/POPULATE.db

####################################################################################################
# Test 2.7: GenbankFeatureImporter parsing-correctness probe
####################################################################################################
INFO "Test 2.7: probing GenbankFeatureImporter (bacterial, eukaryotic multi-seg, origin-crossing, malformed-linear, unknown type)"

mkdir -p $output_dir/importer_probe
python $files/sequence_features/probe_importer.py $output_dir/importer_probe

####################################################################################################
# Test 3: bacterial GenBank end-to-end (anvi-script-process-genbank + anvi-import-genbank-features)
####################################################################################################
INFO "Test 3: bacterial GenBank end-to-end through the CLI"

anvi-script-process-genbank -i $files/sequence_features/bacterial.gb \
                            --output-fasta $output_dir/bact_t3.fa \
                            --output-gene-calls $output_dir/bact_t3.tsv \
                            --output-functions $output_dir/bact_t3_funcs.tsv \
                            --annotation-source 'genbank_test_3' >/dev/null

anvi-gen-contigs-database -f $output_dir/bact_t3.fa \
                          -o $output_dir/BACT_T3.db \
                          --external-gene-calls $output_dir/bact_t3.tsv \
                          --ignore-internal-stop-codons \
                          --project-name "BacterialTest3" \
                          -L 1000 --no-progress $thread_controller >/dev/null

anvi-import-genbank-features -c $output_dir/BACT_T3.db \
                             -i $files/sequence_features/bacterial.gb \
                             --source-name 'bact_t3' >/dev/null

assert_query $output_dir/BACT_T3.db "SELECT COUNT(*) FROM contigs_sequence_features WHERE feature_type='gene'" "4" "test 3: 4 gene rows"
assert_query $output_dir/BACT_T3.db "SELECT COUNT(*) FROM contigs_sequence_features WHERE feature_type='CDS'" "4" "test 3: 4 CDS rows"
assert_query $output_dir/BACT_T3.db "SELECT COUNT(*) FROM CDS_features WHERE translation IS NOT NULL" "4" "test 3: 4 CDS translations stored"
assert_query $output_dir/BACT_T3.db "SELECT COUNT(*) FROM feature_relationships WHERE relationship='derives_from'" "4" "test 3: 4 derives_from relationships"
# external gene calls match GenBank coords exactly → all 4 genes should have GCIDs
assert_query $output_dir/BACT_T3.db "SELECT COUNT(*) FROM contigs_sequence_features WHERE feature_type='gene' AND gene_callers_id IS NOT NULL" "4" "test 3: all 4 gene rows have gene_callers_id reconciled from external gene calls"
# zero-based half-open: TST_001's gene is at GenBank 100..500 → start=99, stop=500
assert_query $output_dir/BACT_T3.db "SELECT start FROM contigs_sequence_features WHERE external_id='TST_001' AND feature_type='gene'" "99" "test 3: zero-based start for TST_001"
assert_query $output_dir/BACT_T3.db "SELECT stop FROM contigs_sequence_features WHERE external_id='TST_001' AND feature_type='gene'" "500" "test 3: half-open stop for TST_001"
assert_query $output_dir/BACT_T3.db "SELECT direction FROM contigs_sequence_features WHERE external_id='TST_002' AND feature_type='gene'" "r" "test 3: reverse direction encoded as 'r'"
# product qualifier preserved on CDS
assert_query $output_dir/BACT_T3.db "SELECT value FROM feature_qualifiers WHERE feature_id IN (SELECT feature_id FROM contigs_sequence_features WHERE external_id='TST_001' AND feature_type='CDS') AND key='product'" "Test protein A" "test 3: product qualifier preserved on CDS"

####################################################################################################
# Test 4: eukaryotic multi-segment CDS end-to-end
####################################################################################################
INFO "Test 4: eukaryotic multi-segment CDS end-to-end"

anvi-gen-contigs-database -f $files/sequence_features/eukaryotic.fa \
                          -o $output_dir/EUK_T4.db \
                          --skip-gene-calling \
                          --project-name "EukaryoticTest4" \
                          -L 1000 --no-progress $thread_controller >/dev/null

anvi-import-genbank-features -c $output_dir/EUK_T4.db \
                             -i $files/sequence_features/eukaryotic.gb \
                             --source-name 'euk_t4' >/dev/null

# forward CDS: 3 segments
assert_query $output_dir/EUK_T4.db "SELECT COUNT(*) FROM contigs_sequence_features WHERE feature_type='CDS' AND direction='f'" "3" "test 4: forward CDS has 3 segments"
# all three share the same feature_group_id
assert_query $output_dir/EUK_T4.db "SELECT COUNT(DISTINCT feature_group_id) FROM contigs_sequence_features WHERE feature_type='CDS' AND direction='f'" "1" "test 4: forward CDS segments share a feature_group_id"
# segment_order 0 has the leftmost start
assert_query $output_dir/EUK_T4.db "SELECT start FROM contigs_sequence_features WHERE feature_type='CDS' AND direction='f' AND segment_order=0" "99" "test 4: forward CDS segment_order=0 is leftmost"
# reverse CDS: segment_order=0 must be GENOMIC RIGHTMOST (start=3999)
assert_query $output_dir/EUK_T4.db "SELECT start FROM contigs_sequence_features WHERE feature_type='CDS' AND direction='r' AND segment_order=0" "3999" "test 4: reverse CDS segment_order=0 is genomic rightmost"
# canonical group traversal returns segments in segment_order
assert_query $output_dir/EUK_T4.db "SELECT GROUP_CONCAT(segment_order, ',') FROM (SELECT segment_order FROM contigs_sequence_features WHERE COALESCE(feature_group_id, feature_id) = (SELECT feature_id FROM contigs_sequence_features WHERE feature_type='CDS' AND direction='f' AND segment_order=0) ORDER BY segment_order)" "0,1,2" "test 4: canonical-group query returns segments in segment_order"
# translation is stored only on segment_order=0
assert_query $output_dir/EUK_T4.db "SELECT COUNT(*) FROM CDS_features cf JOIN contigs_sequence_features csf USING(feature_id) WHERE cf.translation IS NOT NULL AND csf.direction='f' AND csf.feature_type='CDS'" "1" "test 4: forward CDS translation only on segment_order=0"
# canonical-parent convention: 3-segment CDS → 3 rows in feature_relationships pointing to the canonical mRNA
assert_query $output_dir/EUK_T4.db "SELECT COUNT(*) FROM feature_relationships fr JOIN contigs_sequence_features csf ON fr.child_feature_id = csf.feature_id WHERE csf.feature_type='CDS' AND csf.direction='f'" "3" "test 4: 3 CDS→mRNA relationship rows (one per CDS segment)"
assert_query $output_dir/EUK_T4.db "SELECT COUNT(DISTINCT fr.parent_feature_id) FROM feature_relationships fr JOIN contigs_sequence_features csf ON fr.child_feature_id = csf.feature_id WHERE csf.feature_type='CDS' AND csf.direction='f'" "1" "test 4: all forward CDS segments share one canonical mRNA parent"
# external_id replicates across every segment of a multi-segment feature
assert_query $output_dir/EUK_T4.db "SELECT COUNT(DISTINCT external_id) FROM contigs_sequence_features WHERE feature_type='CDS' AND direction='f'" "1" "test 4: external_id replicates across all forward CDS segments"

####################################################################################################
# Test 5: origin-crossing feature on a circular contig
####################################################################################################
INFO "Test 5: origin-crossing feature on a circular contig"

anvi-gen-contigs-database -f $files/sequence_features/circular.fa \
                          -o $output_dir/CIRC_T5.db \
                          --skip-gene-calling \
                          --project-name "CircularTest5" \
                          -L 1000 --no-progress $thread_controller >/dev/null

anvi-import-genbank-features -c $output_dir/CIRC_T5.db \
                             -i $files/sequence_features/circular.gb \
                             --source-name 'circ_t5' >/dev/null

assert_query $output_dir/CIRC_T5.db "SELECT COUNT(*) FROM contigs_sequence_features" "2" "test 5: origin-crossing stored as exactly 2 rows"
assert_query $output_dir/CIRC_T5.db "SELECT start FROM contigs_sequence_features WHERE segment_order=0" "899" "test 5: segment_order=0 is [899,1000)"
assert_query $output_dir/CIRC_T5.db "SELECT stop  FROM contigs_sequence_features WHERE segment_order=0" "1000" "test 5: segment_order=0 stop is 1000"
assert_query $output_dir/CIRC_T5.db "SELECT start FROM contigs_sequence_features WHERE segment_order=1" "0" "test 5: segment_order=1 is [0,200)"
assert_query $output_dir/CIRC_T5.db "SELECT stop  FROM contigs_sequence_features WHERE segment_order=1" "200" "test 5: segment_order=1 stop is 200"
assert_query $output_dir/CIRC_T5.db "SELECT COUNT(DISTINCT feature_group_id) FROM contigs_sequence_features" "1" "test 5: both segments share feature_group_id"
assert_query $output_dir/CIRC_T5.db "SELECT COUNT(DISTINCT external_id) FROM contigs_sequence_features" "1" "test 5: both segments share external_id (locus_tag)"

####################################################################################################
# Test 6: malformed origin-crossing on a LINEAR contig
####################################################################################################
INFO "Test 6: malformed origin-crossing on a linear contig is warn-and-skip"

anvi-gen-contigs-database -f $files/sequence_features/linear_malformed.fa \
                          -o $output_dir/LIN_T6.db \
                          --skip-gene-calling \
                          --project-name "LinearMalformedTest6" \
                          -L 1000 --no-progress $thread_controller >/dev/null

anvi-import-genbank-features -c $output_dir/LIN_T6.db \
                             -i $files/sequence_features/linear_malformed.gb \
                             --source-name 'lin_t6' 2> $output_dir/lin_t6_stderr.txt > /dev/null

# malformed feature absent; healthy feature present
assert_query $output_dir/LIN_T6.db "SELECT COUNT(*) FROM contigs_sequence_features WHERE external_id='LINEAR_BAD_001'" "0" "test 6: malformed origin-crosser absent"
assert_query $output_dir/LIN_T6.db "SELECT COUNT(*) FROM contigs_sequence_features WHERE external_id='LINEAR_OK_001'" "1" "test 6: healthy feature still imported"
grep -q "origin-crossing" $output_dir/lin_t6_stderr.txt || \
    grep -q "Skipping a feature" $output_dir/lin_t6_stderr.txt || \
    { echo "ASSERT FAILED [test 6: warning was emitted on stderr]"; exit 1; }
echo "  OK [test 6: warning emitted on stderr]"

####################################################################################################
# Test 7: --force and multi-source coexistence (also verifies determinism)
####################################################################################################
INFO "Test 7: --force, multi-source coexistence, and feature_id determinism"

anvi-gen-contigs-database -f $files/sequence_features/bacterial.fa \
                          -o $output_dir/BACT_T7.db \
                          --skip-gene-calling \
                          --project-name "ForceMultiSourceTest7" \
                          -L 1000 --no-progress $thread_controller >/dev/null

# first import: source alpha
anvi-import-genbank-features -c $output_dir/BACT_T7.db \
                             -i $files/sequence_features/bacterial.gb \
                             --source-name 'alpha' >/dev/null
N_ALPHA_BEFORE=$(sqlite3 $output_dir/BACT_T7.db "SELECT COUNT(*) FROM contigs_sequence_features WHERE source='alpha'")
IDS_ALPHA_BEFORE=$(sqlite3 $output_dir/BACT_T7.db "SELECT GROUP_CONCAT(feature_id, ',') FROM (SELECT feature_id FROM contigs_sequence_features WHERE source='alpha' ORDER BY feature_id)")

# re-import alpha without --force should fail
if anvi-import-genbank-features -c $output_dir/BACT_T7.db -i $files/sequence_features/bacterial.gb --source-name 'alpha' >/dev/null 2>&1; then
    echo "ASSERT FAILED [test 7: re-import without --force unexpectedly succeeded]"; exit 1
fi
echo "  OK [test 7: re-import without --force fails as expected]"

# re-import alpha with --force should succeed and the set of feature_ids should be IDENTICAL (determinism)
anvi-import-genbank-features -c $output_dir/BACT_T7.db -i $files/sequence_features/bacterial.gb --source-name 'alpha' --force >/dev/null
N_ALPHA_AFTER=$(sqlite3 $output_dir/BACT_T7.db "SELECT COUNT(*) FROM contigs_sequence_features WHERE source='alpha'")
IDS_ALPHA_AFTER=$(sqlite3 $output_dir/BACT_T7.db "SELECT GROUP_CONCAT(feature_id, ',') FROM (SELECT feature_id FROM contigs_sequence_features WHERE source='alpha' ORDER BY feature_id)")
[ "$N_ALPHA_BEFORE" = "$N_ALPHA_AFTER" ] || { echo "ASSERT FAILED [test 7: alpha row count changed after --force]"; exit 1; }
[ "$IDS_ALPHA_BEFORE" = "$IDS_ALPHA_AFTER" ] || { echo "ASSERT FAILED [test 7: feature_ids differ after --force; not deterministic]"; exit 1; }
echo "  OK [test 7: --force replays produce identical feature_id set]"

# now add a second source (no --force needed) and verify both sources coexist
anvi-import-genbank-features -c $output_dir/BACT_T7.db -i $files/sequence_features/bacterial.gb --source-name 'beta' >/dev/null
assert_query $output_dir/BACT_T7.db "SELECT COUNT(*) FROM contigs_sequence_features WHERE source='alpha'" "$N_ALPHA_BEFORE" "test 7: alpha source survives addition of beta"
assert_query $output_dir/BACT_T7.db "SELECT COUNT(*) FROM contigs_sequence_features WHERE source='beta'"  "$N_ALPHA_BEFORE" "test 7: beta source imported with the same row count"
assert_query $output_dir/BACT_T7.db "SELECT COUNT(DISTINCT source) FROM contigs_sequence_features" "2" "test 7: two distinct sources coexist in contigs_sequence_features"

# --force on alpha must not affect beta
anvi-import-genbank-features -c $output_dir/BACT_T7.db -i $files/sequence_features/bacterial.gb --source-name 'alpha' --force >/dev/null
assert_query $output_dir/BACT_T7.db "SELECT COUNT(*) FROM contigs_sequence_features WHERE source='beta'" "$N_ALPHA_BEFORE" "test 7: --force on alpha leaves beta untouched"

####################################################################################################
# Test 8: unknown feature type registers as non-builtin
####################################################################################################
INFO "Test 8: unknown feature type registers as non-builtin"

anvi-gen-contigs-database -f $files/sequence_features/unknown_type.fa \
                          -o $output_dir/UNK_T8.db \
                          --skip-gene-calling \
                          --project-name "UnknownTypeTest8" \
                          -L 1000 --no-progress $thread_controller >/dev/null

anvi-import-genbank-features -c $output_dir/UNK_T8.db \
                             -i $files/sequence_features/unknown_type.gb \
                             --source-name 'unk_t8' >/dev/null

assert_query $output_dir/UNK_T8.db "SELECT COUNT(*) FROM feature_types WHERE feature_type='regulatory'" "1" "test 8: regulatory registered in feature_types"
assert_query $output_dir/UNK_T8.db "SELECT is_builtin FROM feature_types WHERE feature_type='regulatory'" "0" "test 8: regulatory has is_builtin=0"
assert_query $output_dir/UNK_T8.db "SELECT has_dedicated_table FROM feature_types WHERE feature_type='regulatory'" "0" "test 8: regulatory has has_dedicated_table=0"
assert_query $output_dir/UNK_T8.db "SELECT COUNT(*) FROM contigs_sequence_features WHERE feature_type='regulatory'" "1" "test 8: regulatory feature imported"

####################################################################################################
# Test 9: GCID mismatch when gene calls don't align with GenBank
####################################################################################################
INFO "Test 9: GCID mismatch warning when gene calls diverge from GenBank"

# Use pyrodigal_gv (anvi'o's default) on the bacterial fasta — its predictions will
# differ from the GenBank-derived coordinates, so most genes won't reconcile.
anvi-gen-contigs-database -f $files/sequence_features/bacterial.fa \
                          -o $output_dir/BACT_T9.db \
                          --project-name "GCIDMismatchTest9" \
                          -L 1000 --no-progress $thread_controller >/dev/null 2>&1 || \
    { echo "(pyrodigal_gv may not be available; skipping test 9)"; SKIP_T9=1; }

if [ -z "$SKIP_T9" ]; then
    anvi-import-genbank-features -c $output_dir/BACT_T9.db \
                                 -i $files/sequence_features/bacterial.gb \
                                 --source-name 'gcid_test' 2> $output_dir/gcid_stderr.txt > /dev/null

    # all features should still be imported regardless of GCID match outcome
    assert_query $output_dir/BACT_T9.db "SELECT COUNT(*) FROM contigs_sequence_features" "8" "test 9: all 8 features imported despite gene-call divergence"
    # at least one warning should be emitted on stderr
    grep -q "No gene call matches" $output_dir/gcid_stderr.txt || { echo "ASSERT FAILED [test 9: expected at least one GCID mismatch warning on stderr]"; exit 1; }
    echo "  OK [test 9: GCID mismatch warning emitted on stderr]"
fi

####################################################################################################
# Test 10: error paths produce clean ConfigError messages (not stack traces)
####################################################################################################
INFO "Test 10: error paths produce ConfigError"

# 10a: nonexistent contigs db path
if anvi-import-genbank-features -c /no/such/path.db -i $files/sequence_features/bacterial.gb 2>$output_dir/t10a.txt >/dev/null; then
    echo "ASSERT FAILED [test 10a: nonexistent contigs DB should fail]"; exit 1
fi
grep -q "Traceback" $output_dir/t10a.txt && { echo "ASSERT FAILED [test 10a: stack trace leaked to stderr]"; exit 1; }
echo "  OK [test 10a: nonexistent contigs DB fails cleanly]"

# 10b: too-old contigs DB (the un-migrated v24 fixture)
cp $files/sequence_features/CONTIGS-v24.db $output_dir/V24_RAW.db
if anvi-import-genbank-features -c $output_dir/V24_RAW.db -i $files/sequence_features/bacterial.gb 2>$output_dir/t10b.txt >/dev/null; then
    echo "ASSERT FAILED [test 10b: v24 db should fail without migration]"; exit 1
fi
grep -q "Traceback" $output_dir/t10b.txt && { echo "ASSERT FAILED [test 10b: stack trace leaked]"; exit 1; }
echo "  OK [test 10b: too-old contigs DB fails cleanly]"

# 10c: invalid --source-name
if anvi-import-genbank-features -c $output_dir/BACT_T3.db -i $files/sequence_features/bacterial.gb --source-name 'bad name!' 2>$output_dir/t10c.txt >/dev/null; then
    echo "ASSERT FAILED [test 10c: invalid source-name should fail]"; exit 1
fi
grep -q "Traceback" $output_dir/t10c.txt && { echo "ASSERT FAILED [test 10c: stack trace leaked]"; exit 1; }
echo "  OK [test 10c: invalid --source-name fails cleanly]"

# 10d: GenBank LOCUS that doesn't match any contig in the database
anvi-gen-contigs-database -f $files/sequence_features/eukaryotic.fa \
                          -o $output_dir/MISMATCH.db \
                          --skip-gene-calling \
                          --project-name "MismatchTest10d" \
                          -L 1000 --no-progress $thread_controller >/dev/null
if anvi-import-genbank-features -c $output_dir/MISMATCH.db -i $files/sequence_features/bacterial.gb 2>$output_dir/t10d.txt >/dev/null; then
    echo "ASSERT FAILED [test 10d: unmatched LOCUS should fail]"; exit 1
fi
grep -q "Traceback" $output_dir/t10d.txt && { echo "ASSERT FAILED [test 10d: stack trace leaked]"; exit 1; }
echo "  OK [test 10d: unmatched LOCUS fails cleanly]"

# 10e: re-import without --force when source already exists (also covered by test 7 but we re-check here for clarity)
if anvi-import-genbank-features -c $output_dir/BACT_T3.db -i $files/sequence_features/bacterial.gb --source-name 'bact_t3' 2>$output_dir/t10e.txt >/dev/null; then
    echo "ASSERT FAILED [test 10e: re-import without --force should fail]"; exit 1
fi
grep -q "Traceback" $output_dir/t10e.txt && { echo "ASSERT FAILED [test 10e: stack trace leaked]"; exit 1; }
echo "  OK [test 10e: re-import without --force fails cleanly]"

# 10f: a different source name on the same DB succeeds (the inverse of 10e — multi-source coexistence)
anvi-import-genbank-features -c $output_dir/BACT_T3.db -i $files/sequence_features/bacterial.gb --source-name 'bact_t3_other' >/dev/null
assert_query $output_dir/BACT_T3.db "SELECT COUNT(*) FROM contigs_sequence_features WHERE source='bact_t3_other'" "8" "test 10f: different --source-name on same DB coexists"

####################################################################################################
echo
echo "Sequence-features component tests 1, 2, 2.5, 2.6, 2.7, and 3–10 PASSED."

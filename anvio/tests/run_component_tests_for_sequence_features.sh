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

assert_query $output_dir/MIGRATED.db "SELECT COUNT(*) FROM feature_types" "6" "feature_types has 6 builtin rows"
assert_query $output_dir/MIGRATED.db "SELECT COUNT(*) FROM feature_types WHERE is_builtin=1" "6" "all 6 feature_types are builtin"
assert_query $output_dir/MIGRATED.db "SELECT GROUP_CONCAT(feature_type, ',') FROM (SELECT feature_type FROM feature_types ORDER BY feature_type)" "CDS,exon,gene,intron,mRNA,transcript" "feature_types contains the six expected names"
assert_query $output_dir/MIGRATED.db "SELECT has_dedicated_table FROM feature_types WHERE feature_type='CDS'" "1" "CDS feature type is flagged has_dedicated_table"
assert_query $output_dir/MIGRATED.db "SELECT has_dedicated_table FROM feature_types WHERE feature_type='gene'" "0" "gene feature type is not flagged has_dedicated_table"
assert_query $output_dir/MIGRATED.db "SELECT is_builtin FROM feature_types WHERE feature_type='transcript'" "1" "transcript feature type is a builtin"
assert_query $output_dir/MIGRATED.db "SELECT has_dedicated_table FROM feature_types WHERE feature_type='transcript'" "0" "transcript feature type does not have a dedicated table"

# the two synthesis-related columns and the partial index on derivation
assert_query $output_dir/MIGRATED.db "SELECT COUNT(*) FROM pragma_table_info('contigs_sequence_features') WHERE name='derivation'" "1" "contigs_sequence_features has the derivation column"
assert_query $output_dir/MIGRATED.db "SELECT COUNT(*) FROM pragma_table_info('contigs_sequence_features') WHERE name='derived_from_feature_id'" "1" "contigs_sequence_features has the derived_from_feature_id column"

for tbl in contigs_sequence_features feature_relationships feature_qualifiers CDS_features; do
    assert_query $output_dir/MIGRATED.db "SELECT COUNT(*) FROM $tbl" "0" "$tbl is empty after migration"
done

# invariant: untouched tables retain their pre-migration row counts
assert_query $output_dir/MIGRATED.db "SELECT COUNT(*) FROM contigs_basic_info" "$PRE_CONTIGS" "contigs_basic_info row count is preserved"
assert_query $output_dir/MIGRATED.db "SELECT COUNT(*) FROM splits_basic_info"  "$PRE_SPLITS"  "splits_basic_info row count is preserved"
assert_query $output_dir/MIGRATED.db "SELECT COUNT(*) FROM genes_in_contigs"   "$PRE_GENES"   "genes_in_contigs row count is preserved"

# verify the indexes anvi'o expects on the new tables actually exist
for idx in idx_pk_contigs_sequence_features idx_pk_feature_types idx_pk_feature_relationships idx_pk_feature_qualifiers idx_pk_CDS_features idx_contigs_sequence_features_contig_start idx_contigs_sequence_features_contig_stop idx_contigs_sequence_features_feature_type idx_contigs_sequence_features_group_id idx_contigs_sequence_features_gcid idx_contigs_sequence_features_derivation idx_feature_relationships_parent; do
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

assert_query $output_dir/FRESH.db "SELECT COUNT(*) FROM feature_types" "6" "feature_types has 6 builtin rows in fresh db"
assert_query $output_dir/FRESH.db "SELECT GROUP_CONCAT(feature_type, ',') FROM (SELECT feature_type FROM feature_types ORDER BY feature_type)" "CDS,exon,gene,intron,mRNA,transcript" "fresh db feature_types contains the six expected builtins"
for tbl in contigs_sequence_features feature_relationships feature_qualifiers CDS_features; do
    assert_query $output_dir/FRESH.db "SELECT COUNT(*) FROM $tbl" "0" "$tbl is empty in fresh db"
done

# fresh DB has the synthesis-related columns and partial index
assert_query $output_dir/FRESH.db "SELECT COUNT(*) FROM pragma_table_info('contigs_sequence_features') WHERE name='derivation'" "1" "fresh db contigs_sequence_features has the derivation column"
assert_query $output_dir/FRESH.db "SELECT COUNT(*) FROM pragma_table_info('contigs_sequence_features') WHERE name='derived_from_feature_id'" "1" "fresh db contigs_sequence_features has the derived_from_feature_id column"
assert_query $output_dir/FRESH.db "SELECT COUNT(*) FROM sqlite_master WHERE type='index' AND name='idx_contigs_sequence_features_derivation'" "1" "fresh db has the partial index on derivation"

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
# canonical-parent convention: 3-segment CDS → 3 rows in feature_relationships pointing to the canonical
# literal mRNA, and synthesis adds 3 more pointing to the canonical synthesized transcript (the
# literal mRNA's children are re-parented to the synthesized transcript that shadows it).
assert_query $output_dir/EUK_T4.db "SELECT COUNT(*) FROM feature_relationships fr JOIN contigs_sequence_features csf ON fr.child_feature_id = csf.feature_id JOIN contigs_sequence_features pcsf ON fr.parent_feature_id = pcsf.feature_id WHERE csf.feature_type='CDS' AND csf.direction='f' AND pcsf.feature_type='mRNA' AND pcsf.derivation IS NULL" "3" "test 4: 3 CDS→literal-mRNA relationship rows (one per CDS segment)"
assert_query $output_dir/EUK_T4.db "SELECT COUNT(*) FROM feature_relationships fr JOIN contigs_sequence_features csf ON fr.child_feature_id = csf.feature_id WHERE csf.feature_type='CDS' AND csf.direction='f'" "6" "test 4: 6 forward-CDS relationship rows total (3 to literal mRNA + 3 to synthesized transcript)"
assert_query $output_dir/EUK_T4.db "SELECT COUNT(DISTINCT fr.parent_feature_id) FROM feature_relationships fr JOIN contigs_sequence_features csf ON fr.child_feature_id = csf.feature_id JOIN contigs_sequence_features pcsf ON fr.parent_feature_id = pcsf.feature_id WHERE csf.feature_type='CDS' AND csf.direction='f' AND pcsf.feature_type='mRNA' AND pcsf.derivation IS NULL" "1" "test 4: all forward CDS segments share one canonical literal mRNA parent"
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

# literal origin-crossing gene → 2 segment rows; synthesis adds 2 transcript segments + 2 exons → 6 total
assert_query $output_dir/CIRC_T5.db "SELECT COUNT(*) FROM contigs_sequence_features WHERE derivation IS NULL" "2" "test 5: origin-crossing literal gene stored as exactly 2 rows"
assert_query $output_dir/CIRC_T5.db "SELECT COUNT(*) FROM contigs_sequence_features" "6" "test 5: total = 2 literal gene segments + 2 synthesized transcript segments + 2 synthesized exons"
assert_query $output_dir/CIRC_T5.db "SELECT start FROM contigs_sequence_features WHERE segment_order=0 AND feature_type='gene'" "899" "test 5: literal gene segment_order=0 is [899,1000)"
assert_query $output_dir/CIRC_T5.db "SELECT stop  FROM contigs_sequence_features WHERE segment_order=0 AND feature_type='gene'" "1000" "test 5: literal gene segment_order=0 stop is 1000"
assert_query $output_dir/CIRC_T5.db "SELECT start FROM contigs_sequence_features WHERE segment_order=1 AND feature_type='gene'" "0" "test 5: literal gene segment_order=1 is [0,200)"
assert_query $output_dir/CIRC_T5.db "SELECT stop  FROM contigs_sequence_features WHERE segment_order=1 AND feature_type='gene'" "200" "test 5: literal gene segment_order=1 stop is 200"
assert_query $output_dir/CIRC_T5.db "SELECT COUNT(DISTINCT feature_group_id) FROM contigs_sequence_features WHERE feature_type='gene'" "1" "test 5: literal gene segments share feature_group_id"
assert_query $output_dir/CIRC_T5.db "SELECT COUNT(DISTINCT external_id) FROM contigs_sequence_features" "1" "test 5: every row (literal + synthesized) shares external_id (locus_tag)"
# lone-gene-case synthesis with an origin-crossing source produces a multi-segment transcript mirroring the source
assert_query $output_dir/CIRC_T5.db "SELECT COUNT(*) FROM contigs_sequence_features WHERE feature_type='transcript' AND derivation='gene'" "2" "test 5: origin-crossing gene → 2-segment synthesized transcript"
assert_query $output_dir/CIRC_T5.db "SELECT start FROM contigs_sequence_features WHERE segment_order=0 AND feature_type='transcript'" "899" "test 5: synthesized transcript seg0 mirrors literal seg0 (start=899)"
assert_query $output_dir/CIRC_T5.db "SELECT COUNT(*) FROM contigs_sequence_features WHERE feature_type='exon' AND derivation='gene'" "2" "test 5: origin-crossing gene → 2 synthesized exons (one per gene segment)"

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

# malformed feature absent; healthy literal feature present (synthesized rows inherit external_id,
# so we constrain to derivation IS NULL to count only the literal gene row).
assert_query $output_dir/LIN_T6.db "SELECT COUNT(*) FROM contigs_sequence_features WHERE external_id='LINEAR_BAD_001'" "0" "test 6: malformed origin-crosser absent"
assert_query $output_dir/LIN_T6.db "SELECT COUNT(*) FROM contigs_sequence_features WHERE external_id='LINEAR_OK_001' AND derivation IS NULL" "1" "test 6: healthy literal feature still imported"
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

    # all literal features should still be imported regardless of GCID match outcome (synthesized
    # rows are excluded by the derivation filter)
    assert_query $output_dir/BACT_T9.db "SELECT COUNT(*) FROM contigs_sequence_features WHERE derivation IS NULL" "8" "test 9: all 8 literal features imported despite gene-call divergence"
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
assert_query $output_dir/BACT_T3.db "SELECT COUNT(*) FROM contigs_sequence_features WHERE source='bact_t3_other' AND derivation IS NULL" "8" "test 10f: different --source-name on same DB coexists (8 literal rows)"

####################################################################################################
# Test 11: bacterial synthesis (gene + CDS, no mRNA, no exon) — the CDS-only synthesis case
####################################################################################################
INFO "Test 11: bacterial synthesis (CDS-only case — gene + CDS, no mRNA / no exon)"

# reuse BACT_T3.db built earlier (4 genes + 4 CDSs, single-source 'bact_t3')
# every gene → 1 synthesized transcript (derivation='CDS') and 1 synthesized exon (derivation='CDS')
assert_query $output_dir/BACT_T3.db "SELECT COUNT(*) FROM contigs_sequence_features WHERE source='bact_t3' AND feature_type='transcript' AND derivation='CDS'" "4" "test 11: 4 synthesized transcripts (derivation='CDS'), one per gene"
assert_query $output_dir/BACT_T3.db "SELECT COUNT(*) FROM contigs_sequence_features WHERE source='bact_t3' AND feature_type='exon' AND derivation='CDS'" "4" "test 11: 4 synthesized exons (derivation='CDS'), one per gene"

# TST_001: gene at [99, 500), CDS at [99, 500), all forward. Verify the synthesized transcript
# inherits the gene's coords/direction/external_id and the synthesized exon inherits the CDS's.
assert_query $output_dir/BACT_T3.db "SELECT start, stop, direction, external_id FROM contigs_sequence_features WHERE source='bact_t3' AND external_id='TST_001' AND feature_type='transcript' AND derivation='CDS'" "99|500|f|TST_001" "test 11: TST_001 synthesized transcript matches gene coords and inherits external_id"
assert_query $output_dir/BACT_T3.db "SELECT start, stop, direction FROM contigs_sequence_features WHERE source='bact_t3' AND external_id='TST_001' AND feature_type='exon' AND derivation='CDS'" "99|500|f" "test 11: TST_001 synthesized exon matches CDS coords"

# derived_from_feature_id on the synthesized transcript points to the canonical CDS feature_id
CANON_CDS_TST001=$(sqlite3 $output_dir/BACT_T3.db "SELECT feature_id FROM contigs_sequence_features WHERE source='bact_t3' AND external_id='TST_001' AND feature_type='CDS' AND derivation IS NULL")
SYNTH_T_TST001_REF=$(sqlite3 $output_dir/BACT_T3.db "SELECT derived_from_feature_id FROM contigs_sequence_features WHERE source='bact_t3' AND external_id='TST_001' AND feature_type='transcript' AND derivation='CDS'")
[ "$CANON_CDS_TST001" = "$SYNTH_T_TST001_REF" ] || { echo "ASSERT FAILED [test 11: TST_001 transcript derived_from = canonical CDS fid]"; exit 1; }
echo "  OK [test 11: TST_001 synthesized transcript points to the canonical CDS feature_id]"

# new feature_relationships rows: synthesized-transcript→gene (part_of) and synthesized-exon→synthesized-transcript (part_of)
# (BACT_T3.db has two sources after test 10f, so we restrict to source='bact_t3' for cleanliness)
assert_query $output_dir/BACT_T3.db "SELECT COUNT(*) FROM feature_relationships fr JOIN contigs_sequence_features csf ON fr.child_feature_id = csf.feature_id JOIN contigs_sequence_features pcsf ON fr.parent_feature_id = pcsf.feature_id WHERE csf.derivation='CDS' AND csf.feature_type='transcript' AND csf.source='bact_t3' AND pcsf.feature_type='gene' AND pcsf.derivation IS NULL AND fr.relationship='part_of'" "4" "test 11: 4 synthesized-transcript → literal-gene part_of relationships (source=bact_t3)"
assert_query $output_dir/BACT_T3.db "SELECT COUNT(*) FROM feature_relationships fr JOIN contigs_sequence_features csf ON fr.child_feature_id = csf.feature_id JOIN contigs_sequence_features pcsf ON fr.parent_feature_id = pcsf.feature_id WHERE csf.derivation='CDS' AND csf.feature_type='exon' AND csf.source='bact_t3' AND pcsf.derivation='CDS' AND pcsf.feature_type='transcript' AND fr.relationship='part_of'" "4" "test 11: 4 synthesized-exon → synthesized-transcript part_of relationships (source=bact_t3)"

# the literal CDS now has two child-side entries in feature_relationships: the original
# `derives_from gene` link (from the parent-relationship resolution pass), and a new `part_of`
# link to the synthesized transcript (the CDS-only synthesis case links every CDS segment to
# the synthesized transcript it acted as the source for).
assert_query $output_dir/BACT_T3.db "SELECT COUNT(*) FROM feature_relationships fr JOIN contigs_sequence_features csf ON fr.child_feature_id = csf.feature_id WHERE csf.feature_type='CDS' AND csf.external_id='TST_001' AND csf.source='bact_t3'" "2" "test 11: TST_001 CDS has 2 child-side relationship rows (derives_from gene + part_of synthesized transcript)"

####################################################################################################
# Test 12: eukaryotic synthesis from mRNA (gene + mRNA(join) + CDS(join), no literal exons)
####################################################################################################
INFO "Test 12: eukaryotic synthesis from mRNA (transcript-source-child case — mRNA child, no literal exons)"

anvi-gen-contigs-database -f $files/sequence_features/eukaryotic_no_exon.fa \
                          -o $output_dir/EUKNX_T12.db \
                          --skip-gene-calling \
                          --project-name "EukaryoticNoExonTest12" \
                          -L 1000 --no-progress $thread_controller >/dev/null

anvi-import-genbank-features -c $output_dir/EUKNX_T12.db \
                             -i $files/sequence_features/eukaryotic_no_exon.gb \
                             --source-name 'eukno_t12' >/dev/null

# 2 genes; each gene → 1 synthesized transcript (derivation='mRNA'), N synthesized exons per mRNA segment
assert_query $output_dir/EUKNX_T12.db "SELECT COUNT(*) FROM contigs_sequence_features WHERE feature_type='transcript' AND derivation='mRNA'" "2" "test 12: 2 synthesized transcripts (one per literal mRNA)"
# 3 segments forward mRNA + 2 segments reverse mRNA = 5 synthesized exons total
assert_query $output_dir/EUKNX_T12.db "SELECT COUNT(*) FROM contigs_sequence_features WHERE feature_type='exon' AND derivation='mRNA'" "5" "test 12: 5 synthesized exons (3 forward mRNA segments + 2 reverse mRNA segments)"

# synthesized exons match mRNA segment coords exactly (NOT CDS coords — confirms mRNA priority)
# forward mRNA segments: [99, 500), [900, 1500), [2000, 2500)
assert_query $output_dir/EUKNX_T12.db "SELECT GROUP_CONCAT(start || '-' || stop, ',') FROM (SELECT start, stop FROM contigs_sequence_features WHERE feature_type='exon' AND derivation='mRNA' AND direction='f' ORDER BY start)" "99-500,900-1500,2000-2500" "test 12: forward synthesized exons match mRNA segments (not CDS)"
# reverse mRNA segments (transcription order: rightmost first): [3999, 4500), [2999, 3500)
assert_query $output_dir/EUKNX_T12.db "SELECT GROUP_CONCAT(start || '-' || stop, ',') FROM (SELECT start, stop FROM contigs_sequence_features WHERE feature_type='exon' AND derivation='mRNA' AND direction='r' ORDER BY start)" "2999-3500,3999-4500" "test 12: reverse synthesized exons match mRNA segments (not CDS)"

# synthesized exons are independent single-segment rows (feature_group_id IS NULL, segment_order IS NULL)
assert_query $output_dir/EUKNX_T12.db "SELECT COUNT(*) FROM contigs_sequence_features WHERE feature_type='exon' AND derivation='mRNA' AND feature_group_id IS NULL AND segment_order IS NULL" "5" "test 12: every synthesized exon is single-segment (NULL feature_group_id and segment_order)"

# direction inheritance verified by the forward/reverse splits above
# each CDS segment has a new part_of relationship to the synthesized transcript (canonical-parent convention)
assert_query $output_dir/EUKNX_T12.db "SELECT COUNT(*) FROM feature_relationships fr JOIN contigs_sequence_features csf ON fr.child_feature_id = csf.feature_id JOIN contigs_sequence_features pcsf ON fr.parent_feature_id = pcsf.feature_id WHERE csf.feature_type='CDS' AND pcsf.derivation='mRNA' AND pcsf.feature_type='transcript' AND fr.relationship='part_of'" "5" "test 12: 5 CDS-segment → synthesized-transcript part_of relationships (3+2)"

####################################################################################################
# Test 13: eukaryotic synthesis from CDS when mRNA absent (CDS-only case with multi-segment CDS)
####################################################################################################
INFO "Test 13: eukaryotic synthesis from CDS when mRNA absent (CDS-only case — gene + CDS, no mRNA, no exon)"

anvi-gen-contigs-database -f $files/sequence_features/eukaryotic_no_mrna.fa \
                          -o $output_dir/EUKNM_T13.db \
                          --skip-gene-calling \
                          --project-name "EukaryoticNoMRNATest13" \
                          -L 1000 --no-progress $thread_controller >/dev/null

anvi-import-genbank-features -c $output_dir/EUKNM_T13.db \
                             -i $files/sequence_features/eukaryotic_no_mrna.gb \
                             --source-name 'eukno_t13' >/dev/null

# 1 gene; 1 synthesized transcript with derivation='CDS' and coords equal to GENE (not CDS)
assert_query $output_dir/EUKNM_T13.db "SELECT COUNT(*) FROM contigs_sequence_features WHERE feature_type='transcript' AND derivation='CDS'" "1" "test 13: one synthesized transcript with derivation='CDS'"
assert_query $output_dir/EUKNM_T13.db "SELECT start, stop FROM contigs_sequence_features WHERE feature_type='transcript' AND derivation='CDS'" "99|2500" "test 13: synthesized transcript spans the gene's coords (not the CDS's)"

# derived_from_feature_id of the transcript points to the canonical CDS feature_id
CANON_CDS_NMR=$(sqlite3 $output_dir/EUKNM_T13.db "SELECT feature_id FROM contigs_sequence_features WHERE feature_type='CDS' AND derivation IS NULL AND segment_order=0")
SYNTH_T_NMR_REF=$(sqlite3 $output_dir/EUKNM_T13.db "SELECT derived_from_feature_id FROM contigs_sequence_features WHERE feature_type='transcript' AND derivation='CDS'")
[ "$CANON_CDS_NMR" = "$SYNTH_T_NMR_REF" ] || { echo "ASSERT FAILED [test 13: transcript derived_from = canonical CDS]"; exit 1; }
echo "  OK [test 13: synthesized transcript's derived_from_feature_id is the canonical CDS feature_id]"

# 3 synthesized exons, one per CDS segment, with coords matching each CDS segment
assert_query $output_dir/EUKNM_T13.db "SELECT COUNT(*) FROM contigs_sequence_features WHERE feature_type='exon' AND derivation='CDS'" "3" "test 13: 3 synthesized exons (one per CDS segment)"
assert_query $output_dir/EUKNM_T13.db "SELECT GROUP_CONCAT(start || '-' || stop, ',') FROM (SELECT start, stop FROM contigs_sequence_features WHERE feature_type='exon' AND derivation='CDS' ORDER BY start)" "200-700,1100-1600,2000-2400" "test 13: synthesized exon coords match CDS segments"

# each CDS segment has TWO child-side relationships: original derives_from gene + new part_of synthesized_transcript
assert_query $output_dir/EUKNM_T13.db "SELECT COUNT(*) FROM feature_relationships fr JOIN contigs_sequence_features csf ON fr.child_feature_id = csf.feature_id WHERE csf.feature_type='CDS' AND csf.derivation IS NULL" "6" "test 13: each of 3 CDS segments has 2 relationships (derives_from gene + part_of synth transcript)"

####################################################################################################
# Test 14: pseudogene synthesis from mRNA when CDS absent (transcript-source-child case, no CDS)
####################################################################################################
INFO "Test 14: pseudogene synthesis from mRNA when CDS absent"

anvi-gen-contigs-database -f $files/sequence_features/pseudogene.fa \
                          -o $output_dir/PSEUDO_T14.db \
                          --skip-gene-calling \
                          --project-name "PseudogeneTest14" \
                          -L 1000 --no-progress $thread_controller >/dev/null

anvi-import-genbank-features -c $output_dir/PSEUDO_T14.db \
                             -i $files/sequence_features/pseudogene.gb \
                             --source-name 'pseudo_t14' >/dev/null

# 1 synthesized transcript with derivation='mRNA' and coords = mRNA spanning extent
assert_query $output_dir/PSEUDO_T14.db "SELECT COUNT(*) FROM contigs_sequence_features WHERE feature_type='transcript' AND derivation='mRNA'" "1" "test 14: 1 synthesized transcript (derivation='mRNA')"
assert_query $output_dir/PSEUDO_T14.db "SELECT start, stop FROM contigs_sequence_features WHERE feature_type='transcript' AND derivation='mRNA'" "99|2500" "test 14: synthesized transcript spans the mRNA's extent"
# 3 synthesized exons matching mRNA segments
assert_query $output_dir/PSEUDO_T14.db "SELECT COUNT(*) FROM contigs_sequence_features WHERE feature_type='exon' AND derivation='mRNA'" "3" "test 14: 3 synthesized exons (one per mRNA segment)"
assert_query $output_dir/PSEUDO_T14.db "SELECT GROUP_CONCAT(start || '-' || stop, ',') FROM (SELECT start, stop FROM contigs_sequence_features WHERE feature_type='exon' AND derivation='mRNA' ORDER BY start)" "99-500,900-1500,2000-2500" "test 14: synthesized exon coords match mRNA segments"

####################################################################################################
# Test 15: non-mRNA RNA synthesis (gene + ncRNA, gene + tRNA)
####################################################################################################
INFO "Test 15: non-mRNA RNA synthesis (ncRNA + tRNA)"

anvi-gen-contigs-database -f $files/sequence_features/ncrna_trna.fa \
                          -o $output_dir/NCRNA_T15.db \
                          --skip-gene-calling \
                          --project-name "NcRNATRnaTest15" \
                          -L 1000 --no-progress $thread_controller >/dev/null

anvi-import-genbank-features -c $output_dir/NCRNA_T15.db \
                             -i $files/sequence_features/ncrna_trna.gb \
                             --source-name 'ncrna_t15' >/dev/null

# PARENT_RULES extension: every transcript-source RNA type (not just mRNA) now links to a
# gene parent with `part_of`, so ncRNA and tRNA each have a `part_of` link to their gene
assert_query $output_dir/NCRNA_T15.db "SELECT COUNT(*) FROM feature_relationships fr JOIN contigs_sequence_features csf ON fr.child_feature_id = csf.feature_id JOIN contigs_sequence_features pcsf ON fr.parent_feature_id = pcsf.feature_id WHERE csf.feature_type='ncRNA' AND pcsf.feature_type='gene' AND pcsf.derivation IS NULL AND fr.relationship='part_of'" "1" "test 15: ncRNA has part_of relationship to its gene parent"
assert_query $output_dir/NCRNA_T15.db "SELECT COUNT(*) FROM feature_relationships fr JOIN contigs_sequence_features csf ON fr.child_feature_id = csf.feature_id JOIN contigs_sequence_features pcsf ON fr.parent_feature_id = pcsf.feature_id WHERE csf.feature_type='tRNA' AND pcsf.feature_type='gene' AND pcsf.derivation IS NULL AND fr.relationship='part_of'" "1" "test 15: tRNA has part_of relationship to its gene parent"

# synthesis: one transcript per RNA child, derivation matches the source RNA type
assert_query $output_dir/NCRNA_T15.db "SELECT COUNT(*) FROM contigs_sequence_features WHERE feature_type='transcript' AND derivation='ncRNA'" "1" "test 15: 1 synthesized transcript with derivation='ncRNA'"
assert_query $output_dir/NCRNA_T15.db "SELECT COUNT(*) FROM contigs_sequence_features WHERE feature_type='transcript' AND derivation='tRNA'" "1" "test 15: 1 synthesized transcript with derivation='tRNA'"
# transcripts span the RNA's extent and inherit external_id
assert_query $output_dir/NCRNA_T15.db "SELECT start, stop, external_id FROM contigs_sequence_features WHERE feature_type='transcript' AND derivation='ncRNA'" "99|400|NC_001" "test 15: ncRNA-derived transcript spans the ncRNA's extent and inherits external_id"
assert_query $output_dir/NCRNA_T15.db "SELECT start, stop, external_id FROM contigs_sequence_features WHERE feature_type='transcript' AND derivation='tRNA'" "799|900|TR_001" "test 15: tRNA-derived transcript spans the tRNA's extent and inherits external_id"

# synthesis: one exon per RNA child (since RNAs are single-segment, one exon = RNA span)
assert_query $output_dir/NCRNA_T15.db "SELECT COUNT(*) FROM contigs_sequence_features WHERE feature_type='exon' AND derivation='ncRNA'" "1" "test 15: 1 synthesized exon with derivation='ncRNA'"
assert_query $output_dir/NCRNA_T15.db "SELECT COUNT(*) FROM contigs_sequence_features WHERE feature_type='exon' AND derivation='tRNA'" "1" "test 15: 1 synthesized exon with derivation='tRNA'"

####################################################################################################
# Test 16: lone gene synthesis (gene only, no children) — the lone-gene synthesis case
####################################################################################################
INFO "Test 16: lone gene synthesis (no children)"

anvi-gen-contigs-database -f $files/sequence_features/lone_gene.fa \
                          -o $output_dir/LONE_T16.db \
                          --skip-gene-calling \
                          --project-name "LoneGeneTest16" \
                          -L 1000 --no-progress $thread_controller >/dev/null

anvi-import-genbank-features -c $output_dir/LONE_T16.db \
                             -i $files/sequence_features/lone_gene.gb \
                             --source-name 'lone_t16' >/dev/null

# 1 gene + 1 transcript (derivation='gene') + 1 exon (derivation='gene') = 3 total
assert_query $output_dir/LONE_T16.db "SELECT COUNT(*) FROM contigs_sequence_features" "3" "test 16: lone gene produces 3 rows total (gene + transcript + exon)"
assert_query $output_dir/LONE_T16.db "SELECT start, stop, direction FROM contigs_sequence_features WHERE feature_type='transcript' AND derivation='gene'" "199|800|f" "test 16: synthesized transcript spans the gene"
assert_query $output_dir/LONE_T16.db "SELECT start, stop, direction FROM contigs_sequence_features WHERE feature_type='exon' AND derivation='gene'" "199|800|f" "test 16: synthesized exon spans the gene"

# relationships: gene ← transcript ← exon
assert_query $output_dir/LONE_T16.db "SELECT COUNT(*) FROM feature_relationships fr JOIN contigs_sequence_features csf ON fr.child_feature_id = csf.feature_id JOIN contigs_sequence_features pcsf ON fr.parent_feature_id = pcsf.feature_id WHERE csf.derivation='gene' AND csf.feature_type='transcript' AND pcsf.feature_type='gene' AND pcsf.derivation IS NULL" "1" "test 16: synthesized transcript → literal gene (part_of)"
assert_query $output_dir/LONE_T16.db "SELECT COUNT(*) FROM feature_relationships fr JOIN contigs_sequence_features csf ON fr.child_feature_id = csf.feature_id JOIN contigs_sequence_features pcsf ON fr.parent_feature_id = pcsf.feature_id WHERE csf.derivation='gene' AND csf.feature_type='exon' AND pcsf.derivation='gene' AND pcsf.feature_type='transcript'" "1" "test 16: synthesized exon → synthesized transcript (part_of)"

####################################################################################################
# Test 17: literal exons present (gene + mRNA(join) + exon × N + CDS(join))
####################################################################################################
INFO "Test 17: literal exons present (reuse eukaryotic fixture; no synthesized exons)"

# reuse EUK_T4.db built earlier
# the synthesized transcript exists (transcripts are always synthesized for uniformity)
assert_query $output_dir/EUK_T4.db "SELECT COUNT(*) FROM contigs_sequence_features WHERE feature_type='transcript' AND derivation='mRNA'" "2" "test 17: 2 synthesized transcripts (one per gene), derivation='mRNA'"
# no synthesized exons since literal exons exist
assert_query $output_dir/EUK_T4.db "SELECT COUNT(*) FROM contigs_sequence_features WHERE feature_type='exon' AND derivation IS NOT NULL" "0" "test 17: zero synthesized exons (literal exons present)"
# total exon rows = literal exon count (5 segments: 3 forward + 2 reverse)
assert_query $output_dir/EUK_T4.db "SELECT COUNT(*) FROM contigs_sequence_features WHERE feature_type='exon'" "5" "test 17: 5 literal exon rows total (all derivation IS NULL)"
assert_query $output_dir/EUK_T4.db "SELECT COUNT(*) FROM contigs_sequence_features WHERE feature_type='exon' AND derivation IS NULL" "5" "test 17: all 5 exon rows are literal (derivation IS NULL)"

# each literal exon has TWO entries in feature_relationships as a child:
# (1) part_of literal mRNA (from the parent-relationship resolution pass),
# (2) part_of synthesized transcript (from synthesis's re-parenting of literal children).
assert_query $output_dir/EUK_T4.db "SELECT COUNT(*) FROM feature_relationships fr JOIN contigs_sequence_features csf ON fr.child_feature_id = csf.feature_id WHERE csf.feature_type='exon' AND fr.relationship='part_of'" "10" "test 17: 10 exon part_of rows (5 to literal mRNA + 5 to synthesized transcript)"
assert_query $output_dir/EUK_T4.db "SELECT COUNT(*) FROM feature_relationships fr JOIN contigs_sequence_features csf ON fr.child_feature_id = csf.feature_id JOIN contigs_sequence_features pcsf ON fr.parent_feature_id = pcsf.feature_id WHERE csf.feature_type='exon' AND pcsf.feature_type='mRNA' AND pcsf.derivation IS NULL" "5" "test 17: 5 exon → literal mRNA part_of rows"
assert_query $output_dir/EUK_T4.db "SELECT COUNT(*) FROM feature_relationships fr JOIN contigs_sequence_features csf ON fr.child_feature_id = csf.feature_id JOIN contigs_sequence_features pcsf ON fr.parent_feature_id = pcsf.feature_id WHERE csf.feature_type='exon' AND pcsf.derivation='mRNA' AND pcsf.feature_type='transcript'" "5" "test 17: 5 exon → synthesized transcript part_of rows"

####################################################################################################
# Test 18: multi-isoform gene (alternative splicing — two mRNAs with distinct locus_tags, shared exon)
####################################################################################################
INFO "Test 18: multi-isoform gene with shared exon coordinate"

anvi-gen-contigs-database -f $files/sequence_features/multi_isoform.fa \
                          -o $output_dir/ISO_T18.db \
                          --skip-gene-calling \
                          --project-name "MultiIsoformTest18" \
                          -L 1000 --no-progress $thread_controller >/dev/null

anvi-import-genbank-features -c $output_dir/ISO_T18.db \
                             -i $files/sequence_features/multi_isoform.gb \
                             --source-name 'iso_t18' >/dev/null

# both literal mRNAs survive (verifying the 9-field hash with external_id disambiguates them
# despite the shared exon coordinate at [99, 500))
assert_query $output_dir/ISO_T18.db "SELECT COUNT(DISTINCT external_id) FROM contigs_sequence_features WHERE feature_type='mRNA' AND derivation IS NULL" "2" "test 18: both literal mRNAs survive (distinct external_id values)"
assert_query $output_dir/ISO_T18.db "SELECT COUNT(*) FROM contigs_sequence_features WHERE feature_type='mRNA' AND derivation IS NULL AND external_id='ISO_A_001'" "3" "test 18: isoform A literal mRNA has 3 segments"
assert_query $output_dir/ISO_T18.db "SELECT COUNT(*) FROM contigs_sequence_features WHERE feature_type='mRNA' AND derivation IS NULL AND external_id='ISO_B_001'" "3" "test 18: isoform B literal mRNA has 3 segments"

# two synthesized transcripts, one per literal mRNA, different external_id and different derived_from
assert_query $output_dir/ISO_T18.db "SELECT COUNT(*) FROM contigs_sequence_features WHERE feature_type='transcript' AND derivation='mRNA'" "2" "test 18: 2 synthesized transcripts (one per isoform)"
assert_query $output_dir/ISO_T18.db "SELECT COUNT(DISTINCT derived_from_feature_id) FROM contigs_sequence_features WHERE feature_type='transcript' AND derivation='mRNA'" "2" "test 18: synthesized transcripts have distinct derived_from_feature_id values"
assert_query $output_dir/ISO_T18.db "SELECT COUNT(DISTINCT external_id) FROM contigs_sequence_features WHERE feature_type='transcript' AND derivation='mRNA'" "2" "test 18: synthesized transcripts inherit distinct external_id values"

# each transcript gets its own synthesized exon set: 3+3 = 6 synthesized exons
assert_query $output_dir/ISO_T18.db "SELECT COUNT(*) FROM contigs_sequence_features WHERE feature_type='exon' AND derivation='mRNA'" "6" "test 18: 6 synthesized exons total (3 per isoform)"
# the shared exon coordinate at [99, 500) appears in BOTH isoforms — two distinct synthesized exon rows
assert_query $output_dir/ISO_T18.db "SELECT COUNT(*) FROM contigs_sequence_features WHERE feature_type='exon' AND derivation='mRNA' AND start=99 AND stop=500" "2" "test 18: shared exon coordinate [99,500) produces 2 distinct synthesized exon rows"
assert_query $output_dir/ISO_T18.db "SELECT COUNT(DISTINCT external_id) FROM contigs_sequence_features WHERE feature_type='exon' AND derivation='mRNA' AND start=99 AND stop=500" "2" "test 18: the two shared-coord exons have distinct external_id values (one per isoform)"
assert_query $output_dir/ISO_T18.db "SELECT COUNT(DISTINCT feature_id) FROM contigs_sequence_features WHERE feature_type='exon' AND derivation='mRNA'" "6" "test 18: 6 distinct synthesized-exon feature_ids (no deduplication across isoforms)"

####################################################################################################
# Test 19: determinism with synthesis (extends test 7's logic to verify synthesized rows are deterministic)
####################################################################################################
INFO "Test 19: determinism with synthesis"

# test 7 already verified determinism by comparing the ORDER-BY-feature_id GROUP_CONCAT of feature_ids
# across --force replays. With synthesis on, those replays include synthesized rows as well, so
# test 7's determinism check now covers the full set. The assertions below explicitly verify
# the synthesized rows ARE in BACT_T7's alpha set and that their feature_ids are stable.
N_SYNTH_BEFORE_FORCE=$(sqlite3 $output_dir/BACT_T7.db "SELECT COUNT(*) FROM contigs_sequence_features WHERE source='alpha' AND derivation IS NOT NULL")
IDS_SYNTH_BEFORE_FORCE=$(sqlite3 $output_dir/BACT_T7.db "SELECT GROUP_CONCAT(feature_id, ',') FROM (SELECT feature_id FROM contigs_sequence_features WHERE source='alpha' AND derivation IS NOT NULL ORDER BY feature_id)")

# one more --force replay; synthesized feature_ids should be identical
anvi-import-genbank-features -c $output_dir/BACT_T7.db -i $files/sequence_features/bacterial.gb --source-name 'alpha' --force >/dev/null
N_SYNTH_AFTER_FORCE=$(sqlite3 $output_dir/BACT_T7.db "SELECT COUNT(*) FROM contigs_sequence_features WHERE source='alpha' AND derivation IS NOT NULL")
IDS_SYNTH_AFTER_FORCE=$(sqlite3 $output_dir/BACT_T7.db "SELECT GROUP_CONCAT(feature_id, ',') FROM (SELECT feature_id FROM contigs_sequence_features WHERE source='alpha' AND derivation IS NOT NULL ORDER BY feature_id)")

[ "$N_SYNTH_BEFORE_FORCE" = "$N_SYNTH_AFTER_FORCE" ] || { echo "ASSERT FAILED [test 19: synthesized row count changed after --force]"; exit 1; }
[ "$IDS_SYNTH_BEFORE_FORCE" = "$IDS_SYNTH_AFTER_FORCE" ] || { echo "ASSERT FAILED [test 19: synthesized feature_ids differ after --force; not deterministic]"; exit 1; }
[ "$N_SYNTH_BEFORE_FORCE" -gt 0 ] || { echo "ASSERT FAILED [test 19: expected synthesized rows under source='alpha']"; exit 1; }
echo "  OK [test 19: synthesized feature_ids are deterministic across --force replays ($N_SYNTH_AFTER_FORCE synthesized rows)]"

####################################################################################################
echo
echo "Sequence-features component tests 1, 2, 2.5, 2.6, 2.7, and 3–19 PASSED."

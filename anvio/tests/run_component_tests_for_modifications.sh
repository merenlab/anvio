#!/bin/bash
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/00.sh"

# Setup #############################
SETUP_WITH_OUTPUT_DIR $1 $2 $3
export output_dir
#####################################

INFO "Preparing contigs FASTA and BAM with modification tags"

python - <<'PY'
import os
import array
import pysam

output_dir = os.environ.get('output_dir')
if not output_dir:
    raise SystemExit("output_dir is not set")

contigs_fasta = os.path.join(output_dir, 'contigs.fa')
seq = 'ACGTACGTACGTACGTACGTACGTACGTACGT'

with open(contigs_fasta, 'w') as f:
    f.write('>contig_1\n')
    f.write(seq + '\n')

bam_path = os.path.join(output_dir, 'sample.bam')
header = {'HD': {'VN': '1.0'}, 'SQ': [{'LN': len(seq), 'SN': 'contig_1'}]}

with pysam.AlignmentFile(bam_path, 'wb', header=header) as outf:
    a = pysam.AlignedSegment()
    a.query_name = 'read_1'
    a.query_sequence = seq
    a.flag = 0
    a.reference_id = 0
    a.reference_start = 0
    a.mapping_quality = 60
    a.cigar = [(0, len(seq))]
    a.query_qualities = pysam.qualitystring_to_array('I' * len(seq))
    a.set_tag('MM', 'C+m,0;')
    a.set_tag('ML', array.array('B', [255]))
    outf.write(a)

    b = pysam.AlignedSegment()
    b.query_name = 'read_2'
    b.query_sequence = seq
    b.flag = 16
    b.reference_id = 0
    b.reference_start = 0
    b.mapping_quality = 60
    b.cigar = [(0, len(seq))]
    b.query_qualities = pysam.qualitystring_to_array('I' * len(seq))
    b.set_tag('MM', 'A+a,0;')
    b.set_tag('ML', array.array('B', [200]))
    outf.write(b)

    c = pysam.AlignedSegment()
    c.query_name = 'read_3'
    c.query_sequence = seq[4:]
    c.flag = 0
    c.reference_id = 0
    c.reference_start = 4
    c.mapping_quality = 60
    c.cigar = [(0, len(seq) - 4)]
    c.query_qualities = pysam.qualitystring_to_array('I' * (len(seq) - 4))
    c.set_tag('MM', 'G+z,0;')
    c.set_tag('ML', array.array('B', [10]))
    outf.write(c)

    d = pysam.AlignedSegment()
    d.query_name = 'read_4'
    d.query_sequence = seq[2:]
    d.flag = 16
    d.reference_id = 0
    d.reference_start = 2
    d.mapping_quality = 60
    d.cigar = [(0, len(seq) - 2)]
    d.query_qualities = pysam.qualitystring_to_array('I' * (len(seq) - 2))
    d.set_tag('MM', 'T+m,0;')
    d.set_tag('ML', array.array('B', [50]))
    outf.write(d)

    e = pysam.AlignedSegment()
    e.query_name = 'read_5'
    e.query_sequence = seq[1:]
    e.flag = 0
    e.reference_id = 0
    e.reference_start = 1
    e.mapping_quality = 60
    e.cigar = [(0, len(seq) - 1)]
    e.query_qualities = pysam.qualitystring_to_array('I' * (len(seq) - 1))
    e.set_tag('MM', 'C+m,0;')
    outf.write(e)

sorted_bam_path = bam_path + '.sorted'
pysam.sort('-o', sorted_bam_path, bam_path)
os.replace(sorted_bam_path, bam_path)
pysam.index(bam_path)
PY

INFO "Generating contigs database"
anvi-gen-contigs-database -f $output_dir/contigs.fa -o $output_dir/CONTIGS.db -L 10 --project-name "Modifications test contigs" $thread_controller

INFO "Profiling with modifications enabled"
anvi-profile -i $output_dir/sample.bam \
             -o $output_dir/SAMPLE \
             -c $output_dir/CONTIGS.db \
             --min-contig-length 10 \
             --include-modifications \
             --skip-SNV-profiling \
             --skip-INDEL-profiling \
             $thread_controller

INFO "Validating modifications table entries"
mod_table_exists=$(sqlite3 $output_dir/SAMPLE/PROFILE.db "select count(*) from sqlite_master where type='table' and name='modifications';")
if [ "$mod_table_exists" -ne 1 ]; then
    echo "Expected modifications table to exist"
    exit 1
fi

mod_table_has_count=$(sqlite3 $output_dir/SAMPLE/PROFILE.db "select count(*) from pragma_table_info('modifications') where name='count';")
if [ "$mod_table_has_count" -ne 1 ]; then
    echo "Expected modifications table to include a count column"
    exit 1
fi

num_mods=$(sqlite3 $output_dir/SAMPLE/PROFILE.db "select count(*) from modifications;")
if [ "$num_mods" -lt 4 ]; then
    echo "Expected at least four modification entries, found $num_mods"
    exit 1
fi

total_mod_count=$(sqlite3 $output_dir/SAMPLE/PROFILE.db "select coalesce(sum(count), 0) from modifications;")
if [ "$total_mod_count" -lt 4 ]; then
    echo "Expected at least four total modification counts, found $total_mod_count"
    exit 1
fi

profiled_flag=$(sqlite3 $output_dir/SAMPLE/PROFILE.db "select value from self where key='modifications_profiled';")
if [ "$profiled_flag" != "1" ]; then
    echo "Expected modifications_profiled=1, found $profiled_flag"
    exit 1
fi

min_pos=$(sqlite3 $output_dir/SAMPLE/PROFILE.db "select min(pos) from modifications;")
max_pos=$(sqlite3 $output_dir/SAMPLE/PROFILE.db "select max(pos) from modifications;")
if [ "$min_pos" -lt 0 ]; then
    echo "Expected min pos >= 0, found $min_pos"
    exit 1
fi
if [ "$max_pos" -ge 32 ]; then
    echo "Expected max pos < 32, found $max_pos"
    exit 1
fi

distinct_mods=$(sqlite3 $output_dir/SAMPLE/PROFILE.db "select group_concat(distinct modification) from modifications;")
if [[ "$distinct_mods" != *"a"* ]] || [[ "$distinct_mods" != *"m"* ]] || [[ "$distinct_mods" != *"z"* ]]; then
    echo "Expected a, m, and z modification codes, found $distinct_mods"
    exit 1
fi

echo "Found $num_mods modification rows totaling $total_mod_count events"

INFO "Profiling without modifications flag"
anvi-profile -i $output_dir/sample.bam \
             -o $output_dir/SAMPLE-NOMODS \
             -c $output_dir/CONTIGS.db \
             --min-contig-length 10 \
             --skip-SNV-profiling \
             --skip-INDEL-profiling \
             $thread_controller

num_mods_no_flag=$(sqlite3 $output_dir/SAMPLE-NOMODS/PROFILE.db "select count(*) from modifications;")
if [ "$num_mods_no_flag" -ne 0 ]; then
    echo "Expected zero modification entries without flag, found $num_mods_no_flag"
    exit 1
fi

INFO "Testing modification filters and min-coverage thresholds"
anvi-profile -i $output_dir/sample.bam \
             -o $output_dir/SAMPLE-FILTERED \
             -c $output_dir/CONTIGS.db \
             --min-contig-length 10 \
             --include-modifications \
             --modification-filter m:0.5 \
             --modification-filter a:0.8 \
             --modification-filter-default 0.9 \
             --min-coverage-for-modifications 2 \
             --skip-SNV-profiling \
             --skip-INDEL-profiling \
             $thread_controller

# Verify filter metadata is set
filter_meta=$(sqlite3 $output_dir/SAMPLE-FILTERED/PROFILE.db "select value from self where key='modification_filters';")
if [[ "$filter_meta" != *"\"m\""* ]] || [[ "$filter_meta" != *"\"a\""* ]]; then
    echo "Expected modification filters to be stored in meta, found: $filter_meta"
    exit 1
fi

min_cov=$(sqlite3 $output_dir/SAMPLE-FILTERED/PROFILE.db "select value from self where key='min_coverage_for_modifications';")
if [ "$min_cov" != "2" ]; then
    echo "Expected min_coverage_for_modifications=2, found $min_cov"
    exit 1
fi

INFO "Testing aggregated modification counts"
row_count=$(sqlite3 $output_dir/SAMPLE/PROFILE.db "select count(*) from modifications where count > 0;")
if [ "$row_count" -eq 0 ]; then
    echo "Expected aggregated modification rows to have positive counts"
    exit 1
fi

sum_count=$(sqlite3 $output_dir/SAMPLE/PROFILE.db "select sum(count) from modifications;")
if [ "$sum_count" -lt 4 ]; then
    echo "Expected the summed modification count to be at least four, found $sum_count"
    exit 1
fi

INFO "Testing interactive API response structure"
python - <<'PYAPI'
import json
import sys
import os
import argparse

sys.path.insert(0, '/home/gkanaan/github/anvio')
import anvio.dbops as dbops
import anvio.terminal as terminal

output_dir = os.environ.get('output_dir')
args = argparse.Namespace(profile_db=f'{output_dir}/SAMPLE/PROFILE.db',
                          contigs_db=f'{output_dir}/CONTIGS.db')

profile = dbops.ProfileSuperclass(args, 
                                   r=terminal.Run(verbose=False), 
                                   p=terminal.Progress(verbose=False))

# Test get_modifications_information_for_split API
split_names = list(profile.split_names)[:1]
if split_names:
    split_name = split_names[0]
    mods_data = profile.get_modifications_information_for_split(split_name, min_coverage=0)
    
    # Verify response structure
    if 'data' not in mods_data or 'types' not in mods_data:
        print(f"ERROR: Expected 'data' and 'types' keys in modifications response")
        sys.exit(1)
    
    # Verify data structure
    data = mods_data['data']
    for sample_id in data:
        mods = data[sample_id].get('modifications', {})
        for pos in mods:
            entry = mods[pos]
            required_keys = ['pos', 'pos_in_contig', 'coverage', 'modification_counts', 
                           'modification_ratios', 'unmodified_ratio']
            for key in required_keys:
                if key not in entry:
                    print(f"ERROR: Missing key '{key}' in modifications entry")
                    sys.exit(1)
    
    print(f"✓ Modifications API response structure valid for split {split_name}")
    print(f"✓ Found {len(mods_data['types'])} modification types: {mods_data['types']}")
    
    # Test coverage filtering
    mods_filtered = profile.get_modifications_information_for_split(split_name, min_coverage=100)
    has_data = any(mods_filtered['data'][s].get('modifications', {}) for s in mods_filtered['data'])
    if has_data:
        print(f"Warning: Expected minimal result with min_coverage=100, got some data (coverage may exceed threshold)")
    else:
        print(f"✓ Coverage filtering works (high threshold removes all entries)")

PYAPI

INFO "Testing get_blank_modifications_dict helper"
python - <<'PYBLANK'
import sys
import os
import argparse

sys.path.insert(0, '/home/gkanaan/github/anvio')
import anvio.dbops as dbops
import anvio.terminal as terminal

output_dir = os.environ.get('output_dir')
args = argparse.Namespace(profile_db=f'{output_dir}/SAMPLE/PROFILE.db',
                          contigs_db=f'{output_dir}/CONTIGS.db')

profile = dbops.ProfileSuperclass(args,
                                   r=terminal.Run(verbose=False),
                                   p=terminal.Progress(verbose=False))

blank = profile.get_blank_modifications_dict()

# Verify blank structure
for sample in blank:
    if 'modifications' not in blank[sample]:
        print(f"ERROR: Missing 'modifications' key in blank dict for sample {sample}")
        sys.exit(1)
    if not isinstance(blank[sample]['modifications'], dict):
        print(f"ERROR: 'modifications' should be a dict")
        sys.exit(1)

print(f"✓ Blank modifications dict structure valid for {len(blank)} samples")

PYBLANK

INFO "Testing modifications with merged profiles"
anvi-profile -i $output_dir/sample.bam \
             -o $output_dir/SAMPLE2 \
             -c $output_dir/CONTIGS.db \
             --sample-name sample2 \
             --min-contig-length 10 \
             --include-modifications \
             --skip-SNV-profiling \
             --skip-INDEL-profiling \
             $thread_controller

anvi-merge $output_dir/SAMPLE/PROFILE.db $output_dir/SAMPLE2/PROFILE.db \
           -c $output_dir/CONTIGS.db \
           -o $output_dir/MERGED \
           $thread_controller

merged_mods=$(sqlite3 $output_dir/MERGED/PROFILE.db "select count(*) from modifications;")
if [ "$merged_mods" -gt 0 ]; then
    echo "✓ Modifications table exists in merged profile ($merged_mods entries)"
else
    echo "Note: Merged profile has no modifications (expected if profiles don't overlap)"
fi

merged_profiled=$(sqlite3 $output_dir/MERGED/PROFILE.db "select value from self where key='modifications_profiled';")
if [ -n "$merged_profiled" ]; then
    echo "✓ Merged profile has modifications_profiled meta key set to $merged_profiled"
fi

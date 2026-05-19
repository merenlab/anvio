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

num_mods=$(sqlite3 $output_dir/SAMPLE/PROFILE.db "select count(*) from modifications;")
if [ "$num_mods" -lt 4 ]; then
    echo "Expected at least four modification entries, found $num_mods"
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

distinct_bases=$(sqlite3 $output_dir/SAMPLE/PROFILE.db "select group_concat(distinct base) from modifications;")
if [[ "$distinct_bases" != *"A"* ]] || [[ "$distinct_bases" != *"C"* ]] || [[ "$distinct_bases" != *"G"* ]] || [[ "$distinct_bases" != *"T"* ]]; then
    echo "Expected A, C, G, and T bases, found $distinct_bases"
    exit 1
fi

distinct_mods=$(sqlite3 $output_dir/SAMPLE/PROFILE.db "select group_concat(distinct modification) from modifications;")
if [[ "$distinct_mods" != *"a"* ]] || [[ "$distinct_mods" != *"m"* ]] || [[ "$distinct_mods" != *"z"* ]]; then
    echo "Expected a, m, and z modification codes, found $distinct_mods"
    exit 1
fi

echo "Found $num_mods modification entries"

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

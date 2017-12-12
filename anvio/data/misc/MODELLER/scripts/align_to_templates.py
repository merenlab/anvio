# Step 2: prepare an alignment of all template structures and the
#         target sequence
#
# Align all of the best template structures detected in the previous step
# and compare them. Then align the target sequence with this block of
# aligned structures to generate an alignment suitable for modeling.
import sys
TARGET_PIR            = sys.argv[1]
TEMPLATE_IDS_FILE     = sys.argv[2]
ALIGNMENT_OUT_PIR     = sys.argv[3]
ALIGNMENT_OUT_PAP     = sys.argv[4]
PROTEIN_FAMILY_MATRIX = sys.argv[5]

# read in file written by MODELLER.run_compare()
template_ids = [tuple(x.strip().split("\t")) for x in open(TEMPLATE_IDS_FILE).readlines()]

from modeller import *

env = environ()
env.io.atom_files_directory = ['./TEMPLATE_PDBS']

# Create an alignment of the 'A' chains of 1clf, 1dur, 1fca and 2fdn
aln = alignment(env)
for (pdb, chain) in template_ids:
    m = model(env, file=pdb, model_segment=('FIRST:'+chain, 'LAST:'+chain))
    aln.append_model(m, atom_files=pdb, align_codes=pdb+chain)

# Structurally align all four templates and compare them
aln.malign()
aln.malign3d()
aln.compare_structures()
aln.id_table(matrix_file=PROTEIN_FAMILY_MATRIX)
env.dendrogram(matrix_file=PROTEIN_FAMILY_MATRIX, cluster_cut=-1.0)

# Align the target sequence to the previously-aligned structures
align_block = len(aln)
aln.append(file=TARGET_PIR)
aln.align2d(align_block=align_block, max_gap_length=50)
aln.write(file=ALIGNMENT_OUT_PIR, alignment_format='PIR')
aln.write(file=ALIGNMENT_OUT_PAP, alignment_format='PAP')

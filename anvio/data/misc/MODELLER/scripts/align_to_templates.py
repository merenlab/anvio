# Step 2: prepare an alignment of all template structures and the
#         target sequence
#
# Align all of the best template structures detected in the previous step
# and compare them. Then align the target sequence with this block of
# aligned structures to generate an alignment suitable for modeling.
import sys

TARGET_PIR = sys.argv[1]
TARGET_ID = sys.argv[2]
TEMPLATE_IDS_FILE = sys.argv[3]
ALIGNMENT_OUT_PIR = sys.argv[4]
ALIGNMENT_OUT_PAP = sys.argv[5]
PROTEIN_FAMILY_MATRIX = sys.argv[6]

# read in file written by MODELLER.run_compare()
template_ids = [
    tuple(x.strip().split("\t")) for x in open(TEMPLATE_IDS_FILE).readlines()
]

from modeller import *

env = environ()
env.io.atom_files_directory = ["./%s_TEMPLATE_PDBS" % TARGET_ID]

aln = alignment(env)
for pdb, chain in template_ids:
    m = model(env, file=pdb, model_segment=("FIRST:" + chain, "LAST:" + chain))
    aln.append_model(m, atom_files=pdb, align_codes=pdb + chain)

# Structurally align all four templates and compare them
aln.malign()
# gap_penalties_3d[1] is a gap extension penalty, say 1.75. Pairs of positions are considered
# equivalent when they have their selected atoms at most 2 times gap_penalties_3d[1] angstroms apart
# in the current superposition. If too strict, alignment can sometimes fail, so the distance to be
# considered equivalent increases iteratively until the procedure does not fail.
d = 1.75
while d < 10.00:
    try:
        aln.malign3d(gap_penalties_3d=(0.0, d))
        break
    except:
        d += 1.00

aln.compare_structures()
aln.id_table(matrix_file=PROTEIN_FAMILY_MATRIX)

# if there is only one template this will fail
if len(template_ids) > 1:
    env.dendrogram(matrix_file=PROTEIN_FAMILY_MATRIX, cluster_cut=-1.0)

# Align the target sequence to the previously-aligned structures
align_block = len(aln)
aln.append(file=TARGET_PIR)
aln.align2d(align_block=align_block, max_gap_length=50)
aln.write(file=ALIGNMENT_OUT_PIR, alignment_format="PIR")
aln.write(file=ALIGNMENT_OUT_PAP, alignment_format="PAP")

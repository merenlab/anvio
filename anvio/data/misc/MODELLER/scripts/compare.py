# Step 2: prepare an alignment of all template structures and the
#         target sequence
#
# Align all of the best template structures detected in the previous step
# and compare them. Then align the target sequence with this block of
# aligned structures to generate an alignment suitable for modeling.


from modeller import *


env = environ()
env.io.atom_files_directory = ['../atom_files']

# Create an alignment of the 'A' chains of 1clf, 1dur, 1fca and 2fdn
aln = alignment(env)
for (pdb, chain) in (('1clf', 'A'), ('1dur', 'A'), ('1fca', 'A'),
                     ('2fdn', 'A')):
    m = model(env, file=pdb, model_segment=('FIRST:'+chain, 'LAST:'+chain))
    aln.append_model(m, atom_files=pdb, align_codes=pdb+chain)

# Structurally align all four templates and compare them
aln.malign()
aln.malign3d()
aln.compare_structures()
aln.id_table(matrix_file='family.mat')
env.dendrogram(matrix_file='family.mat', cluster_cut=-1.0)

# Align the target sequence to the previously-aligned structures
align_block = len(aln)
aln.append(file='my.pir')
aln.align2d(align_block=align_block, max_gap_length=50)
aln.write(file='alignment.ali')

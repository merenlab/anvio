"""
Positional arguments:
1. INPUT - Target protein sequence (.pir format)
2. INTPUT - Database (.bin format)
3. OUTPUT - Alignment table
"""
import sys
TARGET = sys.argv[1]
DATABASE = sys.argv[2]
ALIGNMENT_TABLE = sys.argv[3]

from modeller import *
import os

log.verbose()
env = environ()

# Read in the database of PDB chains clustered at 95% sequence identity
print(os.path.isfile(DATABASE))
sdb = sequence_db(env)
sdb.read(seq_database_file=DATABASE, seq_database_format='BINARY', chains_list='ALL')

# Read in the target sequence in PIR alignment format
aln = alignment(env)
aln.append(file=TARGET, alignment_format='PIR', align_codes='ALL')

# Convert the target sequence from alignment to profile format
prf = aln.to_profile()

# Scan sequence database to pick up homologous sequences
prf.build(sdb, matrix_offset=-450, rr_file='${LIB}/blosum62.sim.mat',
          gap_penalties_1d=(-500, -50), n_prof_iterations=1,
          check_profile=False, max_aln_evalue=0.01)

# Write out the profile in text format
prf.write(file=ALIGNMENT_TABLE, profile_format='TEXT')

# Step 1: sequence search
#
# Search for sequences related to 1fdx.chn in the whole representative
# list of PDB chains and align them automatically with the target sequence.
# This will only work if pdb_95.pir is first downloaded from the Modeller
# website into this directory.
# The output is a file search.prf listing potential templates.

from modeller import *
import sys

log.verbose()
env = environ()

# Read in the database of PDB chains clustered at 95% sequence identity
sdb = sequence_db(env)
try:
    sdb.read(seq_database_file='pdb_95.pir', seq_database_format='PIR',
             chains_list='ALL', minmax_db_seq_len=(30, 4000),
             clean_sequences=True)
except IOError:
    errmsg = str(sys.exc_info()[1]) + """
   Could not open the PDB 95 database. This database is *not* included in the
   Modeller release, since it is frequently updated. To run this example, you
   first need to download pdb_95.pir.gz from the "Data file downloads" page on
   the Modeller website. (However, you can still run the rest of the
   examples in this directory without that database.)
"""
    print(errmsg)
    sys.exit(0)

# Read in the target sequence in PIR alignment format
aln = alignment(env)
aln.append(file='my.pir', alignment_format='PIR', align_codes='ALL')

# Convert the target sequence from alignment to profile format
prf = aln.to_profile()

# Scan sequence database to pick up homologous sequences
prf.build(sdb, matrix_offset=-450, rr_file='${LIB}/blosum62.sim.mat',
          gap_penalties_1d=(-500, -50), n_prof_iterations=1,
          check_profile=False, max_aln_evalue=0.01)

# Write out the profile in text format
prf.write(file='search.prf', profile_format='TEXT')

"""
Positional arguments:
1. INPUT - PIR database 
2. OUTPUT - BINARY database
"""
import sys
PIR = sys.argv[1]
BIN = sys.argv[2]

import os
from modeller import *

log.verbose()
env = environ()

#-- Read in the sequence database
sdb = sequence_db(env)
sdb.read(seq_database_file=PIR, seq_database_format='PIR',
         chains_list='ALL', minmax_db_seq_len=(30, 4000), clean_sequences=True)

#-- Write the sequence database in binary form
sdb.write(seq_database_file=BIN, seq_database_format='BINARY',
          chains_list='ALL')

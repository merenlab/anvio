# Step 5: adds chain identifier to a single model structure
import sys
from modeller import *

ENV = sys.argv[1]
NAME = sys.argv[2]

# load the optimized model and add chain ID, overwrite
env = environ(rand_seed=-12312)
env.io.atom_files_directory = [ENV]
mdl = model(env, file=NAME)
mdl.rename_segments(segment_ids=("A"))
mdl.write(file=NAME)

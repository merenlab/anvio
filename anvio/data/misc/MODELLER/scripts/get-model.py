# Step 4: model building
import sys
ALIGNMENT          = sys.argv[1]
TARGET_ID          = sys.argv[2]
TEMPLATE_IDS_FILE  = sys.argv[3]
NUM_MODELS         = int(sys.argv[4])
DEVIATION          = float(sys.argv[5])
VERY_FAST          = int(sys.argv[6])

# read in file written by MODELLER.run_align_to_templates(); creates tuple, e.g (1durA, 2fdnB, ...)
template_ids = tuple(["".join(x.strip().split("\t")) for x in open(TEMPLATE_IDS_FILE).readlines()])

from modeller import *
from modeller.automodel import *    # Load the automodel class
import os

log.verbose()
env = environ(rand_seed=-12312)

# directories for input atom files
env.io.atom_files_directory = ['./%s_TEMPLATE_PDBS' % TARGET_ID]

a = automodel(env,
              alnfile=ALIGNMENT,            # alignment filename
              knowns=template_ids,          # codes of the templates
              sequence=TARGET_ID,           # code of the target
              assess_methods=assess.GA341)  # request GA341 assessment

# prepare for an extremely fast optimization
if VERY_FAST:
    a.very_fast()

# index of the first model
a.starting_model = 1
# index of the last model (i.e. the total number of models)
a.ending_model = NUM_MODELS
# has to >0 if more than 1 model
a.deviation = DEVIATION

# do the homology modelling
a.make()

# Create the average model (unoptimized: cluster.ini, optimized: cluster.opt).
a.cluster()

# rename the model structures. we do this here instead of in the anvi'o driver because we have
# direct access to the file names (they don't have to be guessed).  Model files are renamed here
# according to a logic that is assumed by the anvi'o driver. Mess this with and you mess with
# anvi'o.
MODEL_PDBS = "%s_MODEL_PDBS" % TARGET_ID
try:
    os.mkdir(MODEL_PDBS)
except FileExistsError:
    pass

models = a.outputs
for ind, model in enumerate(models):
    os.rename(model['name'], os.path.join(MODEL_PDBS, "%s_model_%d.pdb" % (TARGET_ID, ind+1)))

os.rename("cluster.opt", os.path.join(MODEL_PDBS, TARGET_ID+"_model_avg.pdb"))

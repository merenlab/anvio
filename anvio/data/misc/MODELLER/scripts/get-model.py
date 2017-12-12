# Step 4: model building
import sys
ALIGNMENT          = sys.argv[1]
TARGET_ID          = sys.argv[2]
TEMPLATE_IDS_FILE  = sys.argv[3]
NUM_MODELS         = int(sys.argv[4])
DEVIATION          = float(sys.argv[5])

# read in file written by MODELLER.run_align_to_templates()
# creates tuple, e.g (1durA, 2fdnB, ...)
template_ids = tuple(["".join(x.strip().split("\t")) for x in open(TEMPLATE_IDS_FILE).readlines()])

from modeller import *
from modeller.automodel import *    # Load the automodel class

log.verbose()
env = environ(rand_seed=-12312)  # To get different models from another script

# directories for input atom files
env.io.atom_files_directory = ['./TEMPLATE_PDBS']

a = automodel(env,
              alnfile=ALIGNMENT,            # alignment filename
              knowns=template_ids,          # codes of the templates
              sequence=TARGET_ID,           # code of the target
              assess_methods=assess.GA341)  # request GA341 assessment
a.starting_model= 1                 # index of the first model
a.ending_model  = NUM_MODELS        # index of the last model
                                    # (determines how many models to calculate)
a.deviation = DEVIATION             # has to >0 if more than 1 model

a.make()                            # do homology modelling

# Step 4: model building
#
# This script should produce two models, 1fdx_my.B99990001.pdb and
# 1fdx_my.B99990002.pdb.

from modeller import *
from modeller.automodel import *    # Load the automodel class

log.verbose()
env = environ(rand_seed=-12312)  # To get different models from another script

# directories for input atom files
env.io.atom_files_directory = ['../atom_files']

a = automodel(env,
              alnfile='alignment.ali',      # alignment filename
              knowns=('1durA', '2fdnA'),    # codes of the templates
              sequence='1fdx_my',           # code of the target
              assess_methods=assess.GA341)  # request GA341 assessment
a.starting_model= 1                 # index of the first model
a.ending_model  = 2                 # index of the last model
                                    # (determines how many models to calculate)
a.deviation = 4.0                   # has to >0 if more than 1 model

a.make()                            # do homology modelling

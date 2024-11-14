# Step 4: model building
import sys

ALIGNMENT = sys.argv[1]
TARGET_ID = sys.argv[2]
TEMPLATE_IDS_FILE = sys.argv[3]
NUM_MODELS = int(sys.argv[4])
DEVIATION = float(sys.argv[5])
VERY_FAST = int(sys.argv[6])
MODEL_INFO = sys.argv[7]

# read in file written by MODELLER.run_align_to_templates(); creates tuple, e.g (1durA, 2fdnB, ...)
template_ids = tuple(
    ["".join(x.strip().split("\t")) for x in open(TEMPLATE_IDS_FILE).readlines()]
)

from modeller import *
from modeller.automodel import *  # Load the automodel class

log.verbose()
env = environ(rand_seed=-12312)

# directories for input atom files
env.io.atom_files_directory = ["./%s_TEMPLATE_PDBS" % TARGET_ID]

a = automodel(
    env,
    alnfile=ALIGNMENT,  # alignment filename
    knowns=template_ids,  # codes of the templates
    sequence=TARGET_ID,  # code of the target
    assess_methods=(assess.GA341, assess.DOPE),
)

a.set_output_model_format("PDB")

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

###################################################################################################

"""
`a` has the important attribute `outputs` that is a list of dictionaries giving
information about each model. I write this information into a text-delimited file
for further analysis with anvi'o
"""

models = a.outputs

model_info = {}
ignore = ["pdfterms", "failure"]
for model in models:
    for key in model:
        if key in ignore:
            continue
        if key not in model_info.keys():
            model_info[key] = []
        if key == "GA341 score":
            model_info[key].append(model[key][0])
        else:
            model_info[key].append(model[key])

model_info["num"].extend([NUM_MODELS])
model_info["molpdf"].extend([""])
model_info["GA341 score"].extend([""])
model_info["DOPE score"].extend([""])

columns = ["num", "name", "molpdf", "GA341 score", "DOPE score"]
f = open(MODEL_INFO, "w")
f.write("\t".join([column.replace(" ", "_") for column in columns]) + "\n")

for i in range(NUM_MODELS):
    line = []
    for column in columns:
        line.append(str(model_info[column][i]))
    f.write("\t".join(line) + "\n")
f.close()

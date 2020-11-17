# Load libraries
#---------------
# packages <- c("tidyverse", "ggpubr", "odseq")
library("tidyverse")
# suppressMessages(lapply(packages, library, character.only = TRUE))

# Set WD
#-------
setwd("~/github/ACE_SO_metagenomes/analyses/SCGs_test")

# Get 
args <- commandArgs()
cat(args, sep = "\n")

# Load MSA
#---------



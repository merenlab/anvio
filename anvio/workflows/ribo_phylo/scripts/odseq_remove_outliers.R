#! ~/miniconda3/envs/anvio-master/bin/r Rscript

# Load libraries
#---------------
packages <- c("tidyverse", "odseq", "optparse", "msa")
suppressMessages(lapply(packages, library, character.only = TRUE))

# Set WD
#-------
setwd("~/github/ACE_SO_metagenomes/analyses/SCGs_test")

# Parse command line arugments
option_list <- list(
  make_option(
    c("-msa"),
    type="character",
    help="help"
  )
)

op <- OptionParser(
  option_list=option_list,
  description="Use OD-SEQ to remove outliers from MSAs",
)


args <- parse_args(op)

# Load MSA
#---------

#! ~/miniconda3/envs/anvio-master/bin/r Rscript

# Load libraries
#---------------
packages <- c("dplyr", "odseq", "optparse", "msa", "readr")
suppressMessages(lapply(packages, library, character.only = TRUE))

# Parse arguments
#----------------
option_list <- list(
  make_option(
    c("-m", "--msa"),
    type="character",
    help="help"
  ),
  make_option(
    c("-o", "--output"),
    type="character",
    help="help"
  )
)

op <- OptionParser(
  option_list=option_list,
  description="Here is what the program does",
)
args <- parse_args(op)


# Load fasta
args$msa <- "~/github/ACE_SO_metagenomes/analyses/SCGs_test/05_MSA/MSA_1/Ribosomal_L16/Ribosomal_L16_aligned.fasta"
msa_all <- readAAStringSet(args$msa)

# Align with Muscle
msa_all_muscle <- msa(msa_all, method = "Muscle")

msa_outliers_removed <- odseq(msa_all_muscle)
msa_outliers_removed

args$output


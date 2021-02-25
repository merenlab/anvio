library(ape)
library(fs)
library(tidyverse)
library(phylotools)


packages <- c("tidyverse","ape")
suppressMessages(lapply(packages, library, character.only = TRUE))

setwd("~/github/ACE_SO_metagenomes/analyses/test_case")
# Make index of names and splits 
#-------------------------------
#####
SCG = "Ribosomal_L16"
####

# make paths
#-----------
# misc data
misc_data_suffix <- stringr::str_c(SCG, "_all_misc_data_final.tsv")
misc_data_path <- path("RIBO_PHYLO_WORKFLOW/07_MISC_DATA/", SCG, misc_data_suffix)

# reformat file
reformat_file_suffix <- stringr::str_c(SCG, "-reformat-report.txt")
reformat_file_path <- path("PROFILE_SCGs/02_FASTA/", SCG, reformat_file_suffix)

# gene_calls
gene_calls_path <- path("PROFILE_SCGs/07_SUMMARY/", SCG, "bin_by_bin/EVERYTHING/EVERYTHING-gene_calls.txt")

# Mapping data
mapping_data_path <- path("PROFILE_SCGs/07_SUMMARY/", SCG, "bin_by_bin/EVERYTHING/EVERYTHING-gene_coverages.txt")

# Output file
outfile_path_suffix <- stringr::str_c(SCG, "_all_misc_data_final_coverage.tsv")
outfile_path <- path("RIBO_PHYLO_WORKFLOW/07_MISC_DATA/", SCG, outfile_path_suffix)


# load data
#----------
misc_data <- read_tsv(misc_data_path)
reformat_file <- read_tsv(reformat_file_path, col_names = FALSE) %>% rename(contig = X1, name = X2)
gene_calls <- read_tsv(gene_calls_path)
mapping_data <- read_tsv(mapping_data_path)
# Load in tree and splits data
tree <- read.tree("RIBO_PHYLO_WORKFLOW/06_TREES/Ribosomal_L16/Ribosomal_L16.nwk")

# Make index of names
#--------------------
index <- mapping_data %>%
  left_join(gene_calls) %>%
  left_join(reformat_file) %>%
  select(gene_callers_id, contig, name)

# convert
index$contig_new <- paste(index$contig, "_split_00001", sep="")

# replace tips

dat <- index %>% select(name, contig_new) %>% as.data.frame()
ntree <- sub.taxa.label(tree, dat) # Bless this persons heart: https://rdrr.io/cran/phylotools/man/sub.taxa.label.html


# Export
write.tree(ntree, file = "RIBO_PHYLO_WORKFLOW/06_TREES/Ribosomal_L16/Ribosomal_L16_test.nwk")

# gold: https://tbradley1013.github.io/2018/06/19/subsetting-phylogenetic-trees/

# Now change the misc data to have the same names as the splits
misc <- read_tsv("RIBO_PHYLO_WORKFLOW/07_MISC_DATA/Ribosomal_L16/Ribosomal_L16_all_misc_data_final_coverage_project.tsv")

index_tibble <- index %>% select(name, contig_new)

misc_new <- misc %>% left_join(index_tibble) %>% select(-name) %>% rename(name = contig_new) %>% relocate(name)

misc_new %>% write_tsv("RIBO_PHYLO_WORKFLOW/07_MISC_DATA/Ribosomal_L16/Ribosomal_L16_all_misc_data_final_coverage_project_final.tsv")

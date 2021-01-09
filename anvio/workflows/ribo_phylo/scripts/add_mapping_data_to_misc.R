# Load libraries
#---------------
packages <- c("tidyverse")
suppressMessages(lapply(packages, library, character.only = TRUE))

setwd("~/scratch/OSD")

# misc_data <- read_tsv("RIBO_PHYLO_WORKFLOW/07_MISC_DATA/Ribosomal_L16/Ribosomal_L16_all_misc_data_final.tsv")
# mapping_data <- read_tsv("profile_SCGs/07_SUMMARY/Ribosomal_l16/bin_by_bin/EVERYTHING/EVERYTHING-gene_coverages.txt")
# reformat_file <- read_tsv("profile_SCGs/02_FASTA/Ribosomal_l16/Ribosomal_l16-reformat-report.txt", col_names = FALSE) %>%
#   rename(contig = X1, name = X2)
# 
# misc_data <- misc_data %>%
#   filter(new_header != "new_header")
# 
# index <- external_gene_calls %>%
#   select(gene_callers_id, contig) %>%
#   left_join(mapping_data) %>%
#   left_join(reformat_file)
# 
# misc_plus_cov <- external_gene_calls %>%
#   select(gene_callers_id, contig) %>%
#   left_join(index) %>%
#   select(-contig, -gene_callers_id)
# 
# misc_plus_cov_2 <- misc_data %>% left_join(misc_plus_cov)

# For now I am hacking this and just col_binding the two datasets. 
misc_data <- read_tsv("RIBO_PHYLO_WORKFLOW/07_MISC_DATA/Ribosomal_L16/Ribosomal_L16_all_misc_data_final.tsv")
mapping_data <- read_tsv("profile_SCGs/07_SUMMARY/Ribosomal_l16/bin_by_bin/EVERYTHING/EVERYTHING-gene_coverages.txt")

misc_data <- misc_data %>%
  filter(new_header != "new_header")

misc_plus_cov_2 <- misc_data %>%
  bind_cols(mapping_data) %>%
  select(-gene_callers_id)

write_tsv(misc_plus_cov_2, "08_MISC_DATA/Ribosomal_L16/Ribosomal_L16_all_misc_data_final_coverage.tsv")



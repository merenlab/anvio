# Load libraries
#---------------
packages <- c("tidyverse")
suppressMessages(lapply(packages, library, character.only = TRUE))

misc_data <- read_tsv("08_MISC_DATA/Ribosomal_L16/Ribosomal_L16_all_misc_data_final.tsv")
mapping_data <- read_tsv("profile_SCGs/07_SUMMARY/g01/bin_by_bin/EVERYTHING/EVERYTHING-gene_coverages.txt")
external_gene_calls <- read_tsv("profile_SCGs/02_FASTA/g01/g01-external-gene-calls.txt")
reformat_file <- read_tsv("profile_SCGs/02_FASTA/g01/g01-reformat-report.txt", col_names = FALSE) %>%
  rename(contig = X1, name = X2)

index <- external_gene_calls %>% 
  select(gene_callers_id, contig) %>%
  left_join(mapping_data) %>%
  left_join(reformat_file)
  
misc_plus_cov <- external_gene_calls %>%
  select(gene_callers_id, contig) %>%
  left_join(index) %>%
  select(-contig, -gene_callers_id) 

misc_plus_cov_2 <- misc_data %>% left_join(misc_plus_cov)

write_tsv(misc_plus_cov_2, "08_MISC_DATA/Ribosomal_L16/Ribosomal_L16_all_misc_data_final_coverage.tsv")

# Load libraries
#---------------
packages <- c("tidyverse")
suppressMessages(lapply(packages, library, character.only = TRUE))

setwd("/Users/mschechter/github/ACE_SO_metagenomes/analyses/test_case")

misc_data <- read_tsv("RIBO_PHYLO_WORKFLOW/07_MISC_DATA/Ribosomal_L16/Ribosomal_L16_all_misc_data_final.tsv")
mapping_data <- read_tsv("profile_SCGs/07_SUMMARY/Ribosomal_l16/bin_by_bin/EVERYTHING/EVERYTHING-gene_coverages.txt")
reformat_file <- read_tsv("profile_SCGs/02_FASTA/Ribosomal_l16/Ribosomal_l16-reformat-report.txt", col_names = FALSE) %>%
  rename(contig = X1, name = X2)
gene_calls <- read_tsv("profile_SCGs/07_SUMMARY/Ribosomal_l16/bin_by_bin/EVERYTHING/EVERYTHING-gene_calls.txt")

# Make index of names
index <- mapping_data %>% left_join(gene_calls) %>% left_join(reformat_file) %>% select(gene_callers_id, contig, name)

# Make final table
final <- mapping_data %>% 
          left_join(index) %>%
          left_join(misc_data %>% rename(name = new_header)) %>%
          select(-gene_callers_id, -contig) %>%
          select(name, sample, contig_db_type, t_domain, t_phylum,t_order    ,    t_family    ,      t_genus ,   t_species , has_genomic_SCG_in_cluster , ERR771092)
write_tsv(final, "RIBO_PHYLO_WORKFLOW/07_MISC_DATA/Ribosomal_L16/Ribosomal_L16_all_misc_data_final_coverage.tsv")




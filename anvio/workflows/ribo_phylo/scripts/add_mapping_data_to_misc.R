# Load libraries
#---------------
packages <- c("tidyverse", "fs")
suppressMessages(lapply(packages, library, character.only = TRUE))

add_coverage_to_metadata <- function(SCG) {
  
  #####
  # SCG = "Ribosomal_l17"
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
  
  # Make index of names
  #--------------------
  index <- mapping_data %>% 
    left_join(gene_calls) %>% 
    left_join(reformat_file) %>% 
    select(gene_callers_id, contig, name)

  # Make final table
  misc_data %>%
    rename(name = new_header) %>%
    left_join(index) %>%
    inner_join(mapping_data) %>%
    select(-gene_callers_id, -contig) %>%
    write_tsv(outfile_path)
}

add_coverage_to_metadata(SCG = "Ribosomal_L17")
add_coverage_to_metadata(SCG = "Ribosomal_L16")

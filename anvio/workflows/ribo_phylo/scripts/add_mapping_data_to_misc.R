# Load libraries
#---------------
packages <- c("tidyverse", "fs", "optparse")
suppressMessages(lapply(packages, library, character.only = TRUE))

option_list = list(
  make_option(c("-s", "--SCG"), action="store", default=NA, type='character',
              help="name of scg"),
  make_option(c("--ribophylopath"), action="store", default=NA, type='character',
              help="choose 'Num_references' or 'clustering_threshold'"),
  make_option(c("--profilepath"), action="store", default=NA, type='character',
              help="input data"),
  make_option(c("--genecalls"), action="store", default=NA, type='character',
              help="input data"),
  make_option(c("-o", "--output-file"), action="store", default=NA, type='character',
              help="output file.") 
)

opt <- parse_args(OptionParser(option_list=option_list))


# For dev purposes
#-----------------
opt$SCG <- "Ribosomal_S2"
opt$ribophylopath <- "RIBO_PHYLO_WORKFLOW_94"
opt$profilepath <- "PROFILE_SCGs_94"
#----------------



add_coverage_to_metadata <- function(SCG) {
  
  #####
  # SCG = "Ribosomal_L16"
  ####
  
  # make paths
  #-----------
  # misc data
  misc_data_suffix <- stringr::str_c(opt$SCG, "_all_misc_data_final.tsv")
  misc_data_path <- path(opt$ribophylopath, "07_MISC_DATA", opt$SCG , misc_data_suffix)
  
  # reformat file
  reformat_file_suffix <- stringr::str_c(opt$SCG, "-reformat-report.txt")
  reformat_file_path <- path(opt$profilepath, "02_FASTA/", opt$SCG, reformat_file_suffix)
  
  # gene_calls
  gene_calls_path <- path(opt$profilepath, "07_SUMMARY", opt$SCG, "bin_by_bin/EVERYTHING/EVERYTHING-gene_calls.txt")
  
  # need to up load the tree toooooo
  
  # Mapping data
  mapping_data_path <- path(opt$profilepath, "07_SUMMARY/", opt$SCG, "bin_by_bin/EVERYTHING/EVERYTHING-gene_coverages.txt")
  
  # Output file
  outfile_path_suffix <- stringr::str_c(opt$SCG, "_all_misc_data_final_1.tsv")
  outfile_path <- path(opt$ribophylopath, "07_MISC_DATA/", opt$SCG, outfile_path_suffix)
  
  
  # load data
  #----------
  misc_data <- read_tsv(misc_data_path)
  reformat_file <- read_tsv(reformat_file_path, col_names = FALSE) %>% rename(contig = X1, name = X2)
  gene_calls <- read_tsv(gene_calls_path)
  # FIXME: gene_calls has some duplicates in the contig column IF some of the reference sequences when
  # being processed into contigDBs split them into two gene calls, you suck Prodigal! I could fix this
  # by providing an external-gene-calls file with the reference SCGs.
  # gene_calls <- gene_calls[!duplicated(gene_calls$contig), ]
  mapping_data <- read_tsv(mapping_data_path)
  
  # Make index of names
  #--------------------
  index <- mapping_data %>%
    inner_join(gene_calls) %>%
    left_join(reformat_file) %>%
    select(gene_callers_id, contig, name)
  
  # Make final table
  final <- misc_data %>%
    rename(name = new_header) %>%
    left_join(index) %>%
    # inner_join(mapping_data) %>%
    select(-gene_callers_id) %>%
    rename(orig_name = name, name = contig) %>%
    relocate(name) %>%
    filter(!is.na(name))
  
  final %>%
    write_tsv(outfile_path)
}

add_coverage_to_metadata(SCG = "Ribosomal_S2")
add_coverage_to_metadata(SCG = "Ribosomal_L16")

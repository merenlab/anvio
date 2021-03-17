packages <- c("tidyverse", "ape", "fs", "phylotools", "optparse")
suppressMessages(lapply(packages, library, character.only = TRUE))

option_list = list(
  make_option(c("-s", "--SCG"), action="store", default=NA, type='character',
              help="name of scg"),
  make_option(c("--ribophylopath"), action="store", default=NA, type='character',
              help="choose 'Num_references' or 'clustering_threshold'"),
  make_option(c("--profilepath"), action="store", default=NA, type='character',
              help="input data")
)

opt <- parse_args(OptionParser(option_list=option_list))

#####
# SCG = "Ribosomal_S2"
# ribophylopath = "RIBO_PHYLO_WORKFLOW_90"
# profilepath = "PROFILE_SCGs_competitive_90"
####

replace_tree_names <- function(SCG, ribophylopath, profilepath){
  # make paths
  #-----------
  # misc data
  misc_data_path_suffix <- stringr::str_c(SCG, "_all_misc_data_final.tsv")
  misc_data_path <- path(ribophylopath, "07_MISC_DATA/", SCG, misc_data_path_suffix)
  
  # reformat file
  reformat_file_suffix <- stringr::str_c(SCG, "-reformat-report.txt")
  reformat_file_path <- path(profilepath, "02_FASTA/", SCG, reformat_file_suffix)
  
  # gene_calls
  gene_calls_path <- path(profilepath, "07_SUMMARY/", SCG, "bin_by_bin/EVERYTHING/EVERYTHING-gene_calls.txt")
  
  # Mapping data
  mapping_data_path <-  path(profilepath, "07_SUMMARY/", SCG, "bin_by_bin/EVERYTHING/EVERYTHING-gene_coverages.txt")
  
  # tree path
  tree_path_suffix <- stringr::str_c(SCG, ".nwk")
  tree_path <-  path(ribophylopath, "06_TREES/", SCG, tree_path_suffix)
  
  # Output file
  outfile_path_suffix <- stringr::str_c(SCG, "_all_misc_data_final_coverage.tsv")
  outfile_path <- path(ribophylopath, "07_MISC_DATA/", SCG, outfile_path_suffix)
  
  # tree out path
  tree_path_out_suffix <- stringr::str_c(SCG, "_renamed.nwk")
  tree_path_out <- path(ribophylopath, "06_TREES/", SCG, tree_path_out_suffix)
  
  # # Output file
  # outfile_path_suffix <- stringr::str_c(opt$SCG, "_all_misc_data_final_1.tsv")
  # outfile_path <- path(opt$ribophylopath, "07_MISC_DATA/", opt$SCG, outfile_path_suffix)
  
  # Load data
  #----------
  misc_data <- read_tsv(misc_data_path)
  reformat_file <- read_tsv(reformat_file_path, col_names = FALSE) %>% rename(contig = X1, name = X2)
  gene_calls <- read_tsv(gene_calls_path)
  mapping_data <- read_tsv(mapping_data_path)
  tree <- read.tree(tree_path)
  
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
  ntree <- phylotools::sub.taxa.label(tree, dat) # Bless this persons heart: https://rdrr.io/cran/phylotools/man/sub.taxa.label.html
  
  # Export tree with new tip names
  #-------------------------------
  write.tree(ntree, file = tree_path_out)
  # gold: https://tbradley1013.github.io/2018/06/19/subsetting-phylogenetic-trees/
  
  # misc path
  #----------
  
  # convert
  index <- mapping_data %>%
    inner_join(gene_calls) %>%
    left_join(reformat_file) %>%
    select(gene_callers_id, contig, name)
  
  index$contig_new <- paste(index$contig, "_split_00001", sep="")
  
  # Make final table
  final <- misc_data %>%
    rename(name = new_header) %>%
    left_join(index) %>%
    # inner_join(mapping_data) %>%
    select(-gene_callers_id) %>%
    rename(orig_name = name, name = contig_new) %>%
    relocate(name) %>%
    filter(!is.na(name))
  
  final %>%
    write_tsv(outfile_path)
}

replace_tree_names(SCG = opt$SCG, ribophylopath = opt$ribophylopath, profilepath = opt$profilepath)


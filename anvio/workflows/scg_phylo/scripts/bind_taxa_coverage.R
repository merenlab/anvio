# This script binds the SCG taxonomy misc data with the SCG read recruitment data

packages <- c("tidyverse", "ape", "fs", "phylotools", "optparse")

new_packages <- packages[!(packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)

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

bind_taxa_coverage <- function(SCG, ribophylopath, profilepath){
  
  print(opt$SCG)
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
  
  # external-gene-calls
  external_gene_calls_path_suffix <- stringr::str_c(SCG, "_external_gene_calls_all_renamed.tsv")
  external_gene_calls_path <-  path(ribophylopath, "03_NR_FASTAS/", SCG, external_gene_calls_path_suffix)
  
  # Output file
  outfile_path_suffix <- stringr::str_c(SCG, "_all_misc_data_final_coverage.tsv")
  outfile_path <- path(ribophylopath, "07_MISC_DATA/", SCG, outfile_path_suffix)
  
  # tree out path
  tree_path_out_suffix <- stringr::str_c(SCG, "_renamed.nwk")
  tree_path_out <- path(ribophylopath, "06_TREES/", SCG, tree_path_out_suffix)
  
  # Load data
  #----------
  
  ######
  # misc_data_path <- "RIBO_PHYLO_WORKFLOW/07_MISC_DATA/Ribosomal_L16/Ribosomal_L16_all_misc_data_final.tsv"
  # reformat_file_path <- "METAGENOMICS_WORKFLOW/02_FASTA/Ribosomal_L16/Ribosomal_L16-reformat-report.txt"
  # gene_calls_path <- "METAGENOMICS_WORKFLOW/07_SUMMARY/Ribosomal_L16/bin_by_bin/EVERYTHING/EVERYTHING-gene_calls.txt"
  # mapping_data_path <- "METAGENOMICS_WORKFLOW/07_SUMMARY/Ribosomal_L16/bin_by_bin/EVERYTHING/EVERYTHING-gene_coverages.txt"
  # tree_path <- "RIBO_PHYLO_WORKFLOW/06_TREES/Ribosomal_L16/Ribosomal_L16.nwk"
  # external_gene_calls_path <- "RIBO_PHYLO_WORKFLOW/03_NR_FASTAS/Ribosomal_L16/Ribosomal_L16_external_gene_calls_all_renamed.tsv"
  #####
  
  misc_data <- read_tsv(misc_data_path)
  gene_calls <- read_tsv(gene_calls_path)
  mapping_data <- read_tsv(mapping_data_path)
  tree <- read.tree(tree_path)
  external_gene_calls <- read_tsv(external_gene_calls_path)
  
  # Make index of names
  #--------------------
  index <- mapping_data %>%
    # left_join(gene_calls) %>%
    # left_join(reformat_file) %>%
    left_join(external_gene_calls) %>%
    select(gene_callers_id, contig)
  
  # convert
  index$contig_new <- paste(index$contig, "_split_00001", sep="")
  
  # replace tips
  dat <- index %>% select(contig, contig_new) %>% filter(contig %in% tree$tip.label)%>% as.data.frame()
  ntree <- phylotools::sub.taxa.label(tree, dat) # Bless this persons heart: https://rdrr.io/cran/phylotools/man/sub.taxa.label.html
  
  # Export tree with new tip names
  #-------------------------------
  write.tree(ntree, file = tree_path_out)
  # gold: https://tbradley1013.github.io/2018/06/19/subsetting-phylogenetic-trees/
  
  # # Make final table
  final <- misc_data %>%
    rename(contig = new_header) %>%
    left_join(index) %>%
    # inner_join(mapping_data) %>%
    select(-gene_callers_id) %>%
    rename(orig_name = contig, name = contig_new) %>%
    relocate(name) %>%
    filter(!is.na(name))

  final %>%
    write_tsv(outfile_path)
}

bind_taxa_coverage(SCG = opt$SCG, ribophylopath = opt$ribophylopath, profilepath = opt$profilepath)


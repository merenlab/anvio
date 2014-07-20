#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(ape))
suppressPackageStartupMessages(library(ade4))
suppressPackageStartupMessages(library(vegan))
suppressPackageStartupMessages(library(gtools))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(flashClust))

# command line options
option_list <- list(
		make_option(c("-d", "--distance"), default="euclidean",
				help = 'Dissimilarity index for clustering of contigs based on their TNF profiles.
						Default is "%default". Available indices are "manhattan", "euclidean", "canberra",
						"bray", "kulczynski", "jaccard", "gower", "altGower", "morisita", "horn",
						"mountford", "raup" , "binomial", "chao" and "cao"'),
		make_option(c("-m", "--method"), default="ward",
				help = 'Agglomeration method for clustering. Default is "%default". Available
						methods are "ward", "single", "complete", "average", "mcquitty", "median"
						and "centroid"'),
		make_option(c("-o", "--output_file_prefix"), default="unknown-tree.txt",
				help = "Output file [default \"%default\"]")
)

parser <- OptionParser(usage = "TNF_matrix_to_newick.R [options] input_matrix", option_list=option_list,
		description="A script to generate newick formattet trees for clustering of contigs")

arguments <- parse_args(parser, positional_arguments = TRUE)
options <- arguments$options

if(invalid(options$metadata))
	options$metadata <- NA

# check if the positional argument is set
if(length(arguments$args) != 1) {
	cat("Incorrect number of required positional arguments\n\n")
	print_help(parser)
	stop()
} else {
	matrix_file <- arguments$args
}

raw_data <- data.frame(t(data.matrix(read.table(matrix_file, header = TRUE, row.names = 1,sep="\t", check.names = FALSE))), check.names = FALSE)

# generate a data frame to hold old and new contig names.. this is just because that
# fucking hclust2phylog function does not have a check.names parameter. all the shitty
# code that follows could have been avoided otherwise. this will cause a lot of
# performance issues as well.
contig_names <- names(raw_data) 
conversion <- data.frame(old = character(), new = character(), stringsAsFactors=FALSE)
N <- 0
for(contig_name in contig_names){
    N = N + 1
    conversion[N, ] <- c(contig_name, paste('contig', sprintf("%09d", N), sep="_"))
}

# put temoprary names in
names(raw_data) <- conversion$new
raw_data <- data.matrix(raw_data)

dcols<-vegdist(t(raw_data), method=options$distance, na.rm=TRUE)
hc <- hclust(dcols, method=options$method)
phy <- hclust2phylog(hc, add.tools = TRUE)
write.tree(as.phylo(phy), file=options$output_file_prefix)

# now read that tree you just wrote
output_tree <- readChar(options$output_file_prefix, file.info(options$output_file_prefix)$size)

# convert temporary names into originals
for(contig_name in contig_names){
    output_tree <- gsub(conversion[conversion$old == contig_name, ]$new, contig_name, output_tree)
}

updated_output_path <- file(options$output_file_prefix)
writeLines(output_tree, updated_output_path)
close(updated_output_path)

#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(ape))
suppressPackageStartupMessages(library(ade4))
suppressPackageStartupMessages(library(vegan))
suppressPackageStartupMessages(library(gtools))
suppressPackageStartupMessages(library(optparse))

# command line options
option_list <- list(
		make_option(c("-d", "--distance"), default="horn",
				help = 'Dissimilarity index for clustering of contigs based on their TNF profiles.
						Default is "%default". Available indices are "manhattan", "euclidean", "canberra",
						"bray", "kulczynski", "jaccard", "gower", "altGower", "morisita", "horn",
						"mountford", "raup" , "binomial", "chao" and "cao"'),
		make_option(c("-m", "--method"), default="ward",
				help = 'Agglomeration method for clustering. Default is "%default". Available
						methods are "ward", "single", "complete", "average", "mcquitty", "median"
						and "centroid"'),
		make_option(c("-o", "--output_file_prefix"), default="unknown",
				help = "Output file prefix for output files [default \"%default\"]")
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

raw_data <- t(data.matrix(read.table(matrix_file, header = TRUE, row.names = 1,sep="\t")))
dcols<-vegdist(t(raw_data), method=options$distance, na.rm=TRUE)
hc <- hclust(dcols, method=options$method)
phy <- hclust2phylog(hc, add.tools = TRUE)
write.tree(as.phylo(phy), file=paste(options$output_file_prefix, '.txt', sep=''))
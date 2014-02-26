#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(optparse))

option_list <- list()

parser <- OptionParser(usage = "script.R [options] reports_file", option_list=option_list,
        description="Visualize coverage and entropy for contigs")

arguments <- parse_args(parser, positional_arguments = TRUE)
options <- arguments$options

if(length(arguments$args) != 1) {
    cat("Incorrect number of required positional arguments\n\n")
    print_help(parser)
    stop()
} else {
    input_file_path <- arguments$args[1]
}

if(file.access(input_file_path) == -1){
    stop(sprintf("Input file '%s' does not exist", input_file_path))
}


df <- data.frame(read.table(input_file_path, header = TRUE, sep="\t"))

for(contig in levels(df$contig)){
    dfx <- df[df$contig == contig, ]
    p <- ggplot(data=dfx, aes(x=pos, y=entropy)) + geom_point()
    r <- ggplot(data=dfx, aes(x=pos, y=entropy_n)) + geom_point()
    q <- ggplot(data=dfx, aes(x=pos,y=coverage)) + geom_line(stat="identity")# + scale_y_sqrt()
    png(paste(input_file_path, '-', contig, '.png', sep=""), width = (max(dfx$pos) / 20000) * 200, height = 8 * 100)
    grid.arrange(p, r, q)
    dev.off()
}

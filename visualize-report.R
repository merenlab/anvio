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
	print(contig)
    df_contig <- df[df$contig == contig, ]
    dfx <- df_contig[df_contig$entropy_n > 0, ]

    p <- ggplot(data=dfx, aes(x=pos, y=entropy_n, color=competing_nt, group=competing_nt, fill=competing_nt))
    p <- p + geom_bar(stat="identity") + geom_point(aes(size=coverage))
    p <- p + geom_line(size=0.2)
    p <- p + facet_grid(competing_nt ~ .)
    p <- p + theme(legend.position = 'top', strip.text.y = element_blank(), strip.background = element_blank(), plot.margin = unit(c(0,0,0,0) , units = "lines"))
    p <- p + coord_cartesian(xlim = c(-max(df_contig$pos)/500, max(df_contig$pos) + max(df_contig$pos)/500), ylim=c(0,1))
    p <- p + scale_y_continuous(labels=c("", "0.25", "0.50", "0.75", "")) 
    p <- p + scale_x_continuous(labels=c()) 
    p <- p + theme(axis.ticks.y=element_blank(), axis.ticks.x=element_blank())
	p <- p + xlab("") 


    q <- ggplot(data=df_contig, aes(x=pos,y=coverage))
    q <- q + geom_line()
    q <- q + coord_cartesian(xlim = c(-max(df_contig$pos)/500, max(df_contig$pos) + max(df_contig$pos)/500))
    q <- q + theme(plot.margin = unit(c(0,0,0,0), units = "lines"))
    q <- q + theme(axis.ticks.x=element_blank())

    gA <- ggplotGrob(p)
    gB <- ggplotGrob(q)
    maxWidth = grid::unit.pmax(gA$widths[2:5], gB$widths[2:5])
    gA$widths[2:5] <- as.list(maxWidth)
    gB$widths[2:5] <- as.list(maxWidth)


    png(paste(input_file_path, '-', contig, '.png', sep=""), width = (max(dfx$pos) / 10000) * 200, height = 8 * 100)
    grid.arrange(gA, gB, heights=c(7,2))
    dev.off()
}



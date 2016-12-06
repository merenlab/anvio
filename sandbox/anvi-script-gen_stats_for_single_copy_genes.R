#!/usr/bin/env Rscript
#
# visualizes stuff..
#

suppressPackageStartupMessages(library(gtools))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(optparse))

# command line options
option_list <- list(
        make_option(c("--e_value"), default=1e-0,
                help = "e-value to retain hits [default \"%default\"]"),
        make_option(c("--output_prefix"),
                help = "Output file name *prefix* (without the extension)"),
        make_option(c("--alphabetical_order"), action = "store_true", default = FALSE,
                help = "Order genes based on name, instead of frequency"),
        make_option(c("--no_colors"), action = "store_true", default = FALSE,
                help = "Order genes based on name, instead of frequency")
        )

parser <- OptionParser(usage = "./self hits_file genes_file", option_list=option_list,
        description="Visualize stuff..")

arguments <- parse_args(parser, positional_arguments = TRUE)
options <- arguments$options

e_value <- options$e_value

# check if the positional argument is set
if(length(arguments$args) != 2) {
    cat("Incorrect number of required positional arguments\n\n")
    print_help(parser)
    stop()
} else {
    hits <- arguments$args[1]
    genes <- arguments$args[2]
}

if(invalid(options$output_prefix)){
    output_prefix <- paste(hits, '_e_', e_value, '_', sep='')
} else {
    output_prefix <- paste(options$output_prefix, '_e_', e_value, '_', sep='')
}

if(options$alphabetical_order){
    order_genes_by_frequency <- FALSE
} else {
    order_genes_by_frequency <- TRUE
}

if(options$no_colors){
    no_colors <- TRUE
} else {
    no_colors <- FALSE
}


# check if the input file is accessible
if(file.access(hits) == -1){
    stop(sprintf("Specified file '%s' does not exist", hits))
}

#e_value <- 1e-20
#hits <- '~/Infant-MERGED-FINAL-5K-FOR-THE-PAPER/Infant-gut-BINS/E_faecaelis-REF-OG1RF.db.hits'
#genes <- '~/Infant-MERGED-FINAL-5K-FOR-THE-PAPER/Infant-gut-BINS/E_faecaelis-REF-OG1RF.db.genes'
# source <- 'Wu_et_al'
# source <- 'Dupont_et_al'

genes_df <- data.frame(read.table(genes, header = TRUE,sep="\t"))
hits_df <- data.frame(read.table(hits, header = TRUE,sep="\t"))
hits_df <- hits_df[hits_df$e_value < e_value, ]

plots <- list()  # new empty list
i <- 1
for(source in c('Alneberg_et_al', 'Creevey_et_al', 'Campbell_et_al', 'Dupont_et_al')){
	print(source)
    source_genes_df <- genes_df[genes_df$source == source, ]
    source_genes_df$source <- factor(source_genes_df$source)
    source_genes_df$gene <- factor(source_genes_df$gene)

    source_df <- hits_df[hits_df$source == source, ]
    source_df$source <- factor(source_df$source)
    source_df$gene <- factor(source_df$gene)

    x <- data.frame(gene = character(), count = numeric(), stringsAsFactors=FALSE)
    N <- 1
    for(gene in levels(source_genes_df$gene)){
        x[N, ] <- c(gene, nrow(source_df[source_df$gene == gene, ]))
        N <- N + 1
    }

    x$count <- as.numeric(x$count)

    frequencies <- as.data.frame(table(x$count))
    names(frequencies) <- c('num', 'freq')
    frequencies <- frequencies[with(frequencies, order(-freq)), ]

    num_genomes <- as.character(frequencies[1, ]$num)
    percent_agrees <- frequencies[1, ]$freq * 100 / sum(frequencies$freq)

    text = sprintf("%s:\n%.2f%% of %d genes\noccur %s times", source, percent_agrees, nrow(x), num_genomes)
    p <- ggplot() + annotate("text", x = 1, y = 1, size=7, label = text) + theme(line = element_blank(),
            text = element_blank(),
            title = element_blank())

    q <- ggplot(x, aes(x=factor(0), y=count))
    q <- q + geom_violin()
    q <- q + geom_jitter(position = position_jitter(width = .2, height = 0), alpha=0.3)
    q <- q + theme_bw()
    q <- q + theme(axis.text.x = element_blank())
    q <- q + theme(axis.text.y = element_text(size = 12))
    q <- q + theme(axis.ticks.x = element_blank())
    q <- q + coord_cartesian(ylim = c(0, ifelse (max(x$count) > 0, max(x$count), 1)))
    q <- q + labs(x='', y='')

    # Draw bars
    if(order_genes_by_frequency)
        r <- ggplot(x, aes(x=reorder(gene, -count), y=as.integer(count)))
    else
        r <- ggplot(x, aes(x=gene, y=as.integer(count)))

    if(no_colors)
        r <- r + geom_bar(stat='identity', width=.9)
    else
        r <- r + geom_bar(aes(fill=count), stat='identity', width=.9)

    r <- r + theme_bw()
    r <- r + expand_limits(y=0)
    r <- r + coord_cartesian(ylim = c(0, ifelse (max(x$count) > 0, max(x$count), 1)))
    r <- r + theme(axis.text.y = element_text(size = 12))
    r <- r + theme(axis.text.x = element_blank())
    r <- r + geom_text(aes(x=gene, y=ifelse (max(x$count) > 0, max(x$count) * 0.001, 0.001), hjust = 0, vjust=0.5, label=gene), angle=90, size = 2.85)
    r <- r + labs(x='', y='')
    r <- r + scale_y_sqrt()

    plots[[i]] <- p
    i <- i + 1
    plots[[i]] <- q
    i <- i + 1
    plots[[i]] <- r
    i <- i + 1

}
pdf(paste(output_prefix, 'new.pdf', sep=''), width=25, height=12)
grid.arrange(plots[[1]], plots[[2]], plots[[3]],
             plots[[4]], plots[[5]], plots[[6]],
             plots[[7]], plots[[8]], plots[[9]],
             plots[[10]], plots[[11]], plots[[12]],
             widths = c(2, 1, 10))
dev.off()


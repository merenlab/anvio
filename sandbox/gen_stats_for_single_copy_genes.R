#!/usr/bin/env Rscript
#
# visualizes stuff..
#

suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(optparse))

# command line options
option_list <- list(
        make_option(c("--e_value"), default=1e-15,
                help = "e-value to retain hits [default \"%default\"]")
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

# check if the input file is accessible
if(file.access(hits) == -1){
    stop(sprintf("Specified file '%s' does not exist", hits))
}

#e_value <- 1e-15
#hits <- '~/Infant-MERGED-FINAL-5K/tarrak.db.hits'
#genes <- '~/Infant-MERGED-FINAL-5K/tarrak.db.genes'

genes_df <- data.frame(read.table(genes, header = TRUE,sep="\t"))
hits_df <- data.frame(read.table(hits, header = TRUE,sep="\t"))
hits_df <- hits_df[hits_df$e_value < e_value, ]

plots <- list()  # new empty list
i <- 1
for(source in c('Wu_et_al', 'Campbell_et_al', 'Dupont_et_al')){
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

	if(num_genomes == "0"){
    	num_genomes <- as.character(frequencies[2, ]$num)
    	percent_agrees <- frequencies[2, ]$freq * 100 / sum(frequencies$freq)
	}
	
    text = sprintf("%s:\n%.2f%% of %d genes\noccur %s times", source, percent_agrees, nrow(x), num_genomes)
    p <- ggplot() + annotate("text", x = 1, y = 1, size=8, label = text) + theme(line = element_blank(),
            text = element_blank(),
            line = element_blank(),
            title = element_blank())


	
	q <- ggplot(x, aes(x=factor(0), y=count))
	q <- q + geom_violin()
	q <- q + geom_jitter(position = position_jitter(width = .2), alpha=0.3)
	q <- q + theme_bw()
	q <- q + theme(axis.text.x = element_blank())
	q <- q + theme(axis.text.y = element_text(size = 12))
	q <- q + theme(axis.ticks.x = element_blank())
	q <- q + labs(x='', y='')
	#q <- q + scale_y_sqrt()

    r <- ggplot(x, aes(x=reorder(gene, -count), y=count, fill=count))
    r <- r + geom_bar(stat='identity', width=.9)
    r <- r + theme_bw()
    r <- r + theme(axis.text.x = element_blank())
    r <- r + theme(axis.text.y = element_text(size = 12))
    r <- r + theme(axis.ticks.x = element_blank())
    r <- r + labs(x='', y='')
    r <- r + scale_y_sqrt()

    plots[[i]] <- p
    i <- i + 1
    plots[[i]] <- q
    i <- i + 1
    plots[[i]] <- r
    i <- i + 1

}
pdf(paste(hits, '-(', e_value, ').pdf', sep=''), width=25, height=10)
grid.arrange(plots[[1]], plots[[2]], plots[[3]],
		     plots[[4]], plots[[5]], plots[[6]],
			 plots[[7]], plots[[8]], plots[[9]],
			 widths = c(2, 1, 10))
dev.off()


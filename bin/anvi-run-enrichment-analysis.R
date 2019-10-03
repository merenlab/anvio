#!/usr/bin/env Rscript

## A script to read in functional annotation data and perform the
## functional enrichment analysis

## Amy Willis, September 2019

## Usage:
# anvi-get-enriched-functions-per-pan-group --input input.txt --output output.txt

## Examples:
# Rscript anvi-get-enriched-functions-per-pan-group --input Functional_enrichment_input_file_five_groups_with_enriched_groups_column.txt --output five-group-out.txt
# Rscript anvi-get-enriched-functions-per-pan-group --input Functional_enrichment_input_file_two_groups_with_enriched_groups_column.txt --output two-group-out.txt

## TODOs:
# how do anvio users want to do the following steps:
# installing BiocManager and then running `BiocManager::install("qvalue")`

## Load required packages
required_packages <- list("tidyverse", "magrittr", "qvalue", "optparse")
installation_instructions <- list("conda install -c r r-tidyverse", "conda install -c bioconda r-magrittr", "conda install -c bioconda bioconductor-qvalue", "conda install -c conda-forge r-optparse")
names(installation_instructions) <- required_packages
missing_packages <- list()
installation_instructions_for_missing_packages <- list()
i <- 1
for (p in required_packages) {
    if(!p %in% rownames(installed.packages())){
        missing_packages[[i]] <- p
        installation_instructions_for_missing_packages[[i]] <- installation_instructions[[p]]
        i <- i+1
    }
}
if(length(missing_packages) > 0) {
    stop(cat("The following R packages are required in order to run this program, but are missing: ", paste(missing_packages, collapse = ', '), ".\n", "You can install these packages using conda by running the following commands:\n", paste(installation_instructions_for_missing_packages, collapse='\n','\n')))
}
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(magrittr))


# command line options
option_list <- list(
  make_option(c("--output"),
              help = "Output file name"),
  make_option(c("--input"),
              help = "Input file name")
)
parser <- OptionParser(option_list = option_list,
                       description = "Pangenome hypothesis testing")
arguments <- parse_args(parser,
                        positional_arguments = TRUE)
opts <- arguments$options

# check if the input file is accessible
if (file.access(opts$input) == -1) {
  stop(sprintf("Specified input file '%s' does not exist", opts$input))
}

# read in data
df_in <- readr::read_tsv(file = opts$input)

# Find number of groups
# Be smart about finding number of header rows
n_columns_before_data <- (names(df_in) %>% startsWith("p_") %>% which %>% min) - 1
G <- (ncol(df_in) - n_columns_before_data) / 2
# check if the input file makes sense
if (G %% 1 != 0) {
  stop(sprintf("The input file '%s' does not look like it came from anvi-get-enriched-functions-per-pan-group.", opt$input))
}

# Set up functions to run
run_test_no_spread <- function(df) {
  df %>%
    glm(cbind(x, N - x) ~ group, .,
        family = binomial(link = "logit")) %>%
    anova(test="Rao")
}
get_enrichment_score <- function(anova_output) {
  anova_output$Rao[2]
}
get_unadjusted_p_value <- function(anova_output) {
  ifelse(!is.na(anova_output$`Pr(>Chi)`[2]),
         anova_output$`Pr(>Chi)`[2],
         ifelse(anova_output$Rao[2] < 1e-3,
                1,
                stop("glm gave you a not-small test statistic but a missing p-value. This shouldn't happen! But it has, so please log an issue and tag @adw96 ASAP")))
}

# fit the GLMs in a nifty way using nest()
w_models <- df_in %>%
  gather(key = "type", value = "value", -c(1:n_columns_before_data)) %>%
  separate(type, into = c("type", "group"), sep = "_") %>%
  spread(type, value) %>%
  mutate(x = N * p) %>%
  select(-p) %>%
  group_by(function_accession) %>%
  nest %>%
  mutate(model = map(data, run_test_no_spread))

# get the p-values, enrichment_scores, q-values;
# then format the data and save
enrichment_output <- w_models %>%
  transmute(function_accession,
            "unadjusted_p_value" = map_dbl(model, get_unadjusted_p_value),
            "enrichment_score" = map_dbl(model, get_enrichment_score)) %>%
  mutate("adjusted_q_value" = qvalue::qvalue(unadjusted_p_value)$qvalues) %>%
  left_join(df_in, by = "function_accession") %>%
  select(names(df_in)[1],  # e.g. COG_FUNCTION
         enrichment_score,
         unadjusted_p_value,
         adjusted_q_value,
         associated_groups,
         function_accession,
         gene_clusters_ids,
         names(df_in)[n_columns_before_data + (1 : (G * 2))]) %>%
  arrange(desc(enrichment_score))

enrichment_output %>%
  write_tsv(path = opts$output)


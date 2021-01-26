# IF one day I want to use OD-SEQ to filter out outliers in 
# an MSA then the beginnings of that script are below

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("odseq")
library(odseq)
library(msa)
library(tidyverse)

setwd("/Users/mschechter/github/ACE_SO_metagenomes/analyses/test_case")
fas <- "RIBO_PHYLO_WORKFLOW/03_NR_FASTAS/Ribosomal_L16/Ribosomal_L16_all_filtered.faa"
sq <- readAAStringSet(fas)
alig <- msa(sq, method = "Muscle")

# Parameters from tutorial
y <- odseq(alig, threshold = 0.05, distance_metric = "affine", B = 1000)
good_seqs <- length(y[y == FALSE])
outliers <- length(y[y == TRUE])

# Default parameters 
y <- odseq(alig)
good_seqs <- length(y[y == FALSE])
good_seqs_names <- y[y == FALSE]

outliers <- length(y[y == TRUE])
outliers_names <- y[y == TRUE]

# Export the outliers and check them out in my own msa
transform(outliers_names) %>% # danke https://stackoverflow.com/questions/14665417/creating-a-logical-vector-from-data-frame
  as.data.frame() %>% 
  rownames_to_column(var = "header") %>% 
  as_tibble() %>%
  select(-X_data) %>%
  write_tsv("outlier_sequences.tsv")
  



alnSubset <- as(AAMultipleAlignment(unmasked(alig)[outliers_names]),
                "MsaAAMultipleAlignment")

print(alnSubset, show="complete")

print(alig, show=c("alignment"), showNames=TRUE, showConsensus=TRUE, halfNrow=200, nameWidth=20)
print(alig, show = "complete")

##Loading packages
library(corrr)
library(tidyr)
library(readr)
library(dplyr)
library(tidyverse)
library(GenomicDataCommons)
library(limma)

##Importing dataset
edn <- read_tsv("expDataNormal.tsv")
edt <- read_tsv("expDataTumor.tsv")

##Tidying up the data, adding suffixes _c for normal (control) and _t for tumour to the sample
edn <- edn %>%
  dplyr::select(!c(gene_id, gene_type))

colnames(edn) <- paste(colnames(edn),"c",sep="_")

edn <- dplyr::rename(edn, gene_name = gene_name_c)

edt <- edt %>%
  dplyr::select(!c(gene_id, gene_type))

colnames(edt)<-paste(colnames(edt),"t",sep="_")

edt <- dplyr::rename(edt, gene_name = gene_name_t)

##Merging tumour and normal tissue dataset, removing duplicates and making it matrix friendly
edmerged <- inner_join(x = edt, y = edn)

edmerged = edmerged[!duplicated(edmerged$gene_name),]

edmerged <- edmerged %>%
  column_to_rownames(var="gene_name")

##Expression EDA

#correlation
edmerged_corr <- edmerged %>% 
  cor(method = "spearman")

rplot(edmerged_corr, shape = 16) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#PCA
pca_matrix <- edmerged %>% 
  as.matrix() %>% 
  t()

sample_pca <- prcomp(pca_matrix)

pc_scores <- sample_pca$x

pc_scores <- pc_scores %>% 
  as_tibble(rownames = "tissue")

pc_scores %>% 
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point()
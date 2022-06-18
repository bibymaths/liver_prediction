##################################################################################################################################################################
##  This script downloads, and preprocess the TCGA-LIHC project data for transcriptome profiling using TCGA Biolinks for Gene Expression Analysis using DESeq2  ##
################################################################################################################################################################## 
 
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("TCGAbiolinks") 
 
library(TCGAbiolinks) 
 
##################################################################
##             Querying the data in GDC Data Portal             ##
################################################################## 

mirna_query <- GDCquery(project = "TCGA-LIHC",
                        data.category = "Transcriptome Profiling",
                        data.type = "miRNA Expression Quantification",
                        #workflow.type = "BCGSC miRNA Profiling",
                        experimental.strategy = "miRNA-Seq")  
 
#################################################################
##               Downloading the expression data               ##
#################################################################

GDCdownload(mirna_query, method = "api", files.per.chunk = 100,
            directory = "~/Desktop/TCGA/miRNA/") 

miR_df <- GDCprepare(mirna_query, directory = "~/Desktop/TCGA/miRNA/")
  
#################################################################
##      Removing unnecessary coloumns, and keeping counts      ##
#################################################################

rownames(miR_df) <- miR_df$miRNA_ID
miR_df <- miR_df[,-1]
number_cols <- ncol(miR_df)
subset <- seq(from = 1, to = number_cols, by = 3)
miR_df <- miR_df[, subset] 

##################################################################
##                           Case IDs                           ##
################################################################## 

colnames(miR_df) <- gsub(".*_","",colnames(miR_df)) 

##################################################################
##               Matching the coloums to metadata               ##
################################################################## 

miR_meta <- mirna_query[[1]][[1]]
miR_meta <- miR_meta[,c("cases", "tissue.definition")]
rownames(miR_meta) <- colnames(miR_df)
table(rownames(miR_meta) == miR_meta$cases) 


miR_meta$tissue.definition <- as.character(miR_meta$tissue.definition)
table(miR_meta$tissue.definition) 

##################################################################
##               Removing and renaming conditions               ##
################################################################## 

metastatic_key <- miR_meta[which(miR_meta$tissue.definition == "Metastatic"),]
miR_meta <- miR_meta[!miR_meta$tissue.definition == metastatic_key$tissue.definition,]
miR_df <- miR_df[, -grep(paste0(metastatic_key$cases), colnames(miR_df))]

miR_meta$tissue.definition <- gsub("Primary solid Tumor", "Tumor", miR_meta$tissue.definition)
miR_meta$tissue.definition <- gsub("Solid Tissue Normal", "Normal", miR_meta$tissue.definition)
miR_meta$tissue.definition <- as.factor(miR_meta$tissue.definition)
levels(miR_meta$tissue.definition)
colnames(miR_meta) <- c("cases", "Condition") 

##################################################################
##                          Tidying up                          ##
################################################################## 

rm(mirna_query)
rm(subset)
rm(number_cols)
rm(metastatic_key)


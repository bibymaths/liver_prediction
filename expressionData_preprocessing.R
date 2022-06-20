##################################################################################################################################################################
##  This script downloads, and preprocess the TCGA-LIHC project data for transcriptome profiling using TCGA Biolinks for Gene Expression Analysis using DESeq2  ##
################################################################################################################################################################## 
  
#################################################################
##                 INSTALLING/LOADING PACKAGES                 ##
################################################################# 
##-----------------------------------------------------------------------------
##  For BioConductor packages in R, install them using BiocManager::install   -
##-----------------------------------------------------------------------------
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager") 

if (!require("pacman")) install.packages("pacman") 
pacman::p_load(DESeq2, RColorBrewer, pheatmap,  
               ggplot2, biomaRt, PCAtools, DT,  
               IHW, apeglm, EnhancedVolcano,  
               ComplexHeatmap, TCGAbiolinks,  
               SummarizedExperiment, pheatmap, 
               reshape2)
 
##################################################################
##             Querying the data in GDC Data Portal             ##
################################################################## 
mrna_query <- GDCquery(project = "TCGA-LIHC",
                       data.category = "Transcriptome Profiling",
                       data.type = "Gene Expression Quantification",
                       workflow.type = "STAR - Counts",
                       experimental.strategy = "RNA-Seq")

 
#################################################################
##               Downloading the expression data               ##
#################################################################  
GDCdownload(mrna_query, method = "api", files.per.chunk = 100,
            directory = "~/Desktop/TCGA/miRNA/") 

dat <- GDCprepare(query = mrna_query, save = TRUE,  
                  save.filename = "exp.rda",  
                  directory = "~/Desktop/TCGA/miRNA/")

# exp matrix
rna <- as.data.frame(SummarizedExperiment::assay(dat))
# clinical data
clinical <- data.frame(dat@colData)

# Check how many tumor and control sample are there:
## data are stored under "definition" column in the clinical dataset.
table(clinical$definition)

# Also from sample id it is possible to count normal and tumor samples:
table(substr(colnames(rna),14,14))

#The count matrix (rna dataset) and the rows of the column data (clinical dataset) MUST be in the same order.

# to see whether all rows of clinical are present in rna datset
all(rownames(clinical) %in% colnames(rna))

# whether they are in the same order:
all(rownames(clinical) == colnames(rna))

# if not reorder them by:
#rna <- rna[, rownames(clinical)]
#all(rownames(clinical) == colnames(rna))

  
#_______Making_Expression_Object__________#

#We will use the column “definition”, as the grouping variable for gene expression analysis. 

# replace spaces with "_" in levels of definition column
clinical$definition <-  gsub(" ", "_", clinical$definition)

# making the definition column as factor
clinical$definition <- as.factor(clinical$definition)
# relevling factor to ensure tumors would be compared to normal tissue.
levels(clinical$definition)
#
clinical$definition <- relevel(clinical$definition, ref = "Solid_Tissue_Normal")

# Making DESeqDataSet object which stores all experiment data
dds <- DESeqDataSetFromMatrix(countData = rna,
                              colData = clinical,
                              design = ~ definition)
#dds

# prefilteration: it is not necessary but recommended to filter out low expressed genes

keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# data tranfromation
vsd <- vst(dds, blind=FALSE)

# making PC object
p <- pca(assay(vsd), metadata = colData(vsd), removeVar = 0.1)

# create PCA plot for PCA1 and PCA2
biplot(p, colby = "definition", lab = NULL, legendPosition = 'right')


# For all top 10 possible combination 
pairsplot(p,
          components = getComponents(p, c(1:10)),
          triangle = TRUE, trianglelabSize = 12,
          hline = 0, vline = 0,
          pointSize = 0.4,
          gridlines.major = FALSE, gridlines.minor = FALSE,
          colby = 'definition',
          title = 'Pairs plot', plotaxes = FALSE,
          margingaps = unit(c(-0.01, -0.01, -0.01, -0.01), 'cm'))

#' ---
#' title: "Differential Expression Gene Analysis for mRNA data from TCGA-LIHC project"
#' author: "Abhinav Mishra"
#' output:
#'   pdf_document:
#'     latex_engine: lualatex
#' ---

#' This script downloads, and preprocess, and analyzes the TCGA-LIHC project data for transcriptome profiling for Gene Expression Analysis

#################################################################
##                 INSTALLING/LOADING PACKAGES                 ##
################################################################# 
##-----------------------------------------------------------------------------
##  For BioConductor packages in R, install them using BiocManager::install   -
##-----------------------------------------------------------------------------

if (!require("pacman", quietly = TRUE))  
  install.packages("pacman")  
if (!require("BiocManager", quietly = TRUE))  
  install.packages("BiocManager")  
## ETA: 30 minutes
BiocManager::install(c("DESeq2", "biomaRt", "PCAtools",
"IHW", "apeglm", "EnhancedVolcano",
"ComplexHeatmap", "TCGAbiolinks",
"SummarizedExperiment"))

pacman::p_load(DESeq2, RColorBrewer, pheatmap,  
               ggplot2, biomaRt, PCAtools, DT,  
               IHW, apeglm, EnhancedVolcano,  
               ComplexHeatmap, TCGAbiolinks,  
               SummarizedExperiment,reshape2)
 
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
## ETA: 10-15 minutes
GDCdownload(mrna_query, method = "api", files.per.chunk = 100,
            directory = "/Users/abhinavmishra/Documents/Freie/DSLS/Project/TCGA/mRNA") 
## ETA : 2-3 minutes 
# Starting to add information to samples
# => Add clinical information to samples
# => Adding TCGA molecular information from marker papers
# => Information will have prefix 'paper_' 
# Available assays in SummarizedExperiment : 
#   => unstranded
# => stranded_first
# => stranded_second
# => tpm_unstrand
# => fpkm_unstrand
# => fpkm_uq_unstrand
# => Saving file: exp.rda
# => File saved
dat <- GDCprepare(query = mrna_query, save = TRUE,  
                  save.filename = "exp.rda",  
                  directory = "/Users/abhinavmishra/Documents/Freie/DSLS/Project/TCGA/mRNA")

# exp matrix considering only unstranded list
rna <- as.data.frame(SummarizedExperiment::assays(dat)$unstranded)
# clinical data  

df <- data.frame(dat@colData) 
 
#removing the irrelevant definition type = "recurrence" 
clinical <- subset(df, df$definition != "Recurrent Solid Tumor") 

#calculating NA values
colSums(is.na(clinical))

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

#if not reorder them by:
rna <- rna[, rownames(clinical)]
all(rownames(clinical) == colnames(rna))

##################################################################
##    Data transformation, quality control and normalization    ##
################################################################## 
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
  

    
#################################################################
##                    Hiearchial Clustering                    ##
################################################################# 
  
normal_idx <- substr(colnames(assay(vsd)),14,14) == "1"
n_sample <- assay(vsd)[, c(normal_idx) ]
colnames(n_sample) <- paste("NT_", substr(colnames(n_sample),1,12))

# Dissimilarity matrix calculation
sampleDists <- dist(t(n_sample))
sampleDistMatrix <- as.matrix(sampleDists)
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
# heatmap visualization
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

#################################################################
##            Differential Gene Expression Analysis            ##
################################################################# 
 
#_______DE_analysis________#

dds <- DESeq(dds)   
## estimating size factors
## estimating dispersions
## gene-wise dispersion estimates
## mean-dispersion relationship
## final dispersion estimates
## fitting model and testing
## would take some time ~ 15 minutes max

#default method 
#alpha controls FDR
res <- results(dds, lfcThreshold = 1.5, altHypothesis = "greaterAbs", alpha = 0.05 )   
resLA <- results(dds, lfcThreshold=1.5, altHypothesis="lessAbs", alpha = 0.05)
resG <- results(dds, lfcThreshold=1.5, altHypothesis="greater", alpha = 0.05)
resL <- results(dds, lfcThreshold=1.5, altHypothesis="less", alpha = 0.05)
## making easier for gene id conversion
res@rownames <- gsub("\\..*","",res@rownames)

## Let's look at the results table         ##
head(res)  
## Summary of differential gene expression ##
summary(res)   
 
#################################################################
##           Plot Gene Count: top 6 genes by p-value           ##
#################################################################
## compare the normalized counts ##
## between treated and control   ## 
## for our top 6 genes           ## 
 
par(mfrow=c(2,3))

plotCounts(dds, gene="ENSG00000043355.12", intgroup="definition")
plotCounts(dds, gene="ENSG00000187730.9", intgroup="definition")
plotCounts(dds, gene="ENSG00000253898.1", intgroup="definition")
plotCounts(dds, gene="ENSG00000147257.15", intgroup="definition")
plotCounts(dds, gene="ENSG00000164362.21", intgroup="definition")
plotCounts(dds, gene="ENSG00000106031.9", intgroup="definition")

#result Log fold change shrinkage method  
#suitable for logfc based visualization
resLFC <- lfcShrink(dds,  
                    coef=resultsNames(dds)[2],  
                    type="apeglm")  

resLFC.Ordered<-resLFC[with(resLFC,  
                            order(abs(log2FoldChange), 
                                  padj, decreasing = TRUE)), ]

#################################################################
##                        ID Conversion                        ##
#################################################################
res$ensembl <- sapply( strsplit( rownames(res), split="\\+" ), "[", 1 )
ensembl = useMart( "ensembl", dataset = "hsapiens_gene_ensembl" )
genemap <- getBM( attributes = c('ensembl_gene_id','entrezgene_id', 
                                 'hgnc_id','hgnc_symbol','external_gene_name'),
                  filters = 'ensembl_gene_id',
                  values = res$ensembl,
                  mart = ensembl )
idx <- match( res$ensembl, genemap$ensembl_gene_id )
res$entrez <- genemap$entrezgene[ idx ]
res$hgnc_symbol <- genemap$hgnc_symbol[ idx ] 
res$hgnc_id <- genemap$hgnc_id[ idx ]
           
#################################################################
##                   Up/Down Regulated Genes                   ## 
#################################################################  
###  subset the results and then sort it by the log2 fold change  
###  estimate to get the significant genes    
## Sort summary list by p-value            ##
resOrdered <- res[order(res$padj),]
resSig <- subset(resOrdered, padj < 0.1)
###  strongest upregulation       
tail(resSig[ order( resSig$log2FoldChange ), ])  
###  strongest downregulation
head( resSig[ order( resSig$log2FoldChange ), ])
 
##################################################################
##                       Diagnostic Plots                       ##
################################################################## 
# overview for two-group comparison  
plotMA( res, ylim = c(-1, 1) )  
# disperson estimates
plotDispEsts( dds, ylim = c(1e-6, 1e1) )  
# histogram of p-values 
hist( res$pvalue, breaks=20, col="grey" )
# create bins using the quantile function
qs <- c( 0, quantile( res$baseMean[res$baseMean > 0], 0:7/7 ) )
# "cut" the genes into the bins
bins <- cut( res$baseMean, qs )
# rename the levels of the bins using the middle point
levels(bins) <- paste0("~",round(.5*qs[-1] + .5*qs[-length(qs)]))
# calculate the ratio of p values less than .01 for each bin
ratios <- tapply( res$pvalue, bins, function(p) mean( p < .01, na.rm=TRUE ) )
# plot these ratios
barplot(ratios, xlab="mean normalized count", ylab="ratio of small p-values") 
#ploting genes differentially expressed  
 
## MA-plot visualize gene significantly                   ##
## (blue dots with p < 0.05) up-regulated (log2FC > 1.5)  ##
## down-regulated (log2FC < -1.5)                         ##
ylim <- c(-6.5,6.5)
drawLines <- function() abline(h=c(-1.5,1.5),col="dodgerblue",lwd=2)
plotMA(resLFC, ylim=ylim); drawLines() 
 
## Alternative Shrinkage Estimators ## 

resultsNames(dds) 
resNorm <- lfcShrink(dds, coef=2, type="normal")
resAsh <- lfcShrink(dds, coef=2, type="ashr") 
par(mfrow=c(1,3), mar=c(4,4,2,1))
xlim <- c(1,1e5); ylim <- c(-3,3)
plotMA(resLFC, xlim=xlim, ylim=ylim, main="apeglm")
plotMA(resNorm, xlim=xlim, ylim=ylim, main="normal")
plotMA(resAsh, xlim=xlim, ylim=ylim, main="ashr") 
 
## Tests of log2 fold change above or below a threshold       ## 
## greaterAbs - |beta|>x - tests are two-tailed               ##
## lessAbs - |beta|<x - p values are the maximum              ##
## of the upper and lower tests                               ##  
## greater - |beta| >x                                        ##
## less - |beta| <−x                                          ##
par(mfrow=c(2,2),mar=c(2,2,1,1))
ylim <- c(-2.5,2.5) 
drawLines <- function() abline(h=c(-.5,.5),col="dodgerblue",lwd=2)
plotMA(res, ylim=ylim); drawLines()
plotMA(resLA, ylim=ylim); drawLines()
plotMA(resG, ylim=ylim); drawLines()
plotMA(resL, ylim=ylim); drawLines() 

## COUNT OUTLIER DETECTION ## 
W <- res$stat
maxCooks <- apply(assays(dds)[["cooks"]],1,max)
idx <- !is.na(W)
plot(rank(W[idx]), maxCooks[idx], xlab="rank of Wald statistic", 
     ylab="maximum Cook's distance per gene",
     ylim=c(0,5), cex=.4, col=rgb(0,0,0,.3))
m <- ncol(dds)
p <- 3
abline(h=qf(.99, p, m - p))  
##----------------------------------------------------------------
##                      Filtering Criteria                       -
##----------------------------------------------------------------
## plot the -log10 p values from all genes  ##
## over the normalized mean counts          ##
plot(res$baseMean+1, -log10(res$pvalue),
     log="x", xlab="mean of normalized counts",
     ylab=expression(-log[10](pvalue)),
     ylim=c(0,30),cex=.4, col=rgb(0,0,0,.3))
 
###########################################################################
###########################################################################
###                                                                     ###
###                            VOLCANO PLOTS                            ###
###                                                                     ###
###########################################################################
###########################################################################
par(mfrow=c(2,2)) 

with(res, plot(log2FoldChange, -log10(pvalue),  
               pch=20, main="Volcano plot (p < 0.01 & logFC > 2 )",  
               xlim=c(-5,5)))
with(subset(res, padj<.01 ),  
     points(log2FoldChange,-log10(pvalue), pch=20, col="blue"))
with(subset(res, padj<.01 & abs(log2FoldChange)>2),  
     points(log2FoldChange, -log10(pvalue), pch=20, col="red"))  

with(res, plot(log2FoldChange, -log10(pvalue),  
               pch=20, main="Volcano plot (p < 0.05 & logFC > 2 )",  
               xlim=c(-5,5)))
with(subset(res, padj<.05 ),  
     points(log2FoldChange,-log10(pvalue), pch=20, col="blue"))
with(subset(res, padj<.05 & abs(log2FoldChange)>2),  
     points(log2FoldChange, -log10(pvalue), pch=20, col="red"))  

with(res, plot(log2FoldChange, -log10(pvalue),  
               pch=20, main="Volcano plot (p < 0.01 & logFC > 3 )",  
               xlim=c(-5,5)))
with(subset(res, padj<.01 ),  
     points(log2FoldChange,-log10(pvalue), pch=20, col="blue"))
with(subset(res, padj<.01 & abs(log2FoldChange)>3),  
     points(log2FoldChange, -log10(pvalue), pch=20, col="red"))  

with(res, plot(log2FoldChange, -log10(pvalue),  
               pch=20, main="Volcano plot (p < 0.05 & logFC > 3 )",  
               xlim=c(-5,5)))
with(subset(res, padj<.05 ),  
     points(log2FoldChange,-log10(pvalue), pch=20, col="blue"))
with(subset(res, padj<.05 & abs(log2FoldChange)>3),  
     points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
##################################################################
##                        Saving Results                        ##
##################################################################
## Complete Results 
write.csv( as.data.frame(res), file="results.csv" )  
  
## Functions to get the result of x genes and save it in a .csv format
down_gene <- function(x){ 
  down <- head( resSig[ order( resSig$log2FoldChange ), ],x) 
  write.csv( as.data.frame(down), file = sprintf("results_down_%d.csv", x)) 
  sprintf("The file with %d downregulated genes has been created.", x)   
} 
up_gene <- function(x){ 
  up <- tail( resSig[ order( resSig$log2FoldChange ), ],x) 
  write.csv( as.data.frame(up), file = sprintf("results_up_%d.csv", x))
  sprintf("The file with %d upregulated genes has been created.", x)   
}  
 
## top 100, 200, 300 genes 
## change for another number of genes you want
x <- seq.int(100, 300, by = 100)  
## a loop for creating all files at once 
for(i in x ){  
  down_gene(i)
  up_gene(i) 
} 
 


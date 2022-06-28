library(dplyr)
library(BiocManager)
library(TCGAbiolinks)
library(SummarizedExperiment)
library(limma)
library(readr)
library(wateRmelon)
library(ExperimentHub)
library(DMRcate)
library(RColorBrewer)
library(minfi)
library(DMRcate)

setwd("C:/Users/Krissi/Desktop/TCGA_LIHC/")


###################################download methylation data (beta values) from TCGA##############################################

mdna_query <- GDCquery(project = "TCGA-LIHC",
                       data.category = "DNA Methylation",
                       data.type = "Methylation Beta Value",
                       platform = "Illumina Human Methylation 450")


GDCdownload(mdna_query, method = "api", files.per.chunk = 50,
            directory = "biolinksData/")

data <- GDCprepare(query = mdna_query, save = TRUE,  save.filename = "dnaM_client.rda",  directory = "clientData/", summarizedExperiment = TRUE)

##load data
load("dnaM.rda")


met <- as.data.frame(SummarizedExperiment::assay(data))
clinical = colDataPrepare(data@colData$samples)
clinical = clinical[clinical$patient %in%  clinical[clinical$sample_type=="Solid Tissue Normal", "patient"],]
###################################download info data from TCGA##############################################

barcodes = data@colData$samples
barcodes = unlist(lapply(barcodes, substr, start = 1, stop = 12))

clinical_query <- GDCquery(project = "TCGA-LIHC",
                           data.category = "Clinical", 
                           barcode = barcodes)

GDCdownload(clinical_query, method = "api", files.per.chunk = 100,
            directory = "ClinicalData/")

info <- GDCprepare_clinic(query = clinical_query, clinical.info = "patient", directory = "ClinicalData/")
info.stage_event <- GDCprepare_clinic(query = clinical_query, clinical.info = "stage_event", directory = "ClinicalData/")

info = merge(info, info.stage_event, by = 'bcr_patient_barcode')

########################################preprocess methylation data############################################

# get  450k annotation data for hg38
hg38_anno = read_tsv("HM450.hg38.manifest.tsv.gz")

## remove probes with NA
probe.na <- rowSums(is.na(met))

table(probe.na == 0)
#FALSE   TRUE 
#191453  294124  
#choose those without NA values in rows
probe <- probe.na[probe.na == 0]
met <- met[row.names(met) %in% names(probe), ]

##remove probes that match to chromosome  X and Y 
keep <- !(row.names(met) %in% hg38_anno$probeID[hg38_anno$CpG_chrm %in% c("chrX","chrY")])
table(keep)
#FALSE   TRUE 
#5568   288556  
met <- met[keep, ]

## remove SNPs overlapped probe with MAF > 0.01
table(hg38_anno$MASK_snp5_GMAF1p)
#FALSE   TRUE 
#462832  22745  
# probes without snp with MAF > 0.01
no.snp.probe <- hg38_anno$probeID[!(hg38_anno$MASK_snp5_GMAF1p)]

# filter met
met <- met[row.names(met) %in% c(no.snp.probe), ]

#remove no-further needed dataset
rm(probe.na)
rm(probe)
rm(keep)
rm(no.snp.probe)


#intra-sample normalisation procedure, correcting the bias of type-2 probe values
#BMIQ(met, design.v)
##get design.v?

##############################################################EDA##########################################################################

densityPlot(as.matrix(met), sampGroups = clinical$sample_type)

pal <- brewer.pal(8,"Dark2")
limma::plotMDS(met, top=1000, gene.selection="common", 
        col=pal[factor(clinical$sample_type)], labels = NULL, pch = 19)

#####################################differential methylation analysis (loci and region)###################################################
unique(clinical$sample_type)
table(clinical$sample_type)
##Primary Tumor Solid Tissue Normal 
##50                  50

##check ordering
met = met[,row.names(clinical)]

design <- model.matrix(~ clinical$sample_type)
fit <- lmFit(met, design)
fit <- eBayes(fit)
top_10 = topTable(fit)

# plot the top 4 most significantly differentially methylated CpGs 
par(mfrow=c(2,2))
sapply(rownames(top_10)[1:4], function(cpg){
  plotCpg(met, cpg=cpg, pheno=clinical$sample_type, ylab = "Beta values")
})


myannotation <- cpg.annotate("array", as.matrix(met), what = "Beta", arraytype = "450K",
                             analysis.type="differential", design=design, coef=2, fdr = 0.01)

dmrcoutput <- dmrcate(myannotation, lambda=2000, C=2)

results.ranges <- extractRanges(dmrcoutput, genome = "hg38")
DMR.plot(ranges=results.ranges, dmr=1, CpGs=as.matrix(met), what="Beta",
         arraytype = "450", phen.col=clinical$sample_type, genome="hg38")

results.ranges

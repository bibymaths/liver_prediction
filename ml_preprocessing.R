library(SummarizedExperiment)
library(TCGAbiolinks)
library(curatedTCGAData)
library(dplyr)
library(tidyr)
library(readr)

setwd("C:/Users/Krissi/Desktop/TCGA_LIHC/LIHC-project")


###################################################load and filter methylation data#################################################
load("data/dnaM.rda")
meth_data = data 

met <- as.data.frame(SummarizedExperiment::assay(meth_data))
clinical_meth = colDataPrepare(data@colData$samples)
clinical_meth = clinical_meth[clinical_meth$sample_type!="Recurrent Tumor",]

CpGs_sig = read_csv("data/top500_cpgs.csv")

met_sig = met[CpGs_sig$top500_cpgs,]
met_sig = met_sig[,row.names(clinical_meth)]
colnames(met_sig) = unlist(lapply(colnames(met_sig), substr, start = 1, stop = 15))

##################################################load and filter expression data###################################################
load("data/exp.rda")
exp_data = data 

rna <- as.data.frame(SummarizedExperiment::assays(exp_data)$unstranded)

# clinical data  
clinical_exp <- data.frame(exp_data@colData) 

up = read_csv("data/results_up_200.csv")
down = read_csv("data/results_down_200.csv")
sign_genes = rbind(up,down)

rna_sig = rna[sign_genes$...1,]
colnames(rna_sig) = unlist(lapply(colnames(rna_sig), substr, start = 1, stop = 15))

##################################################load and filter CNV data###################################################


lihc_gistic <- curatedTCGAData::curatedTCGAData("LIHC",
                                                assays=("GISTIC_Peaks"),
                                                version = "1.1.38",
                                                dry.run=FALSE)

gistic <- as.data.frame(assay(lihc_gistic))

gistic <- tibble::rownames_to_column(gistic, "position")

gistic <- gistic %>%
  separate(col = position, into = c("seqnames", "Position"), sep = ":") %>%
  separate(col = Position, into = c("Start", "End"), sep = "-")

colnames(gistic) = unlist(lapply(colnames(gistic), substr, start = 1, stop = 15))
rownames(gistic) = paste0(gistic$seqnames,":",gistic$Start, "-", gistic$End)

##################################################merge datasets###################################################

##get intersecting cases for  expression and methylation data (for healthy/cancer classification)
intersecting_samples_exp_meth = intersect(colnames(met_sig), colnames(rna_sig))

##get intersecting cases for  expression, methylation, CNV data (for alive/dead classification and survival analysis)
intersecting_samples = intersect(intersect(colnames(met_sig), colnames(rna_sig)),colnames(gistic))

##filter cases
met_sig = met_sig[,intersecting_samples_exp_meth]
rna_sig = rna_sig[,intersecting_samples_exp_meth]

##merge datasets
features_exp_meth = rbind(rna_sig,met_sig) 
features_exp_meth = t(features_exp_meth)

##filter cases
rna_sig = rna_sig[,intersecting_samples]
gistic =gistic[,intersecting_samples]
met_sig = met_sig[,intersecting_samples]

##merge datasets
features_complete = rbind(rna_sig,met_sig,gistic) 
features_complete = t(features_complete)


##add label columns
features_exp_meth = cbind(features_exp_meth, clinical_meth[match(row.names(features_exp_meth), clinical_meth$sample.aux),]$sample_type)
colnames(features_exp_meth)[901] = "sample_type"

features_exp_meth = cbind(features_exp_meth,clinical_meth[match(row.names(features_exp_meth), clinical_meth$sample.aux),]$vital_status)
colnames(features_exp_meth)[902] = "vital_status"

features_complete = cbind(features_complete,clinical_meth[match(row.names(features_complete), clinical_meth$sample.aux),]$vital_status)
colnames(features_complete)[960] = "vital_status"  

features_complete = cbind(features_complete,clinical_meth[match(row.names(features_complete), clinical_meth$sample.aux),]$days_to_death)
colnames(features_complete)[961] = "days_to_death" 


write.table(as.data.frame(features_exp_meth), "data/features_exp_meth.tsv", sep = "\t", quote = F)
write.table(as.data.frame(features_complete), "data/features_complete.tsv", sep = "\t", quote = F)

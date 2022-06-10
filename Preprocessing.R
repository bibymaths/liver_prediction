setwd("D:/TCGA_LIHC/")

library(tidyr)
library(readr)
library(dplyr)
library(GenomicDataCommons)
        
##load metainformation for expression and methylation Info
expInfo = read_tsv("expressionData_gdc_sample_sheet.tsv")
methInfo = read_tsv("methylationData_gdc_sample_sheet.tsv")
clinInfo = read_tsv("clinical.tsv")

##divide into cancer and control Infoframes
expInfoNormal = expInfo %>% filter(expInfo$`Sample Type`=='Solid Tissue Normal')
expInfoTumor = expInfo %>% filter(expInfo$`Sample Type`=='Primary Tumor')
methInfoNormal = methInfo %>% filter(methInfo$`Sample Type`=='Solid Tissue Normal')
methInfoTumor = methInfo %>% filter(methInfo$`Sample Type`=='Primary Tumor')

##get intersecting case IDs
intersectingCases = intersect(expInfoNormal$`Case ID`, methInfoNormal$`Case ID`)

##subset Infosets by case IDs
expInfoNormal = expInfoNormal %>% filter(expInfoNormal$`Case ID` %in% intersectingCases)
expInfoTumor = expInfoTumor %>% filter(expInfoTumor$`Case ID` %in% intersectingCases)
methInfoNormal = methInfoNormal %>% filter(methInfoNormal$`Case ID` %in% intersectingCases)
methInfoTumor = methInfoTumor %>% filter(methInfoTumor$`Case ID` %in% intersectingCases)

##cache information
gdc_cache()
##downlaod data from GDC into cache
gdcdata(expInfoNormal$'File ID')
gdcdata(expInfoTumor$'File ID')
gdcdata(methInfoNormal$'File ID')
gdcdata(methInfoTumor$'File ID')

##load methylation data into dataframes
methDataTumor = read_tsv(paste0(gdc_cache(), "/", methInfoTumor$`File ID`[1], "/", methInfoTumor$`File Name`[1]), comment = "#", col_names = c("CpG", methInfoTumor$`Case ID`[1]))
methDataNormal = read_tsv(paste0(gdc_cache(), "/", methInfoNormal$`File ID`[1], "/", methInfoNormal$`File Name`[1]), comment = "#", col_names = c("CpG", methInfoNormal$`Case ID`[1]))
for(row in 2:nrow(methInfoTumor)){
  tmp_df_tumor = read_tsv(paste0(gdc_cache(), "/", methInfoTumor$`File ID`[row], "/", methInfoTumor$`File Name`[row]), comment = "#", col_names = c("CpG", methInfoTumor$`Case ID`[row]))
  tmp_df_normal = read_tsv(paste0(gdc_cache(), "/", methInfoNormal$`File ID`[row], "/", methInfoNormal$`File Name`[row]), comment = "#", col_names = c("CpG", methInfoNormal$`Case ID`[row]))
  methDataTumor = merge(methDataTumor, tmp_df_tumor, by = 'CpG', all = TRUE)
  methDataNormal = merge(methDataNormal, tmp_df_normal, by = 'CpG', all = TRUE)
}

##write zipped methylation files
write_tsv(methDataNormal, "./data/methDataNormal.tsv.gz")
write_tsv(methDataNormal, "./data/methDataNormal.tsv.gz")

##load expression data into dataframes
expDataTumor = read_tsv(paste0(gdc_cache(), "/", expInfoTumor$`File ID`[1], "/", expInfoTumor$`File Name`[1]), comment = "#", col_select = c("gene_id", "gene_name", "gene_type", "fpkm_unstranded"))[-(1:4),]
names(expDataTumor)[4] = expInfoTumor$`Case ID`[1]
expDataNormal = read_tsv(paste0(gdc_cache(), "/", expInfoNormal$`File ID`[1], "/", expInfoNormal$`File Name`[1]), comment = "#", col_select = c("gene_id", "gene_name", "gene_type", "fpkm_unstranded"))[-(1:4),]
names(expDataNormal)[4] = expInfoNormal$`Case ID`[1]
for(row in 2:nrow(methInfoTumor)){
  tmp_df_tumor = read_tsv(paste0(gdc_cache(), "/", expInfoTumor$`File ID`[row], "/", expInfoTumor$`File Name`[row]), comment = "#", col_select = c("gene_id", "fpkm_unstranded"))[-(1:4),]
  names(tmp_df_tumor)[2] = expInfoTumor$`Case ID`[row]
  tmp_df_normal = read_tsv(paste0(gdc_cache(), "/", expInfoNormal$`File ID`[row], "/", expInfoNormal$`File Name`[row]), comment = "#", col_select = c("gene_id", "fpkm_unstranded"))[-(1:4),]
  names(tmp_df_normal)[2] = expInfoNormal$`Case ID`[row]
  expDataTumor = merge(expDataTumor, tmp_df_tumor, by = 'gene_id', all = TRUE)
  expDataNormal = merge(expDataNormal, tmp_df_normal, by = 'gene_id', all = TRUE)
}

##write zipped expression files
write_tsv(expDataNormal, "./data/expDataNormal.tsv.gz")
write_tsv(expDataTumor, "./data/expDataTumor.tsv.gz")


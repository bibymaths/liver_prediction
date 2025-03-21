library(dplyr)
library(BiocManager)
library(TCGAbiolinks)
library(SummarizedExperiment)
library(limma)
library(readr)
library(ExperimentHub)
library(DMRcate)
library(RColorBrewer)
library(minfi)
library(stringr)
library(ggplot2)
library(EnhancedVolcano)
library(biomaRt)
library(liftOver)
library(coMET)

setwd("C:/Users/Krissi/Desktop/TCGA_LIHC/LIHC-project")


###################################download methylation data (beta values) from TCGA##############################################

mdna_query <- GDCquery(project = "TCGA-LIHC",
                       data.category = "DNA Methylation",
                       data.type = "Methylation Beta Value",
                       platform = "Illumina Human Methylation 450")


GDCdownload(mdna_query, method = "api", files.per.chunk = 50,
            directory = "biolinksData/")

data <- GDCprepare(query = mdna_query, save = TRUE,  save.filename = "dnaM_client.rda",  directory = "clientData/", summarizedExperiment = TRUE)

##load data
load("data/dnaM.rda")

met <- as.data.frame(SummarizedExperiment::assay(data))
clinical = colDataPrepare(data@colData$samples)
clinical = clinical[clinical$patient %in%  clinical[clinical$sample_type=="Solid Tissue Normal", "patient"],]

remove(data)
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
hg38_anno = read_tsv("data/HM450.hg38.manifest.tsv.gz")

## remove probes with NA
probe.na <- rowSums(is.na(met))

table(probe.na == 0)
#FALSE   TRUE 
#191453  294124  
#choose those without NA values in rows
probe <- probe.na[probe.na == 0]
met <- met[row.names(met) %in% names(probe), ]

##remove probes that match to chromosome  X and Y to remove sex bias
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


table(clinical$sample_type)
##Primary Tumor Solid Tissue Normal 
##50                  50

##check ordering
met = met[,row.names(clinical)]

mval <- t(apply(met, 1, function(x) log2(x/(1-x))))


##############################################################EDA##########################################################################

pdf("plots/densityBeta.pdf")
densityPlot(as.matrix(met), sampGroups = clinical$sample_type)
dev.off()

pdf("plots/MDS.pdf")
pal <- brewer.pal(8,"Dark2")
limma::plotMDS(met, top=1000, gene.selection="common", 
        col=pal[factor(clinical$sample_type)], labels = NULL, pch = 19)
dev.off()

#####################################differential methylation analysis (loci and region)###################################################

patient =  factor(make.names(clinical$patient))
sample_type = factor(clinical$sample_type)

design <- model.matrix(~sample_type)

##fit the model
fit <- lmFit(mval, design)
fit <- eBayes(fit)

##get summary stats
infinite = topTable(fit, number = 288556)

pdf("plots/MA_all.pdf")
ggplot(infinite) + geom_point(aes(x=AveExpr, y=logFC))
dev.off()

pdf("plots/enhancedVolcano.pdf")
EnhancedVolcano(infinite,
                lab = rownames(infinite),
                x = 'logFC',
                y = 'adj.P.Val')
dev.off()


##filter and order by p_value and beta value difference between groups, delta_b > 0.2
diff = data.frame(MeanBeta_diff = rowMeans(met[,sample_type=="Solid Tissue Normal"]) - rowMeans(met[,sample_type=="Primary Tumor"]))
merged_stats = merge(infinite, diff, by = 0)
row.names(merged_stats) = merged_stats[,1]
merged_stats = merged_stats[,-1]
merged_stats = merged_stats %>% filter(P.Value < 0.005 & MeanBeta_diff > 0.2) %>% arrange(P.Value)
head(merged_stats)
cat("Number of differentially methylated loci: ", nrow(merged_stats))

# plot the top 4 most significantly differentially methylated CpGs 
pdf("plots/sigCpG.pdf")
par(mfrow=c(2,2))
sapply(rownames(merged_stats)[1:4], function(cpg){
  plotCpg(met, cpg=cpg, pheno=sample_type, ylab = "Beta values")
})
dev.off()

###annotates to hg19!!!!
myannotation <- cpg.annotate("array", mval, what = "M", arraytype = "450K",
                             analysis.type="differential", design=design, coef=ncol(design), fdr = 0.01)


##get correct coordinates from hg38_anno
true_coord = hg38_anno[,c(1,2,3,5)]
anno_df = data.frame(myannotation@ranges)
anno_df$probeID = names(myannotation@ranges)
merged_df = merge(true_coord, anno_df)
merged_df = merged_df[,c(1,2,3,9,10,11,12,13)]
merged_df$start = merged_df$CpG_beg + 1
merged_df$end = merged_df$CpG_beg + 1
Cpg_IDs = merged_df$probeID
merged_df = merged_df[,-c(1,3)]
colnames(merged_df) = c("seqnames","strand","stat","diff","ind.fdr","is.sig", "start", "end")
true_GRanges = makeGRangesFromDataFrame(merged_df, keep.extra.columns = T)  
myannotation@ranges = true_GRanges
names(myannotation@ranges) = Cpg_IDs

dmrcoutput <- dmrcate(myannotation, lambda=1000, C=2, betacutoff = 0.2)

results.ranges <- extractRanges(dmrcoutput, genome = "hg38")

results.ranges



#####################################################visualization############################################################

library(Gviz)

# setting up a variable for grouping and color

pal <- brewer.pal(8,"Dark2")
groups <- pal[1:length(unique(sample_type))]
names(groups) <- levels(sample_type)

visualize_dmr = function(dmrIndex, gen, Granges, file){
  
  cat(dmrIndex, "\n")

  # coordinates are stored under results.ranges[dmrIndex]
  chrom <- as.character(seqnames(results.ranges[dmrIndex]))
  start <- as.numeric(start(results.ranges[dmrIndex]))
  end <- as.numeric(end(results.ranges[dmrIndex]))
  
  # add 25% extra space to plot
  minbase <- start - (0.25*(end-start))
  maxbase <- end + (0.25*(end-start))
  
  
  # defining CpG islands track
  islandTrack = cpgIslands_UCSC(gen, chrom, start, end, title="CpG Islands UCSC")
  
  #Setting up the ideogram, genome, and RefSeq tracks 
  
  iTrack <- IdeogramTrack(genome = gen, chromosome = chrom, name=paste0(chrom))
  gTrack <- GenomeAxisTrack(col="black", cex=1, name="", fontcolor="black")
  biomTrack <- BiomartGeneRegionTrack(genome = gen,
                                      chromosome = chrom, start = minbase, end = maxbase,
                                      name = "ENSEMBL Genes", collapseTranscripts = "longest")
  
  #Ensure that the methylation data is ordered by chromosome and base position.
  bvalOrd <- met[names(myannotation@ranges),]
  
  #Create the data tracks:
  # create genomic ranges object from methylation data
  cpgData <- GRanges(myannotation@ranges, beta = bvalOrd)
  
  # methylation data track
  methTrack <- DataTrack(range=cpgData, 
                         groups=sample_type,
                         genome = gen,
                         chromosome=chrom,
                         ylim=c(-0.05,1.05),
                         col=pal,
                         type=c("a","p"), 
                         name="DNA Meth.\n(beta value)",
                         background.panel="white", 
                         legend=TRUE, 
                         cex.title=0.8,
                         cex.axis=0.8, 
                         cex.legend=0.8)
  
  # DMR position data track
  dmrTrack <- AnnotationTrack(start=start, end=end, genome=gen, name="DMR", 
                              chromosome=chrom,fill="darkred")
  
  
  # Set up the tracklist and indicate the relative sizes of the different tracks. 
  # Finally, draw the plot using the plotTracks function
  tracks <- list(iTrack, gTrack, methTrack, dmrTrack, islandTrack,biomTrack)
  sizes <- c(2,2,7,1,1,2) # set up the relative sizes of the tracks
  # to save figure and scaling graphic device
  pdf(file)
  plotTracks(tracks, from=minbase, to=maxbase, showTitle=TRUE, add53=TRUE, 
             add35=TRUE, grid=TRUE, lty.grid=3, sizes = sizes, length(tracks), transcriptAnnotation = "symbol")
  dev.off()
  
}

for(i in c(2,3,5)){
  visualize_dmr(i,"hg38", results.ranges, paste0("plots/DMR_new", i, ".pdf"))
}




visualize_cnv = function(cnvIndex, cnvRanges, dmrRanges, gen, file){
  cat(cnvIndex, "\n")


  # coordinates are stored under results.ranges[cnvIndex]
  chrom <- as.character(seqnames(cnvRanges[cnvIndex]))
  start <- as.numeric(start(cnvRanges[cnvIndex]))
  end <- as.numeric(end(cnvRanges[cnvIndex]))
  
  # add 25% extra space to plot
  minbase <- start - (0.25*(end-start))
  maxbase <- end + (0.25*(end-start))
  
  
  # defining CpG islands track
  islandTrack = cpgIslands_UCSC(gen, chrom, start, end, title="CpG Islands UCSC")
  
  #Setting up the ideogram, genome, and RefSeq tracks 
  
  iTrack <- IdeogramTrack(genome = gen, chromosome = chrom, name=paste0(chrom))
  gTrack <- GenomeAxisTrack(col="black", cex=1, name="", fontcolor="black")
  biomTrack <- BiomartGeneRegionTrack(genome = gen,
                                      chromosome = chrom, start = minbase, end = maxbase,
                                      name = "ENSEMBL Genes", collapseTranscripts = "longest")
  
  
  overlaps = subsetByOverlaps(dmrRanges, cnvRanges[cnvIndex])

  # DMR position data track
  dmrTrack <- AnnotationTrack(overlaps, genome=gen, name="DMR", 
                              chromosome=chrom,fill="darkred")
  
  
  # Set up the tracklist and indicate the relative sizes of the different tracks. 
  # Finally, draw the plot using the plotTracks function
  tracks <- list(iTrack, gTrack, dmrTrack, islandTrack,biomTrack)
  sizes <- c(2,2,1,1,2) # set up the relative sizes of the tracks
  # to save figure and scaling graphic device
  pdf(file)
  plotTracks(tracks, from=minbase, to=maxbase, showTitle=TRUE, add53=TRUE, 
             add35=TRUE, grid=TRUE, lty.grid=3, sizes = sizes, length(tracks), transcriptAnnotation = "symbol")
  dev.off()
  
}


##get small cnvs
widths = as.data.frame(gistic_GR@ranges)
small_cnvs = head(row.names(widths[order(widths$width),]))


visualize_cnv(34, gistic_GR, results.ranges, "hg38", "plots/CNV34.pdf")

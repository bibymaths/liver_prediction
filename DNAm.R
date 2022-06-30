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

myannotation <- cpg.annotate("array", mval, what = "M", arraytype = "450K",
                             analysis.type="differential", design=design, coef=ncol(design), fdr = 0.01)

dmrcoutput <- dmrcate(myannotation, lambda=1000, C=2, betacutoff = 0.2)

results.ranges <- extractRanges(dmrcoutput, genome = "hg38")
DMR.plot(ranges=results.ranges, dmr=1, CpGs=as.matrix(met), what="Beta",
         arraytype = "450", phen.col=clinical$sample_type, genome="hg38")

results.ranges



#####################################################visualization############################################################

library(Gviz)
dmr.table <- data.frame(results.ranges)

# setting up a variable for grouping and color

pal <- brewer.pal(8,"Dark2")
groups <- pal[1:length(unique(sample_type))]
names(groups) <- levels(sample_type)

#setting up the genomic region 
gen <- "hg38"
# the index of the DMR that we will plot 
dmrIndex <- 1
# coordinates are stored under results.ranges[dmrIndex]
chrom <- as.character(seqnames(results.ranges[dmrIndex]))
start <- as.numeric(start(results.ranges[dmrIndex]))
end <- as.numeric(end(results.ranges[dmrIndex]))

# add 25% extra space to plot
minbase <- start - (0.25*(end-start))
maxbase <- end + (0.25*(end-start))


# defining CpG islands track
cpg <- read.delim("data/cpgIslandExt.txt.gz",header=F, sep = "\t") 
cpg = cpg[cpg$V2==chrom,]

islandData <- GRanges(seqnames=Rle(cpg[,2]), 
                      ranges=IRanges(start=cpg[,3],
                                     end=cpg[,4]),
                      strand=Rle(strand(rep("*",nrow(cpg)))))



#Setting up the ideogram, genome, and RefSeq tracks 

iTrack <- IdeogramTrack(genome = gen, chromosome = chrom, name=paste0(chrom))
gTrack <- GenomeAxisTrack(col="black", cex=1, name="", fontcolor="black")
rTrack <- UcscTrack(genome=gen, chromosome=chrom, track="NCBI RefSeq", 
                    from=minbase, to=maxbase, trackType="GeneRegionTrack", 
                    rstarts="exonStarts", rends="exonEnds", gene="name", symbol="name2", 
                    strand="strand", name="RefSeq", showId=TRUE, geneSymbol=TRUE)

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

# CpG island track
islandTrack <- AnnotationTrack(range=islandData, genome=gen, name="CpG island", 
                               chromosome=chrom,fill="darkgreen")

# DMR position data track
dmrTrack <- AnnotationTrack(start=start, end=end, genome=gen, name="DMR", 
                            chromosome=chrom,fill="darkred")


# Set up the tracklist and indicate the relative sizes of the different tracks. 
# Finally, draw the plot using the plotTracks function
tracks <- list(iTrack, gTrack, methTrack, dmrTrack, islandTrack,rTrack)
sizes <- c(2,2,7,2,2,1) # set up the relative sizes of the tracks
# to save figure and scaling graphic device
pdf("plots/dmr.pdf")
plotTracks(tracks, from=minbase, to=maxbase, showTitle=TRUE, add53=TRUE, 
           add35=TRUE, grid=TRUE, lty.grid=3, sizes = sizes, length(tracks))
dev.off()

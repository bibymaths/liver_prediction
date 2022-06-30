library(curatedTCGAData)
library(TCGAutils)
library(CNVRanger)

lihc <- curatedTCGAData::curatedTCGAData("LIHC",
                                        assays=c("GISTIC_Peaks", "CNVSNP", "RNASeq2GeneNorm"),
                                        version = "1.1.38",
                                        dry.run=FALSE)
lihc

lihc <- TCGAutils::symbolsToRanges(lihc, unmapped=FALSE)

chr1_22 <- c( "chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9",
              "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18",
              "chr19", "chr20", "chr21", "chr22")

for(i in 1:3) 
{
  rr <- rowRanges(lihc[[i]])
  GenomeInfoDb::genome(rr) <- "NCBI37"
  GenomeInfoDb::seqlevelsStyle(rr) <- "UCSC"
  ind <- as.character(seqnames(rr)) %in% chr1_22
  rowRanges(lihc[[i]]) <- rr
  lihc[[i]] <- lihc[[i]][ind,]
}

lihc <- MultiAssayExperiment::intersectColumns(lihc)
TCGAutils::sampleTables(lihc)

lihc <- TCGAutils::TCGAsplitAssays(lihc, sampleCodes="01")

ind <- grep("CNVSNP", names(lihc))
head( mcols(lihc[[ind]]) )

summary( mcols(lihc[[ind]])$Segment_Mean )

smean <- mcols(lihc[[ind]])$Segment_Mean
state <- round(2^smean * 2)
state[state > 4] <- 4
mcols(lihc[[ind]])$state <- state
lihc[[ind]] <- lihc[[ind]][state != 2,]
mcols(lihc[[ind]]) <- mcols(lihc[[ind]])[,3:1]
table(mcols(lihc[[ind]])$state)

res <- cnvEQTL(cnvrs="01_LIHC_GISTIC_Peaks-20160128", 
               calls="01_LIHC_CNVSNP-20160128", 
               rcounts="01_LIHC_RNASeq2GeneNorm-20160128_ranged", 
               data=lihc, window="1Mbp", verbose=TRUE)

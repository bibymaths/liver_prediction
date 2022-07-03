library(curatedTCGAData)
library(tidyverse)
library(GenomicRanges)
library(Gviz)

lihc_gistic <- curatedTCGAData::curatedTCGAData("LIHC",
                                         assays=("GISTIC_Peaks"),
                                         version = "1.1.38",
                                         dry.run=FALSE)

gistic <- as.data.frame(assay(lihc_gistic))

gistic <- tibble::rownames_to_column(gistic, "position")

gistic <- gistic %>%
  separate(col = position, into = c("seqnames", "Position"), sep = ":") %>%
  separate(col = Position, into = c("Start", "End"), sep = "-")

gistic_GR <- makeGRangesFromDataFrame(gistic)
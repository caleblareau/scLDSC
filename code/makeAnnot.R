# Script takes raw annotation file and some SE object to create 
# annotation files for LD Score Regression

library(GenomicRanges)
library(data.table)
library(SummarizedExperiment)
library(Matrix)
library(diffloop)
BiocParallel::register(BiocParallel::MulticoreParam(2, progressbar = FALSE))

if (basename(getwd()) != "code") setwd("code")

# Load human heme
load("/Volumes/dat/Research/BuenrostroResearch/hemeGWASanno/output/peaksCounts.rda")
types <- readRDS("/Volumes/dat/Research/BuenrostroResearch/hemeGWASanno/output/sampleTypes.rds")
sample_annotation <- DataFrame(type = types)
hemese <- SummarizedExperiment(assays = list(counts = Matrix(counts)), rowRanges = peaks, colData = sample_annotation)


makeAnnotFile <- function(chr, outname, se, mincounts = 1){

  rawfile <- paste0("../raw.annot/raw.", as.character(chr), ".annot.gz")
  
  # Import raw annotation file
  df <- data.frame(fread(input = paste0('zcat < ', rawfile)))
  gdf <- makeGRangesFromDataFrame(data.frame(chr = df$CHR, start = df$BP, end = df$BP))
  
  # Deal with SE oject
  boo <- which(as.character(seqnames(rmchr(rowRanges(se)))) == chr)
  rr <- rmchr(rowRanges(se)[boo])
  counts <- assays(se)[["counts"]][boo,]
  ov <- data.frame(findOverlaps(rr, gdf))
  
  # Boolean Matrix 
  booCountsMat <- counts[ov$queryHits,] >= mincounts
  applyOut <- sapply(1:dim(booCountsMat)[2], function(i){
     as.numeric(1:dim(df)[1] %in% ov$subjectHits[booCountsMat[i,]])
  })
  
  dat <- data.frame(df, applyOut)
  names(dat) <- c(colnames(df), dimnames(counts)[[2]])
  
  gz1 <- gzfile(paste0(outname, as.character(chr), ".annot.gz"), "w")
  write.table(dat, file = gz1, row.names = FALSE, col.names = TRUE, sep = "\t")
  close(gz1)
  return(chr)
}

makeAnnotFile(21, "heme", hemese)
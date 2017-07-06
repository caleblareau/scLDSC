# Script takes raw annotation file and some .bed files and makes
# annotation files for LD Score Regression

# One note that's kinda nice I guess: whether 'chr' is on the
# bed file or not shouldn't... matter (I scrub all the files and remove 
# it beforehand but it shouldn't affect the output)

library(GenomicRanges)
library(data.table)
library(diffloop)

BiocParallel::register(BiocParallel::MulticoreParam(2, progressbar = FALSE))

if (basename(getwd()) != "code") setwd("code")

# Arguments: chr <-> chromosome (each one gets its own file)
# bedfile <-> duh; but is it gzipped or nah? (see below)
# outname <-> prefix name for the output file
# gzipIn <-> whether or not the .bed file is gzipped going in


# Function takes 1 bed file and 1 chromosome for a specific
makeAnnotFileFromBed <- function(chr, bedfile, outname = NULL, gzipIn = TRUE){

  # Figure out outname if not distinctly specified
  if(is.null(outname) & gzipIn) outname <- gsub('.{7}$', '',basename(bedfile))
  if(is.null(outname) & !gzipIn) outname <- gsub('.{3}$', '',basename(bedfile))
  
  # Import raw annotation file
  rawfile <- paste0("../raw.annot/raw.", as.character(chr), ".annot.gz")  # change this !!!!!!!!!!!!
  df <- data.frame(fread(input = paste0('zcat < ', rawfile)))
  gdf <- makeGRangesFromDataFrame(data.frame(chr = df$CHR, start = df$BP, end = df$BP))
  
  # Import bedfile
  if(gzipIn){  beddf <- data.frame(fread(input = paste0('zcat < ', bedfile)))
  } else {  beddf <- data.frame(fread(input = paste0('', bedfile))) } 
  gbed <- makeGRangesFromDataFrame(data.frame(chr = beddf$V1, start = beddf$V2, end = beddf$V3))
  
  # Find overlapping hits and stuff
  boo <- as.numeric(1:length(gdf) %in% queryHits(findOverlaps(rmchr(gdf), rmchr(gbed))))
  
  dat <- data.frame(df, boo)
  names(dat) <- c(colnames(df), outname)

  gz1 <- gzfile(paste0("../1000G_annotation/", outname, as.character(chr), ".annot.gz"), "w") # change this !!!!!!!!!!!!
  write.table(dat, file = gz1, row.names = FALSE, col.names = TRUE, sep = " ", quote = FALSE)
  close(gz1)
  return(chr)
}

baselineBeds <- list.files("../baseline", full.names = TRUE)

sapply(baselineBeds[1:2], function(bedfile){
  sapply(1:22, function(i){
    makeAnnotFileFromBed(i, bedfile)
    paste0(bedfile, "_", as.character(i)) # dummy thing to return
  })
})
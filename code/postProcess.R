# Script takes ldsc output + scMyeloid trajectory and does stuff

library(GenomicRanges)
library(data.table)
library(SummarizedExperiment)
library(Matrix)
#library(diffloop)
#BiocParallel::register(BiocParallel::MulticoreParam(8, progressbar = FALSE))

if (basename(getwd()) != "code") setwd("code")

cells <- read.table("../dataOrder.txt", header = TRUE)
ps <-read.table("/data/aryee/caleb/ldsc/scMyeloid/22jan2017_myel_pseudotime.txt")
m <- merge(cells, ps, by.x = "name",  by.y = "V1", sort = FALSE)
names(m) <- c("name", "type", "pseudotime")
ldscout <- read.table("/data/aryee/caleb/ldsc/scMyeloid/BMI.results", header =TRUE)
baseline <- ldscout[1,]
ldscout <- ldscout[-1,]

m <- data.frame(m, ldscout)

cols <- c("#FFC179", "#FFA300", "#AFAFAF", "#7D7D7D", "#4B4B4B","#00441B", "#00AF99", "#FF5A00", "#46A040", "#C390D4")
ggplot(m, aes(x = pseudotime, y = Enrichment, color = type)) + geom_point() + theme_bw() +  
  labs(list(title = "BMI", x = "Pseudotime", y = "LDSC Enrichment")) + scale_colour_manual(values = cols)

ggplot(m, aes(x = pseudotime, y = Prop._SNPs, color = type)) + geom_point() + theme_bw() +  
  labs(list(title = "Proportion of Variants", x = "Pseudotime", y = "LDSC Enrichment")) + scale_colour_manual(values = cols)


means <- tapply(m$Enrichment, m$type, mean)


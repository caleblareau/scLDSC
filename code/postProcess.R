# Script takes ldsc output + scMyeloid trajectory and does stuff

library(GenomicRanges)
library(data.table)
library(SummarizedExperiment)
library(Matrix)

if (basename(getwd()) != "code") setwd("code")

cells <- read.table("../dataOrder.txt", header = TRUE)
ps <-read.table("/data/aryee/caleb/ldsc/scMyeloid/22jan2017_myel_pseudotime.txt")
names(ps) <- c("name", "pseudotime")
oldGWAS <- read.table("/data/aryee/caleb/evolution_conservation/scHeme/GWASdeviations.txt", header = TRUE)
samples <- read.table("/data/aryee/caleb/evolution_conservation/scHeme/data/scATAC/samples.txt", header = TRUE)
df <- data.frame(samples, oldGWAS)

pdf <- merge(ps, df, by = c("name"))
#mdf <- melt(pdf, id.vars = c("name", "pseudotime", "type"))

m2 <- merge(cells, pdf, by.x = "name",  by.y = "name", sort = FALSE)
ldscout <- read.table("/data/aryee/caleb/ldsc/scMyeloid/RA.results", header =TRUE)
baseline <- ldscout[1,]
ldscout <- ldscout[-1,]

m <- data.frame(m2, ldscout)
#m <- m[m$Prop._SNPs > 0.001,]

cols <- c("#FFC179", "#FFA300", "#AFAFAF", "#7D7D7D", "#4B4B4B","#00441B", "#00AF99", "#FF5A00", "#46A040", "#C390D4")
ggplot(m, aes(x = pseudotime, y = Enrichment, color = type.x)) + geom_point() + theme_bw() +  
  labs(list(title = "RA", x = "Pseudotime", y = "LDSC")) + scale_colour_manual(values = cols)

ggplot(m, aes(y = Enrichment, x = Prop._SNPs, color = type)) + geom_point() + theme_bw() +  
  labs(list(title = "Proportion of Variants", x = "Prop of SNPs", y = "LDSC Enrichment")) + scale_colour_manual(values = cols)


ggplot(m, aes(x = ra, y  = Enrichment, color = type.x)) + geom_point() + theme_bw() +  
  labs(list(title = "LDSC vrersus chromVAR", x = "chromVAR", y = "LDSC Enrichment")) + scale_colour_manual(values = cols)


tapply(m$Enrichment, m$type, mean)


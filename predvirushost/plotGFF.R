#!/usr/bin/env Rscript
list.of.packages <- c("tidyverse", "genoPlotR")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages) > 0) install.packages(new.packages)
library(tidyverse)

library(genoPlotR)
args <- commandArgs(trailingOnly = T)
gff_file <- args[1]
gff <- read.table(gff_file, sep = "\t", fill = T, comment.char = "#", quote = "")
colnames(gff) <- c("sequence", "source", "feature", "start", "end", "score", "strand", "phase", "atrribute")
gff <- gff %>% separate(col = atrribute, into = c("call", "description"), remove = F, extra = "merge", sep = ":")
gff <- gff %>% mutate(unique_names = paste(sequence, row_number(), sep = "_"))
gff$call[is.na(gff$call)] <- "none"
names1 <- gff$unique_names
starts1 <- gff$start
ends1 <- gff$end
strands1 <- gff$strand
gff <- gff %>% mutate(colours = ifelse(strand == "+", "grey", "black"))
gff <- gff %>% mutate(colours = ifelse(call == "archaeal", "blue", ifelse(call == "phage", "red", ifelse(call == "eukaryotic", "green", ifelse(call == "none", "black", colours)))))
cols1 <-gff$colours
df1 <- data.frame(name=names1, start=starts1, end=ends1,
                  strand=strands1, col="black", fill = cols1)
dna_seg1 <- dna_seg(df1)
is.dna_seg(dna_seg1)

svg(filename= paste(gff_file, ".svg", sep = ""), 
    width=30, 
    height=15, 
    pointsize=22)
plot_gene_map(dna_segs=list(dna_seg1))

dev.off()
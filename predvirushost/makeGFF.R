#!/usr/bin/env Rscript
list.of.packages <- c("tidyverse", "getopt")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages) > 0) install.packages(new.packages)
options(warn = -1)
library('getopt')
library(tidyverse)


spec = matrix(c(
  'hmmres', 'i', 1, "character",
  'gff', 'g', 1, "character",
  'genome_lookup', 'l', 1, "character",
  'help' , 'h', 0, "logical",
  'file_type', 'f', 2, "character",
  'number_of_genomes', 'n', 2, "character",
  'path_1', 'p', 2, "character",
  'path_2', 'q', 2, "character",
  'contig_name', 'c', 2, "character"
), byrow=TRUE, ncol=4)


opt = getopt(spec)

if ( !is.null(opt$help) ) {
  cat("cmscanToGffWrapper.R version 1.0\n\n")
  cat("Use cmscanToGffWrapper.R <options> -f <cmscan ouptut file> -g <gff file>\n\n")
  cat("Options:\n")
  cat("  -f <cmscan ouptut file> The file that contains the cmscan output\n")
  cat("  -g <gff file> The file that contains the gff data. Do not inclue the gff file extension\n")
  cat("  -f <file path> The location of the other files and the output file\n")
  cat("  -o <output file name> The name of the output file. Do not inclue the gff file extension. The default is the same as the gca input\n")
  q(status=1)
}

# if ( is.null(opt$cmscanOutput) ) {
#   cat("Error: -f <cmscan ouptut file> is required.\n")
#   q(status=1)
# }
# if ( is.null(opt$gcf) ) {
#   cat("Error: -g <gff file> is required.\n")
#   q(status=1)
# }


test <- F
if(test){
  opt$path_1 <- "~/bin/PredVirusHost/predvirushost"
  opt$path_2 <- "~/bin/PredVirusHost/archaea/"
  opt$file_type <- "REFSEQ"
  opt$contig_name <- "Methanothermobacter_phage_psiM100]"
  opt$gff <- "~/bin/PredVirusHost/archaea/gff_files/Methanothermobacter_phage_psiM100.gff3"
  opt$hmmres <- "~/bin/PredVirusHost/archaea/archaea_all_2.proteins.txt"
  opt$genome_lookup <- "~/bin/PredVirusHost/archaea/archaea_all.genome_lookup.txt"
}


filePath <- opt$path_1
filePath2 <- opt$path_2
print(filePath)
input_type <- opt$file_type
gff_used <- F


if(is.null(opt$number_of_genomes)){ opt$number_of_genomes <- 3}

if(!is.null(opt$contig_name)){
  contig_name <- opt$contig_name
  number_of_genomes <- 1
}else{
  number_of_genomes <- opt$number_of_genomes
}


genome_lookup <- read.table(opt$genome_lookup, sep = " ", header = F, comment.char = "", quote = "", fill = T)
genome_lookup <- genome_lookup%>%dplyr::rename(target.name = V2, genome = V1)
if(opt$gff != "None"){
  gffIn <- read.table(opt$gff, comment.char = "#", quote = "", as.is = T, header = F, sep = "\t", fill = T)
  colnames(gffIn) <- c("sequence", "source", "feature", "start", "end", "score", "strand", "phase", "atrribute")
  gff <- gffIn %>% filter(feature == "CDS")
  fasta <- readLines(paste(filePath2, "/tmp.faa", sep = ""))
  fasta <- data.frame(tmp = fasta)
  
  fasta <- fasta %>% separate(col = tmp, into = c("target.name", "sequence"), sep = "\t", extra = "merge")
  
  sequences <- fasta
  sequences$length <- str_length(fasta$sequence)*3
  
  
  sequences <- sequences %>% filter(!is.na(sequence)) %>% 
    mutate(left = 1) %>% 
    mutate(length.total = cumsum(length)) %>% 
    mutate(right =  length.total) %>% 
    mutate(left = right - length + 1) %>% 
    select(target.name, left, right)
  
  sequences <- sequences %>% left_join(genome_lookup, by = "target.name") 
  sequences <- sequences %>% filter(contig_name == genome)
  
  # sequences <- sequences %>% mutate(gene.length = right - left + 3)
  # gff <- gff %>% mutate(gene.length = end - start)
  
  
  if(nrow(sequences) == nrow(gff)){
    gff <- gff %>% bind_cols(sequences)
    gff_used <- T
  }else{
    print(nrow(sequences))
    print(nrow(gff))
    cat("not using gff\n")
  }
}else{
 fasta <- readLines(paste(filePath2, "/tmp.faa", sep = ""))
  fasta <- data.frame(tmp = fasta)
  
  fasta <- fasta %>% separate(col = tmp, into = c("target.name", "sequence"), sep = "\t", extra = "merge")

  sequences <- fasta
  sequences$length <- str_length(fasta$sequence)*3
  
  
  sequences <- sequences %>% filter(!is.na(sequence)) %>% 
    mutate(left = 1) %>% 
    mutate(length.total = cumsum(length)) %>% 
    mutate(right =  length.total) %>% 
    mutate(left = right - length + 1) %>% 
    select(target.name, left, right)
  
  sequences <- sequences %>% left_join(genome_lookup, by = "target.name") 
  

}



hmmDat <- read.table(opt$hmmres, sep = "\t", as.is = T, header = T)

hmmDat <- hmmDat %>% group_by(target.name) %>% arrange(-score) %>% mutate(id = row_number()) %>% filter(id == 1) %>% select(-id)

if(gff_used){
  cat("gff used\n")
  hmmDat <- hmmDat %>% filter(grepl(pattern = contig_name, x = genome))

  hmmDat <- hmmDat %>% select(target.name, protein_descriptions, score, genera, genomes, call) %>% 
    rename(hmmscore = score)
  
  
  gffTmp <- gff %>% left_join(hmmDat, by = "target.name")
  
  gffTmp$call[is.na(gffTmp$call)] <- "none"
  
  gff <- data.frame(seqname = gffTmp$sequence,
                    source = gffTmp$source,
                    feature = "CDS",
                    start = gffTmp$start,
                    end = gffTmp$end,
                    score = gffTmp$hmmscore,
                    strand = gffTmp$strand,
                    frame = gffTmp$phase,
                    attribute = gffTmp$atrribute,
                    attribute1 = gffTmp$protein_descriptions, 
                    attribute2 = gffTmp$genomes,
                    attribute3 = gffTmp$genera, 
                    attribute4 = gffTmp$call,
                    target.name = gffTmp$target.name)

  gff <- gff %>% mutate(attribute = paste(attribute4, attribute1, attribute2, attribute3, attribute, sep = ":"))%>% 
    select(seqname, source, feature, start, end, score, strand, frame, attribute) 
  
  gffLines <- readLines(opt$gff)
  gffLines <- data.frame(x = gffLines)
  gffLines <- gffLines %>% filter(grepl(pattern = "#", x))
  
  fileConn<-file(paste(filePath2, "/", contig_name, "_contig.gff", sep = ""))
  writeLines(c("##gff-version 3",
               "#!gff-spec-version 1.21",
               "#!processor R script (local)",
               "#!Contigs_from_predvirushost",
               paste("#!genome-build-accession NCBI_Assembly:", "predvirushost_name", sep = ""),
               paste("#!annotation-date ", Sys.Date(), sep = ""),
               "#!annotation-source hmmsearch (VOGs) (local version)",
               as.character(gffLines$x[nrow(gffLines)])), fileConn)
  close(fileConn)
  
  cat(paste(filePath2, "/", contig_name, "_contig.gff\n", sep = ""))
  
  
  write.table(x = gff, file = paste(filePath2, "/", contig_name, "_contig.gff", sep = ""), row.names = F, col.names = F, quote = F, sep = "\t", append = T)
  
  
  
  
  }else{
  cat("gff not used\n")
  gff <- data.frame(seqname = hmmDat$genome,
                    source = "predvirushost",
                    feature = "CDS",
                    start = 1,
                    end = 100,
                    score = hmmDat$score,
                    strand = "+",
                    frame = ".",
                    attribute1 = hmmDat$protein_descriptions, 
                    attribute2 = hmmDat$genomes,
                    attribute3 = hmmDat$genera, 
                    attribute4 = hmmDat$call,
                    target.name = hmmDat$target.name)
  
  
  gff <- gff %>% mutate(attribute = paste(attribute4, attribute1, attribute2, attribute3, sep = ":"))
  
  gff <- gff %>% full_join(sequences, by ="target.name")
  
  gff <- gff %>% mutate(start = left, end = right) %>% 
    arrange(start) %>% 
    mutate( source = "predvirushost",
            feature = "CDS",
            score = ifelse(is.na(score), 0, as.numeric(score)),
            strand = ifelse(is.na(strand), "+", as.character(strand)),
            frame = ".",
            seqname = ifelse(is.na(seqname), as.character(genome), as.character(seqname))) %>% 
    select(seqname, source, feature, start, end, score, strand, frame, attribute) 
  
  
  genomes_ordered <- hmmDat %>% arrange(-count) %>% mutate(genome = as.character(genome))
  
  
  
  if(!is.null(opt$contig_name)){
    genome_selected <- contig_name
    dat <- gff %>% filter(seqname == genome_selected)
    location_shift <- dat$start[1] - 1
    dat$start <- dat$start - location_shift
    dat$end <- dat$end - location_shift
    
    
    fileConn<-file(paste(filePath2, "/", dat$seqname[1], "_contig.gff", sep = ""))
    writeLines(c("##gff-version 3",
                 "#!gff-spec-version 1.21",
                 "#!processor R script (local)",
                 "#!Contigs_from_predvirushost",
                 paste("#!genome-build-accession NCBI_Assembly:", "predvirushost_name", sep = ""),
                 paste("#!annotation-date ", Sys.Date(), sep = ""),
                 "#!annotation-source hmmsearch (VOGs) (local version)",
                 paste("##sequence-region ", dat$seqname[1], " ", dat$start[1], " ", dat$end[nrow(dat)], sep = "")), fileConn)
    close(fileConn)
    
    cat(paste(filePath2, "/", dat$seqname[1], "_contig.gff\n", sep = ""))
    
    write.table(x = dat, file = paste(filePath2, "/", dat$seqname[1], "_contig.gff", sep = ""), row.names = F, col.names = F, quote = F, sep = "\t", append = T)
    
  }else{
    for(i in 1:number_of_genomes){
      genome_selected <- unique(genomes_ordered$genome)[i]
      dat <- gff %>% filter(seqname == genome_selected)
      location_shift <- dat$start[1] - 1
      dat$start <- dat$start - location_shift
      dat$end <- dat$end - location_shift
      
      
      fileConn<-file(paste(filePath2, "/", dat$seqname[1], "_contig.gff", sep = ""))
      writeLines(c("##gff-version 3",
                   "#!gff-spec-version 1.21",
                   "#!processor R script (local)",
                   "#!Contigs_from_predvirushost",
                   paste("#!genome-build-accession NCBI_Assembly:", "predvirushost_name", sep = ""),
                   paste("#!annotation-date ", Sys.Date(), sep = ""),
                   "#!annotation-source hmmsearch (VOGs) (local version)",
                   paste("##sequence-region ", dat$seqname[1], " ", dat$start[1], " ", dat$end[nrow(dat)], sep = "")), fileConn)
      close(fileConn)
      
      cat(paste(filePath2, "/", dat$seqname[1], "_contig.gff\n", sep = ""))
      
      write.table(x = dat, file = paste(filePath2, "/", dat$seqname[1], "_contig.gff", sep = ""), row.names = F, col.names = F, quote = F, sep = "\t", append = T)
      
    }
  }
  
}





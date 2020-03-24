#!/usr/bin/env Rscript
options(warn = -1)

suppressMessages(library(dplyr, quietly = T))


args <- commandArgs(trailingOnly = T)
filePath <- args[1]
filePath2 <- args[2]

if(args[3] == "y"){
discriminant_models_only <- T
#cat("Using discriminant models\n")
}else{
#  cat("Using all models\n")
  discriminant_models_only <- F
}

#print(filePath)

aa <- read.table(paste(filePath2, "/arVOG_res.txt", sep = ""), fill = T)
pp <- read.table(paste(filePath2, "/baPOG_res.txt", sep = ""), fill = T)
ee <- read.table(paste(filePath2, "/euVOG_res.txt", sep = ""), fill = T)

aa <- aa[,1:19]
pp <- pp[,1:19]
ee <- ee[,1:19]

hmmColnames <- c("target.name", "accession", "query.name", "accession.2", "E.value", "score", "bias", "E.value.2", "score.2", "bias.2", "exp", "reg", "clu", "ov", "env", "dom", "rep", "inc", "description.of.target")
colnames(aa) <- hmmColnames
colnames(pp) <- hmmColnames
colnames(ee) <- hmmColnames

genome_lookup <- read.table(paste(filePath2, "/genome_lookup.txt", sep = ""), sep = " ", header = F, comment.char = "", quote = "", fill = T)
genome_lookup <- genome_lookup%>%dplyr::rename(target.name = V2, genome = V1)

archaealModels <- read.table(paste(filePath, "/archaeal_model_scores.txt", sep = ""), sep = "\t", comment.char = "", quote = "", fill = T, as.is = T, header = T)
phageModels <- read.table(paste(filePath, "/phage_model_scores.txt", sep = ""), sep = "\t", comment.char = "", quote = "", fill = T, as.is = T, header = T)
eukaryoticModels <- read.table(paste(filePath, "/eukaryotic_model_scores.txt", sep = ""), sep = "\t", comment.char = "", quote = "", fill = T, as.is = T, header = T)


archaealDescr <- read.table(paste(filePath, "/archaeal_hmm_accessions_descriptions.csv", sep = ""), sep = "\t", comment.char = "", quote = "", fill = T, as.is = T, header = T)
phageDescr <- read.table(paste(filePath, "/phage_hmm_accessions_descriptions.csv", sep = ""), sep = "\t", comment.char = "", quote = "", fill = T, as.is = T, header = T)
eukaryoticDesc <- read.table(paste(filePath, "/eukaryotic_hmm_accessions_descriptions.csv", sep = ""), sep = "\t", comment.char = "", quote = "", fill = T, as.is = T, header = T)

# archaealDescr <- read.csv("~/bin/PredVirusHost/archaeal_hmm_accessions_descriptions.csv", comment.char = "", quote = "", fill = T, as.is = T, header = T)
# phageDescr <- read.csv("~/bin/PredVirusHost/phage_hmm_accessions_descriptions.csv", comment.char = "", quote = "", fill = T, as.is = T, header = T)
# eukaryoticDescr <- read.csv("~/bin/PredVirusHost/eukaryotic_hmm_accessions_descriptions.csv", comment.char = "", quote = "", fill = T, as.is = T, header = T)

eukaryoticDescr <- eukaryoticDescr %>% mutate(genera = NA)

allDescr <- archaealDescr %>% bind_rows(phageDescr, eukaryoticDescr)

colnames(allDescr)[1] <- "query.name"

aa <- aa%>%left_join(genome_lookup, by = "target.name")%>%filter(as.numeric(score) > 30)
pp <- pp%>%left_join(genome_lookup, by = "target.name")%>%filter(as.numeric(score) > 30)
ee <- ee%>%left_join(genome_lookup, by = "target.name")%>%filter(as.numeric(score) > 30)


aa <- aa%>%left_join(archaealModels, by = "query.name")
pp <- pp%>%left_join(phageModels, by = "query.name")
ee <- ee%>%left_join(eukaryoticModels, by = "query.name")

aa <- aa%>%left_join(allDescr, by = "query.name")
pp <- pp%>%left_join(allDescr, by = "query.name")
ee <- ee%>%left_join(allDescr, by = "query.name")

if(discriminant_models_only == T){
  cat("Using discriminant models\n")
aa <- aa%>%mutate(score = ifelse(score.percentage < 0.8, 0, score*score.percentage))
pp <- pp%>%mutate(score = ifelse(score.percentage < 0.8, 0, score*score.percentage))
ee <- ee%>%mutate(score = ifelse(score.percentage < 0.8, 0, score*score.percentage))
}else{
  cat("Using all models\n")
aa <- aa%>%mutate(score = score*score.percentage)
pp <- pp%>%mutate(score = score*score.percentage)
ee <- ee%>%mutate(score = score*score.percentage)
}





aaRes <- aa%>%group_by(genome)%>%summarise(archaeal_score = sum(score))%>%
  mutate(archaeal_score = ifelse(is.na(archaeal_score), 0, archaeal_score))


ppRes <- pp%>%group_by(genome)%>%summarise(phage_score = sum(score))%>%
  mutate(phage_score = ifelse(is.na(phage_score), 0, phage_score))


eeRes <- ee%>%group_by(genome)%>%summarise(eukaryotic_score = sum(score))%>%
  mutate(eukaryotic_score = ifelse(is.na(eukaryotic_score), 0, eukaryotic_score))


res <- aaRes%>%full_join(ppRes, by = "genome")%>%full_join(eeRes, by = "genome")
res <- res%>%
  mutate(phage_score = ifelse(is.na(phage_score), 0, phage_score))%>%
  mutate(archaeal_score = ifelse(is.na(archaeal_score), 0, archaeal_score))%>%
  mutate(eukaryotic_score = ifelse(is.na(eukaryotic_score), 0, eukaryotic_score))%>%
  mutate(call = ifelse(archaeal_score >= phage_score,
                       ifelse(archaeal_score >= eukaryotic_score,
                              ifelse(archaeal_score > 60, "archaeal", "none"),
                              ifelse(eukaryotic_score > 150, "eukaryotic", "none")),
                       ifelse(phage_score >= eukaryotic_score, ifelse(phage_score > 110, "phage", "none"), ifelse(eukaryotic_score > 150, "eukaryotic", "none"))))


res <- res%>%
  mutate(difference.percentage = ifelse(call == "archaeal",
                                        ifelse(phage_score >= eukaryotic_score,
                                               round((archaeal_score - phage_score)/archaeal_score*100, 2),
                                               round((archaeal_score - eukaryotic_score)/archaeal_score*100, 2)
                                        ),
                                        ifelse(call == "phage",ifelse(archaeal_score >= eukaryotic_score,
                                                                      round((phage_score - archaeal_score)/phage_score*100, 2),
                                                                      round((phage_score - eukaryotic_score)/phage_score*100, 2)
                                        ), ifelse(archaeal_score >= phage_score,
                                                  round((eukaryotic_score - archaeal_score)/eukaryotic_score*100, 2),
                                                  round((eukaryotic_score - phage_score)/eukaryotic_score*100, 2)
                                        ))
                                        ))#%>%
  #mutate(call = ifelse(archaeal_score == 0, ifelse(phage_score == 0, ifelse(eukaryotic_score == 0, "none", call),call),call))

proteinCounts <- genome_lookup%>%group_by(genome)%>%summarise(protein.counts = n())

res <- res%>%full_join(proteinCounts, by = "genome")
res <- res%>%mutate(call = ifelse(is.na(call), "none", call),
                    archaeal_score = ifelse(is.na(archaeal_score), 0, archaeal_score),
                    phage_score = ifelse(is.na(phage_score), 0, phage_score),
                    eukaryotic_score = ifelse(is.na(eukaryotic_score), 0, eukaryotic_score),
                    difference.percentage = ifelse(is.na(difference.percentage), 0, difference.percentage))



write.table(res, paste(filePath2, "/scores.txt", sep = ""), col.names = T, row.names = F, quote = F, sep = "\t")

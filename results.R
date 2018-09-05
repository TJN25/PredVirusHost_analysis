library(MASS)
x <- read.table('output_with_species.txt', header=TRUE)
levelsX <- levels(x$Contig)
lenX <- length(x$Contig)
resultsE <- c()
for (i in levelsX){
  evalX <- subset(x, Contig==i, select=c(E.value.2))
  lenC <- length(evalX) 
  adjX <- evalX*lenX
  df = 2*length(adjX)
  e <- 1 - pchisq( -2*sum(log(adjX)), df)
  resultsE[i] <- e
}

resultsB <- c()
for (i in levelsX){
  evalX <- subset(x, Contig==i, select=c(E.value.1))
  lenC <- length(evalX) 
  adjX <- evalX*lenX
  df = 2*length(adjX)
  b <- 1 - pchisq( -2*sum(log(adjX)), df)
  resultsB[i] <- b
}

resultsA <- c()
for (i in levelsX){
  evalX <- subset(x, Contig==i, select=c(E.value))
  lenC <- length(evalX) 
  adjX <- evalX*lenX
  df = 2*length(adjX)
  a <- 1 - pchisq( -2*sum(log(adjX)), df)
  resultsA[i] <- a
}

# Takes the adjusted e-values and for each of the contigs a score for each domain is calculated and printed
# The printing is done in the form contig  arc	bac	euk
dataa <- read.table('scores_calc.txt', header=TRUE, sep='	')
datab <- levels(dataa$Species)
arcmean <- c()
arcmed <- c()
bacmean <- c()
bacmed <- c()
eukmean <- c()
eukmed <- c()
lenA <- c()
for (i in datab){
  a <- subset(dataa, Species==i, select=c(Archaeal))
  alen <- length(a$Archaeal)
  b <- subset(dataa, Species==i, select=c(Bacterial))
  c <- subset(dataa, Species==i, select=c(Eukaryotic))
  arc <- as.numeric(as.character(a$Archaeal))
  bac <- as.numeric(as.character(b$Bacterial))
  euk <- as.numeric(as.character(c$Eukaryotic))
  arcmean[i] <- mean(arc)
  arcmed[i] <- median(arc)
  bacmean[i] <- mean(bac)
  bacmed[i] <- median(bac)
  eukmean[i] <- mean(euk)
  eukmed[i] <- median(euk)
  lenA[i] <- alen
}




rownamesM <- c(levelsX)
resultsX <- c(rownamesM,resultsA,resultsB,resultsE, arcmean, bacmean, eukmean, arcmed, bacmed, eukmed, lenA)
matX <- matrix(resultsX, ncol=11, byrow=FALSE)
write.matrix(format(matX, scientific=FALSE), file ="output_scores.txt", sep="\t")

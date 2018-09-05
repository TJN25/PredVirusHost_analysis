#------------------------------------Functions---------------------------------------------#

transform.df <- function(dat){
  for(i in 1:length(dat[1,])){
    if(substr(colnames(dat)[i],(nchar(colnames(dat)[i])-3),nchar(colnames(dat)[i]) )=='.num'){
      #print(paste(colnames(dat)[i], 'is numeric', sep = ' '))
      dat[,i] <- as.character(dat[,i])
      dat[,i] <- as.numeric(dat[,i])
    }else if(substr(colnames(dat)[i], nchar(colnames(dat)[i])-3, nchar(colnames(dat)[i]))=='.y.n'){
      # print(paste(colnames(dat)[i], 'is logical', sep = ' '))
      dat[,i] <- as.logical(dat[,i])
      #dat[,i] <- as.character(dat[,i])
    }else{
      #print(i)
      dat[,i] <- as.character(dat[,i])
    }
  }
  return(dat)
}
#------------------------------------Unused Functions---------------------------------------------#



#------------------------------------Setup and import data---------------------------------------------#
#setwd('~/Desktop/Scripts/testing/')
library(MASS)
protein.dat <- read.table('output_with_species.txt', header=F, sep = '\t', as.is = T)
colnames(protein.dat) <- c('protein','a.model','a.e.val.num','b.model','b.e.val.num','e.model','e.e.val.num','contig')
protein.dat <- transform.df(protein.dat)
protein.dat <- protein.dat[order(protein.dat$contig),]
contig.dat <- data.frame(contig = unique(protein.dat$contig),
                         call = vector(mode = 'character', length = length(unique(protein.dat$contig))),
                         archaeal.score.num = vector(mode = 'character', length = length(unique(protein.dat$contig))),
                         bacterial.score.num = vector(mode = 'character', length = length(unique(protein.dat$contig))),
                         eukaryotic.score.num = vector(mode = 'character', length = length(unique(protein.dat$contig))),
                         num.of.proteins.num = vector(mode = 'character', length = length(unique(protein.dat$contig))),
                         num.of.proteins.scored = vector(mode = 'character', length = length(unique(protein.dat$contig))),
                         archaeal.e.vals = vector(mode = 'character', length = length(unique(protein.dat$contig))),
                         bacterial.e.vals = vector(mode = 'character', length = length(unique(protein.dat$contig))),
                         eukaryotic.e.vals = vector(mode = 'character', length = length(unique(protein.dat$contig))),
                         archaeal.models = vector(mode = 'character', length = length(unique(protein.dat$contig))),
                         bacterial.models = vector(mode = 'character', length = length(unique(protein.dat$contig))),
                         eukaryotic.models = vector(mode = 'character', length = length(unique(protein.dat$contig))))
contig.dat <- transform.df(contig.dat)
contig.list <- unique(protein.dat$contig)
for(i in 1:length(unique(protein.dat$contig))){
  #print(i)
  contig <- contig.list[i]
  dat <- protein.dat[protein.dat[,8]==contig,]
  dat <- transform.df(dat)
  #print(length(dat[,1]))
  contig.dat[i,6] <- length(dat[,1])
  contig.dat[i,8] <- paste(dat$a.e.val.num, collapse = ',')
  contig.dat[i,9] <- paste(dat$b.e.val.num, collapse = ',')
  contig.dat[i,10] <- paste(dat$e.e.val.num, collapse = ',')
  contig.dat[i,11] <- paste(dat$a.model, collapse = ',')
  contig.dat[i,12] <- paste(dat$b.model, collapse = ',')
  contig.dat[i,13] <- paste(dat$e.model, collapse = ',')
  score.dat.a <- dat[dat[,3]<1e-5,]
  score.dat.a <- transform.df(score.dat.a)
  score.dat.b <- dat[dat[,5]<1e-5,]
  score.dat.b <- transform.df(score.dat.b)
  score.dat.e <- dat[dat[,7]<1e-5,]
  score.dat.e <- transform.df(score.dat.e)
  a.len <- length(score.dat.a[,1])
  b.len <- length(score.dat.b[,1])
  e.len <- length(score.dat.e[,1])
  contig.dat[i,7] <- paste('Archaeal: ',a.len,', Bacterial: ',b.len,', Eukaryotic: ',e.len,  sep = '')
  contig.dat[i,3] <- ifelse(mean(score.dat.a$a.e.val.num)==0, 100,  -log(mean(score.dat.a$a.e.val.num)))*(a.len)^(1/2)
  contig.dat[i,4] <- ifelse(mean(score.dat.b$b.e.val.num)==0, 100,  -log(mean(score.dat.b$b.e.val.num)))*(b.len)^(1/2)
  contig.dat[i,5] <- ifelse(mean(score.dat.e$b.e.val.num)==0, 100,  -log(mean(score.dat.e$e.e.val.num)))*(e.len)^(1/2)
  contig.dat[i,3] <- ifelse(is.na(contig.dat[i,3]), 0, contig.dat[i,3])
  contig.dat[i,4] <- ifelse(is.na(contig.dat[i,4]), 0, contig.dat[i,4])
  contig.dat[i,5] <- ifelse(is.na(contig.dat[i,5]), 0, contig.dat[i,5])
  if(contig.dat[i,3]>contig.dat[i,4] & contig.dat[i,3]>contig.dat[i,4]){
    contig.dat[i,2] <- 'Archaea'
  }
  if(contig.dat[i,4]>contig.dat[i,3] & contig.dat[i,4]>contig.dat[i,5]){
    contig.dat[i,2] <- 'Bacteria'
  }
  if(contig.dat[i,5]>contig.dat[i,3] & contig.dat[i,5]>contig.dat[i,4]){
    contig.dat[i,2] <- 'Eukaryotic'
  }
}

write.table(contig.dat, file ="output_scores.txt", sep="\t",quote = F, row.names = F, col.names = T )
write.table(contig.dat[contig.dat[,2]=='Archaea',], file = 'archaeal_contigs.txt', quote = F, row.names = F, col.names = T)

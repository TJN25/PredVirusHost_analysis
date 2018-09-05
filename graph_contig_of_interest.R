# produces a graph of the -log e-values of the proteins for the species of interest 
# this will generate a barchart for all the species in the file 'output_for_graph_species_of_interest.txt'
# be aware that on metagenomes this will generate a large number of graphs. It would be best to make a
# new file containing only the species you want. This can be done with 'rewrite_species_of_interest.py'
a <- read.table('output_for_graph_species_of_interest.txt', sep='\t')
b <- levels(a$V4)
library(lattice)
for (i in b){
  x <- subset(a, V4==i, select=c(V1,V2,V3))
  exten <- '.pdf' 
  pdf(paste(i,exten, sep=''), width=20, height=10)
  colours <- c('red','blue','yellow')
  y <- barchart(V2~V1,data=x,groups=V3,scales=list(x=list(rot=90,cex=0.8)), ylim=c(0,9.5), col=colours, xlab='Protein ID', ylab="Score", auto.key=list(space='right'),par.settings=list(superpose.polygon=list(col=colours)))
  print(y)
  dev.off() 
}

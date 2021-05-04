############################################################
##  Women First Trial: Guatemala Methylation Analysis (QTL Analysis)         
##  Written by Jessica Murphy 
##  Last edited on October 8, 2020.
##  This script counts the number of significant sites per chromosome
##  for several p-value thresholds based on both the nominal and
##  adjusted p-values. It produces a chromosome by threshold table
##  as well as the nominal and adjusted qqplots for each chromosome.
##  Please send any questions to jessica.murphy@ucdenver.edu
############################################################

setwd("/nfs/storage/math/gross-s2/projects/guatemala/Guatemala-omics/Methylation/QTL/results/")

# https://gettinggeneticsdone.blogspot.com/2010/07/qq-plots-of-p-values-in-r-using-base.html
# define the qqplot function
ggd.qqplot = function(pvector, main=NULL, ...) {
  o = -log10(sort(pvector,decreasing=F))
  e = -log10(1:length(o)/length(o) )
  plot(e,o,pch=19,cex=1, main=main, ...,
       xlab=expression(Expected~~-log[10](italic(p))),
       ylab=expression(Observed~~-log[10](italic(p))),
       xlim=c(0,max(e)), ylim=c(0,max(o)))
  lines(e,e,col="red")
}

# create an empty data frame to store the counts
out = data.frame(nominal5=numeric(), nominal6=numeric(), nominal7=numeric(), nominal8=numeric(), 
                 adjusted.05=numeric(), adjusted.01=numeric(), adjusted.001=numeric())

# filter throught the chromosomes
for (chrom in 1:22){
  
  # initiate a pdf to save the qqplots
  pdf(paste0("Chr", chrom, "_plots.pdf"))
  
  # read in the permutation results
  data = paste0("Chr", chrom, "_permutations.txt")
  d = read.table(data, head=FALSE, stringsAsFactors=FALSE)
  
  # count the number of significant sites per threshold
  out[chrom,] = c(length(which(d$V16<10e-5)), length(which(d$V16<10e-6)), length(which(d$V16<10e-7)),
              length(which(d$V16<10e-8)), length(which(d$V19<0.05)), length(which(d$V19<0.01)),
              length(which(d$V19<0.001)))
  
  # qqplot of nominal permutation p-values
  title16 = paste0("Chr", chrom, " Nominal")
  ggd.qqplot(d$V16, main=title16)
  
  # qqplot of adjusted permutation p-values
  title19 = paste0("Chr", chrom, " Adjusted")
  ggd.qqplot(d$V19, main=title19)
  
  # histogram of adjusted p-values
  hist(d$V19, main=title19, xlab="P-value")
  
  # turn plotting device off
  dev.off()
}

# save the table of significant sites per threshold
write.csv(out, "QTL_results.csv", quote=FALSE, row.names=FALSE, col.names=FALSE)

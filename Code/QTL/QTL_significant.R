############################################################
##  Women First Trial: Guatemala Methylation Analysis (QTL Analysis)         
##  Written by Jessica Murphy 
##  Last edited on October 17, 2020.
##  This script filters the significant CpG sites (adjusted p-value < 0.01)
##  from the permutation results.
##  Please send any questions to jessica.murphy@ucdenver.edu
############################################################

setwd("/nfs/storage/math/gross-s2/projects/guatemala/Guatemala-omics/Methylation/QTL/results/")

# loop through the chromosomes
for (chrom in 1:22){
  
  # read in the permutation results
  # would not recommend using data as a variable name (it's a reserved word in R)
  data = paste0("Chr", chrom, "_permutations.txt")
  d = read.table(data, header=F, stringsAsFactors=F)
  d = d[!is.na(d$V19),]
  
  # save the sites with an adjusted p-value < 0.01)
  name = paste0("Chr", chrom, "_significant.txt")
  write.table(d[which(d$V19<0.01), ], name, quote=F, row.names=F, col.names=F)
}

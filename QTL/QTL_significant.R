
setwd("/nfs/storage/math/gross-s2/projects/guatemala/Guatemala-omics/Methylation/QTL/results/")

for (chrom in 1:22){
  
  data = paste0("Chr", chrom, "_permutations.txt")
  d = read.table(data, header=F, stringsAsFactors=F)
  d=d[!is.na(d$V19),]
  
  name = paste0("Chr", chrom, "_significant.txt")
  
  write.table(d[which(d$V19<0.01), ], name, quote=F, row.names=F, col.names=F)
}
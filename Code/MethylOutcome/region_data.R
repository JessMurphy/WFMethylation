library(dplyr)

setwd("/nfs/storage/math/gross-s2/projects/guatemala/Guatemala-omics/Methylation/EWAS/results/")

MEs = read.table("ME_results.txt", header=T, sep='\t') 
top.int.fdr = read.table("top_int_fdr.txt", header=T, sep='\t') 
top.arm2.fdr = read.table("top_arm2_fdr.txt", header=T, sep='\t') 
top.BMI2.fdr = read.table("top_BMI2_fdr.txt", header=T, sep='\t')  

top.MEs = MEs %>% filter(int.fdr < 0.1 | arm.fdr < 0.1 | BMI.fdr < 0.1 | arm.fdr2 < 0.1 | BMI.fdr2 < 0.1) %>% select(methyl, Chr, pos)
top.int = top.int.fdr %>% filter(int.fdr < 0.05 | arm.fdr < 0.05 | BMI.fdr < 0.05) %>% select(methyl, Chr, pos)
top.arm2 = top.arm2.fdr %>% arrange(arm.p2) %>% filter(arm.fdr2 < 0.05) %>% top_n(-10, arm.fdr2) %>% select(methyl, Chr, pos)
top.BMI2 = top.BMI2.fdr %>% arrange(BMI.p2) %>% filter(BMI.fdr2 < 0.05) %>% select(methyl, Chr, pos)

#top = rbind(top.arm2, top.BMI2, top.int)

region_data <- function(top, name) {
  
  top$pos = as.numeric(top$pos)
  top$num = as.numeric(sapply(strsplit(top$Chr, "r", fixed=T), tail, 1))
  top$start = top$pos-1000
  top$end = top$pos+1000
  
  out = data.frame(rep(1, 105))
  results.keep = c()
  
  # loop through chromosomes
  for (i in 1:nrow(top)){
    
    # load the chromosome results 
    file.name2 = paste0("Chr", top$num[i], "_filtered.txt")
    
    # only chromosome 1 has a header
    if (top$num[i] != 1){
      
      results = read.table(file.name2, sep='\t')
      colnames(results) = c('methyl', 'BMI.p', 'arm.p', 'int.p', 'BMI.p2', 'arm.p2', 'arm2.beta1', 'arm3.beta1',
                            'BMI.beta1', 'arm2xbmi.beta1', 'arm3xbmi.beta1', 'arm2.beta2', 'arm3.beta2', 'BMI.beta2',
                            "NAs", "zeros", "non.zeros", "ones", "non.ones", "outliers", "diff.miss", "ME",
                            "Chr", "pos", "genes", "genes_flanking", "genes_type")
    } else {
      
      results = read.table(file.name2, header=T, sep='\t')
    }
    
    keep = between(results$pos, top$start[i], top$end[i])
    sites = results$methyl[keep]
    bin = rep(top$methyl[i], length(sites))
    
    results2 = results %>% filter(methyl %in% sites)
    results2$bin = bin
    
    results.keep = rbind(results.keep, results2)
    
    print(i)
  }
  
  write.table(results.keep, paste0(name, ".txt"), row.names=F, quote=F, sep='\t')
}

# call function
region_data(top.MEs, "top_ME_regions")
region_data(top.int, "top_int_regions")
region_data(top.arm2, "top_arm2_regions")
region_data(top.BMI2, "top_BMI2_regions")

library(dplyr)

setwd("/nfs/storage/math/gross-s2/projects/guatemala/Guatemala-omics/Methylation/BirthOutcomes/results/")

#top = c("chr2.229974407", "chr10.121593892", "chr13.75326388")

top.LGAZ = read.table("top_LGAZ_1000.txt", header=T, sep='\t') %>% filter(methyl.fdr < 0.05) %>% select(methyl, Chr, pos)
top.WGAZ = read.table("top_WGAZ_1000.txt", header=T, sep='\t') %>% filter(methyl.fdr < 0.05) %>% select(methyl, Chr, pos)
top.ME = read.table("ME_results_LGAZ.txt", header=T, sep='\t') %>% filter(methyl.fdr < 0.1) %>% select(methyl, Chr, pos)

region_data <- function(top, outcome, name) {
  
  top$pos = as.numeric(top$pos)
  top$num = as.numeric(sapply(strsplit(top$Chr, "r", fixed=T), tail, 1))
  top$start = top$pos-1000
  top$end = top$pos+1000
  
  out = data.frame(rep(1, 105))
  results.keep = c()
  
  # loop through chromosomes
  for (i in 1:nrow(top)){
    
    # load the chromosome results 
    file.name2 = paste0("Chr", top$num[i], "_filtered_", outcome, ".txt")
    
    # only chromosome 1 has a header
    if (top$num[i] != 1){
      
      results = read.table(file.name2, sep='\t')
      colnames(results) = c('methyl', 'methyl.p', 'arm2.p', 'arm3.p', 'BMI.p', 'age.p', 'expose.p', 'PC1.p',
                            'methyl.beta', 'arm2.beta', 'arm3.beta', 'BMI.beta', 'age.beta', 'expose.beta', 'PC1.beta',
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
region_data(top.LGAZ, "LGAZ", "top_LGAZ_regions")
region_data(top.WGAZ, "WGAZ", "top_WGAZ_regions")
region_data(top.ME, "LGAZ", "top_ME_regions")


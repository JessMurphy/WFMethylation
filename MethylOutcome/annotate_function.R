############################################################
##  Women First Trial: Guatemala Methylation Analysis         
##  Written by Jessica Murphy 
##  Last edited on January 7, 2021
##  This function annotates each chromosome with overlapping & 
##  flanking genes. The gene file was created using
##  /nfs/storage/math/gross-s2/projects/GencodeV31/Extract_genes.R.
##  Please send any questions to jessica.murphy@ucdenver.edu
############################################################

setwd('/nfs/storage/math/gross-s2/projects/Gencode/')
genes <- read.table('Gencode.v31.annotation.genes.txt', header =  FALSE, sep = '\t')

setwd('/nfs/storage/math/gross-s2/projects/guatemala/Guatemala-omics/Methylation/EWAS/results/')

for(i in 1:22){

   # read in chromosome results 
    file.name = paste0("Chr", i, "_results", ".txt")
    re = read.table(file.name, header = T, sep = '\t')
    
    # split the Chr and pos
    re$Chr = '.' 
    re$pos =  '.'
    re$methyl = as.character(re$methyl)
    
    for(j in 1:nrow(re)){
      temp = strsplit(re$methyl[j], split = "[.]")
      re$Chr[j] = temp[[1]][1]
      re$pos[j] = temp[[1]][2]
    }
    
    # make sure position is numeric
    re$pos = as.numeric(re$pos)
    
    # subset the genes for the specified chromosome
    temp_g = genes[which(genes$V1 == paste0("chr", i)),]
    
    # annotate results with overlapping/flanking genes
    re$genes = '.'
    re$genes_flanking =  '.'
    re$genes_type = '.'
    n = 5000 # flanking region around the CpG
    
    for(k in 1:nrow(re)){
      
      # genes overalapping CpG
      gen = temp_g[which(temp_g$V4 <= re$pos[k] & temp_g$V5 >= re$pos[k]),]
      
      # if one gene, fill in the row
      if(nrow(gen)==1){ 
        re$genes[k] = as.character(gen$V10) 
      }
      
      # if mutliple genes, separate with a comma
      if(nrow(gen)>1){ 
        over = gen$V10
        over1 = paste0(over, collapse=",")
        re$genes[k] = as.character(over1)
      }
      
      # genes within the CpG flanking region
      gen_flank = temp_g[which(temp_g$V4 <= (re$pos[k] + n) & temp_g$V5 >= (re$pos[k] - n)),] 
      
      # if one gene, fill in the row
      if(nrow(gen_flank)==1){
        re$genes_flanking[k] = as.character(gen_flank$V10) 
        re$genes_type[k] = as.character(gen_flank$V2) 
      }
      
      # if multiple genes, separate with a comma
      if(nrow(gen_flank)>1){ 
        flank = gen_flank$V10 
        flank1 = paste0(flank, collapse=",")
        re$genes_flanking[k] = as.character(flank1)
        
        type = gen_flank$V2
        type1 = paste0(type, collapse=",")
        re$genes_type[k] = as.character(type1)
      }
      
      #if(round(k/1000)==(k/1000)) {print(k)}
    }
    
    # save annotated results
    write.table(re, paste0('Chr', i, '_annotated', '.txt'), row.names=F, quote=F, sep = '\t')
    
    print(i)
}

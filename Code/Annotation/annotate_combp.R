############################################################
##  Women First Trial: Guatemala Methylation Analysis       
##  Written by Jessica Murphy 
##  Last edited on January 7, 2021
##  This function annotates each combp region with overlapping & 
##  flanking genes. The gene file was created using Extract_genes.R.
##  Please send any questions to jessica.murphy@ucdenver.edu
############################################################

# read in the gene information
setwd("~/Methylation/code/annotation/")
genes = read.table('Gencode.v31.annotation.genes.txt', header =  FALSE, sep = '\t')

# read in the methylation outcome combp results
setwd("~/Methylation/results/combp/")
arm = read.table("combp_arm_results.regions-p.txt", header=T)
BMI = read.table("combp_BMI_results.regions-p.txt", header=T)
int = read.table("combp_interaction_results.regions-p.txt", header=T)

# read in the methylation outcome ME combp results
arm.ME = read.table("ME_combp_arm_results.regions-p.txt", header=T)
int.ME = read.table("ME_combp_interaction_results.regions-p.txt", header=T)

# read in the birth outcomes combp results
LGAZ = read.table("combp_LGAZ_results.regions-p.txt", header=T)
WGAZ = read.table("combp_WGAZ_results.regions-p.txt", header=T)
HCGAZ = read.table("combp_HCGAZ_results.regions-p.txt", header=T)
WLGAZ = read.table("combp_WLGAZ_results.regions-p.txt", header=T)

# define annotation function
# inputs: results - dataframe of combp results
#         name - character string to name the output file
# output: tab-separated text file of annotated combp results
annotate_combp <- function(results, name) {
  
  # annotate results with overlapping/flanking genes
  results$genes = '.'
  results$genes_flanking =  '.'
  results$genes_type = '.'
  n = 5000 # flanking region around the CpG
  
  for(i in 1:nrow(results)){
    
    # subset the genes for the specified chromosome
    temp_g = genes[which(genes$V1 == results$chrom[i]),]
    
    # genes overalapping CpG (V4 - start of gene, V5 - end of gene)
    gen = temp_g %>% filter(V4 <= results$end[i] & V4 >= results$start[i] | 
                              V5 <= results$end[i] & V5 >= results$start[i] |
                              V4 <= results$start[i] & V5 >= results$end[i])
    # either the start of the gene is within the region, 
    # the end of the gene is within the region, or the
    # gene encompasses the entire region
    
    # if one gene, fill in the row
    if(nrow(gen)==1){ 
      results$genes[i] = as.character(gen$V10) 
    }
    
    # if mutliple genes, separate with a comma
    if(nrow(gen)>1){ 
      over = gen$V10
      over1 = paste0(over, collapse=",")
      results$genes[i] = as.character(over1)
    }
    
    # genes within the CpG flanking region
    gen_flank = temp_g %>% filter(V4 <= (results$end[i] + n) & V4 >= (results$start[i] - n) | 
                                    V5 <= (results$end[i] + n) & V5 >= (results$start[i] - n) |
                                    V4 <= (results$start[i] - n) & V5 >= (results$end[i] + n))
    
    # if one gene, fill in the row
    if(nrow(gen_flank)==1){
      results$genes_flanking[i] = as.character(gen_flank$V10) 
      results$genes_type[i] = as.character(gen_flank$V2) 
    }
    
    # if multiple genes, separate with a comma
    if(nrow(gen_flank)>1){ 
      flank = gen_flank$V10 
      flank1 = paste0(flank, collapse=",")
      results$genes_flanking[i] = as.character(flank1)
      
      type = gen_flank$V2
      type1 = paste0(type, collapse=",")
      results$genes_type[i] = as.character(type1)
    }
    
    print(i)
  }
  
  # save annotated results
  write.table(results, name, row.names=F, quote=F, sep = '\t')
}

# call function for methylation outcome results
annotate_combp(int, "combp_int_results_annotated.txt")
annotate_combp(arm, "combp_arm_results_annotated.txt")
annotate_combp(BMI, "combp_BMI_results_annotated.txt")

# call function for ME methylation outcome results
annotate_combp(int.ME, "ME_combp_int_results_annotated.txt")
annotate_combp(arm.ME, "ME_combp_arm_results_annotated.txt")

# call function for birth outcomes results
annotate_combp(LGAZ, "combp_LGAZ_results_annotated.txt")
annotate_combp(WGAZ, "combp_WGAZ_results_annotated.txt")
annotate_combp(HCGAZ, "combp_HCGAZ_results_annotated.txt")
annotate_combp(WLGAZ, "combp_WLGAZ_results_annotated.txt")

############################################################
##  Women First Trial: Guatemala Methylation Analysis         
##  Written by Jessica Murphy 
##  Last edited on December 15, 2020
##  This function performs an EWAS for each chromosome based on a 
##  specified phenotype (LGAZ, WGAZ, HCGAZ, or WLGAZ). The main model is 
##  phenotype = % methylation + arm + BMI + covariates. The results 
##  provide the p-values and beta estimates for the CpG site, arm, 
##  and BMI. 
##  Please send any questions to jessica.murphy@ucdenver.edu
############################################################

library(dplyr)

# define a function to remove outliers
# https://stackoverflow.com/questions/12866189/calculating-the-outliers-in-r
FindOutliers <- function(data) {
  
  lowerq = quantile(data, na.rm=T)[2]
  upperq = quantile(data, na.rm=T)[4]
  iqr = upperq - lowerq
  
  # identify extreme outliers
  threshold.upper = (iqr * 3) + upperq
  threshold.lower = lowerq - (iqr * 3)
  result = which(data > threshold.upper | data < threshold.lower)
}

# define function to perform EWAS
# data: methylation matrix (samples by sites) for a single chromosome with additional variables in the first 29 columns
# outcome: age-adjusted phenotype ("LGAZ", "WGAZ", "HCGAZ", or "WLGAZ")
# chrom: chromosome number for data
# writes a text file of results (coefficients and p-values)

pheno <- function(data, outcome, chrom){
  
  # read in ME data & store ME site names
  ME.data = read.table("/nfs/storage/math/gross-s2/projects/guatemala/Guatemala-omics/Methylation/EWAS/data/ME_all_data4analysis.txt", header=T)
  MEs = names(ME.data)[15:ncol(ME.data)]
  
  # remove the preterm babies
  preterm = c("631530U", "663040J", "606150A", "616780H", "683170H", "652670K")
  data2 = data %>% filter(!wf_id %in% preterm) %>% filter(!is.na(LGAZ))
  
  out = c()
  
  for(i in 30:length(data2)){ #loop through the CpG sites
    
    # count the number of initial NA values
    NAs = sum(is.na(data2[,i]))
    
    # determine the outliers & set to NA
    outliers = FindOutliers(data2[,i])
    data2[outliers, i] = NA
    
    # count the number of outliers, nonzero reads, and zero reads
    n.outliers = length(outliers)
    zeros = length(which(data2[,i]==0))
    non_zeros = length(which(data2[,i]>0)) 

    # count the number of ones and non-ones
    ones = length(which(data2[,i]==1)) 
    non_ones = length(which(data2[,i]<1))

    # test for differential missingness
    temp = data2[,c("arm", names(data2)[i])]
    temp[,2] = is.na(temp[,2])
    tbl = table(temp$arm, temp[,2])
    
    if(dim(tbl)[2] > 1){
      test = fisher.test(tbl)
      diff_miss = test$p.value
    }
    if(dim(tbl)[2] == 1){
      diff_miss = 1
    }
    
    # define models
    covariates = "+ age + exposure + PC1 +"
    func = paste(outcome, " ~ as.factor(arm) + BMI_calc2", covariates, names(data)[i], sep="")
    
    # fit models
    lmod = lm(func, data=data2)
    
    if (is.na(lmod$coefficients[8])){
      methyl.beta = NA
      methyl.p = NA
    } else {
      methyl.beta = summary(lmod)$coefficients[8,1]
      methyl.p = summary(lmod)$coefficients[8,4]
    }
    
    # beta coefficients
    arm2.beta = summary(lmod)$coefficients[2,1]
    arm3.beta = summary(lmod)$coefficients[3,1]
    BMI.beta = summary(lmod)$coefficients[4,1]
    age.beta = summary(lmod)$coefficients[5,1]
    expose.beta = summary(lmod)$coefficients[6,1]
    PC1.beta = summary(lmod)$coefficients[7,1]
    
    # p-values 
    arm2.p = summary(lmod)$coefficients[2,4]
    arm3.p = summary(lmod)$coefficients[3,4]
    BMI.p = summary(lmod)$coefficients[4,4]
    age.p = summary(lmod)$coefficients[5,4]
    expose.p = summary(lmod)$coefficients[6,4]
    PC1.p = summary(lmod)$coefficients[7,4]
    
    # determine if site is an ME
    CpG = names(data2)[i]
    ME = ifelse(CpG %in% MEs, "ME", ".")
    
    re = c(names(data)[i], methyl.p, arm2.p, arm3.p, BMI.p, age.p, expose.p, PC1.p, 
           methyl.beta, arm2.beta, arm3.beta, BMI.beta, age.beta, expose.beta, PC1.beta,
           NAs, zeros, non_zeros, ones, non_ones, n.outliers, diff_miss, ME) 
    
    out = rbind(out, re)
    
    if(round(i/1000)==(i/1000)) {print(i)}
  }
  
  colnames(out) = c('methyl', 'methyl.p', 'arm2.p', 'arm3.p', 'BMI.p', 'age.p', 'expose.p', 'PC1.p', 
                    'methyl.beta', 'arm2.beta', 'arm3.beta', 'BMI.beta', 'age.beta', 'expose.beta', 'PC1.beta',
                    "NAs", "zeros", "non.zeros", "ones", "non.ones", "outliers", "diff.miss", "ME")
  
  write.table(out, paste0('Chr', chrom, '_results_', outcome, '.txt'), row.names=F, quote=F, sep='\t')

  out = as.data.frame(out)

  print(paste0(length(out$diff_miss < 0.001), " sites with diff.miss p-value < 0.001"))
}

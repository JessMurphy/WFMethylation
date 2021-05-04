############################################################
##  Women First Trial: Guatemala Methylation Analysis (Methylation Outcome)         
##  Written by Jessica Murphy 
##  Last edited on December 10, 2020
##  This function performs a methylation analysis for each chromosome. 
##  The main model is % methylation = arm + BMI + arm*BMI + covariates 
##  and another model without the interaction term is also fit. The  
##  results provide the p-values and beta estimates for arm, BMI, and 
##  the interaction between the two for the main model and the p-values 
##  and beta estimates for just arm and BMI for the secondary model 
##  without the interaction term. 
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
# chrom: chromosome number for data
# writes a text file of results (coefficients and p-values with and without an interaction term)

analysis <- function(data, chrom){
  
  # read in ME data & store ME site names
  ME.data = read.table("/nfs/storage/math/gross-s2/projects/guatemala/Guatemala-omics/Methylation/EWAS/data/ME_all_data4analysis.txt", header=T)
  MEs = names(ME.data)[15:ncol(ME.data)]
  
  # remove the preterm babies
  preterm = c("631530U", "663040J", "606150A", "616780H", "683170H", "652670K")
  data2 = data %>% filter(!wf_id %in% preterm)
  
  # filter data by arm
  arm1 = data2 %>% filter(arm==1)
  arm2 = data2 %>% filter(arm==2)
  arm3 = data2 %>% filter(arm==3)
  
  # empty vector to store results
  out = c()
  
  for(i in 30:length(data2)){ #loop through the CpG sites
    
    # count the number of initial NA values
    NAs = sum(is.na(data2[,i]))
    
    # determine the outliers for arm 1 & set to NA
    out1 = FindOutliers(arm1[,i])
    out_ind1 = which(data2$wf_id %in% arm1$wf_id[out1])
    data2[out_ind1, i] = NA
    
    # determine the outliers for arm 2 & set to NA
    out2 = FindOutliers(arm2[,i])
    out_ind2 = which(data2$wf_id %in% arm2$wf_id[out2])
    data2[out_ind2, i] = NA
    
    # determine the outliers for arm 3 & set to NA
    out3 = FindOutliers(arm3[,i])
    out_ind3 = which(data2$wf_id %in% arm3$wf_id[out3])
    data2[out_ind3, i] = NA
    
    # count the number of outliers, nonzero reads, and zero reads
    outliers = length(out_ind1) + length(out_ind2) + length(out_ind3)
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
 
    # need the main model with interaction (func1)
    # if interaction isn't significant, remove it (func2)
    
    # define models
    covariates = "+ age + sex + exposure + PC1 + refactor1 + refactor2"
    func1 = paste(names(data)[i], "~as.factor(arm)*BMI_calc2", covariates, sep="")
    func2 = paste(names(data)[i], "~as.factor(arm) + BMI_calc2", covariates, sep="")
    func.null.BMI = paste(names(data)[i], "~as.factor(arm)", covariates, sep="")
    func.null.arm = paste(names(data)[i], "~BMI_calc2", covariates, sep="")
    
    # fit models
    temp1 = glm(func1, data=data2)
    temp2 = glm(func2, data=data2)
    temp.null.BMI = glm(func.null.BMI, data=data2)
    temp.null.arm = glm(func.null.arm, data=data2)
    
    # p-values with interaction
    BMI.p = anova(temp.null.BMI, temp1, test="LRT")[2,5]
    arm.p = anova(temp.null.arm, temp1, test="LRT")[2,5]
    int.p = anova(temp2, temp1, test="LRT")[2,5]
    
    # p-values without interaction
    BMI.p2 = anova(temp.null.BMI, temp2, test="LRT")[2,5]
    arm.p2 = anova(temp.null.arm, temp2, test="LRT")[2,5]
    
    # betas with interaction
    arm2.beta1 = summary(temp1)$coefficients[2,1]
    arm3.beta1 = summary(temp1)$coefficients[3,1]
    BMI.beta1 = summary(temp1)$coefficients[4,1]
    arm2xbmi.beta1 = summary(temp1)$coefficients[11,1]
    arm3xbmi.beta1 = summary(temp1)$coefficients[12,1]
    
    # betas without interaction
    arm2.beta2 = summary(temp2)$coefficients[2,1]
    arm3.beta2 = summary(temp2)$coefficients[3,1]
    BMI.beta2 = summary(temp2)$coefficients[4,1]
    
    # determine if site is an ME
    CpG = names(data2)[i]
    ME = ifelse(CpG %in% MEs, "ME", ".")
    
    re = c(CpG, BMI.p, arm.p, int.p, BMI.p2, arm.p2, arm2.beta1, arm3.beta1, BMI.beta1, arm2xbmi.beta1, arm3xbmi.beta1,
             arm2.beta2, arm3.beta2, BMI.beta2, NAs, zeros, non_zeros, ones, non_ones, outliers, diff_miss, ME)
    
    out = rbind(out, re)
    
    if(round(i/1000)==(i/1000)) {print(i)}
  }
  
  colnames(out) = c('methyl', 'BMI.p', 'arm.p', 'int.p', 'BMI.p2', 'arm.p2', 'arm2.beta1', 'arm3.beta1',
                    'BMI.beta1', 'arm2xbmi.beta1', 'arm3xbmi.beta1', 'arm2.beta2', 'arm3.beta2', 'BMI.beta2',
                    "NAs", "zeros", "non.zeros", "ones", "non.ones", "outliers", "diff.miss", "ME")
  
  write.table(out, paste0('Chr', chrom, '_results.txt'), row.names=F, quote=F, sep='\t')

  out = as.data.frame(out)

  print(paste0(length(out$diff_miss < 0.001), " sites with diff.miss p-value < 0.001"))
}
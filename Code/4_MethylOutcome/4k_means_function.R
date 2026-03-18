############################################################
##  Women First Trial: Guatemala Methylation Analysis (Methylation Outcome)      
##  Written by Jessica Murphy 
##  Last edited on May 3, 2021.
##  This script adds the arm contrasts as well as the trends
##  (interaction model) or marginal means (non-interaction 
##  model) to the top results.
##  Please send any questions to jessica.murphy@ucdenver.edu
############################################################

# load libraries
library(dplyr)
library(scales) #scientific
library(emmeans)

setwd("~/RESEARCH/Methylation/results/EWAS/")

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

# read in top results
#ME.results = read.table("ME_results.txt", header=T, sep='\t') #%>% filter(int.fdr<0.1) 
#top.int.fdr = read.table("top_int_fdr.txt", header=T, sep='\t') 
#top.arm2.fdr = read.table("top_arm2_fdr.txt", header=T, sep='\t') %>% filter(arm.fdr2 <= 0.05)
#top.BMI2.fdr = read.table("top_BMI2_fdr.txt", header=T, sep='\t') %>% filter(BMI.fdr2 <= 0.05)

# read in top data
#ME.data = read.table("ME_data.txt", header=T, sep='\t') #%>% select(infant_id:refactor2, ME.results$methyl)
#int.data = read.table("top_int_data.txt", header=T, sep='\t') 
#arm2.data = read.table("top_arm2_data.txt", header=T, sep='\t') 
#BMI2.data = read.table("top_BMI2_data.txt", header=T, sep='\t') 

# define merge_means function
# inputs: results - a dataframe of top results 
#         data - a dataframe of methylation data for the top sites
#         name - character string of the output file name 
#         int - interaction term? ("yes" or "no")
# output: a tab-separated text file of the updated top results with the emmeans/emtrends info

merge_means <- function(results, data, name, int) {
  
  # remove the preterm babies
  preterm = c("631530U", "663040J", "606150A", "616780H", "683170H", "652670K")
  data2 = data %>% filter(!wf_id %in% preterm)
  
  # filter data by arm
  arm1 = data2 %>% filter(arm==1)
  arm2 = data2 %>% filter(arm==2)
  arm3 = data2 %>% filter(arm==3)
  
  for(i in 30:length(data2)){ #loop through the CpG sites
    
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
    
  } # end CpG site loop
  
  # convert categorical variables to factors
  data2$arm = as.factor(data2$arm)
  data2$sex = as.factor(data2$sex)
  data2$exposure = as.factor(data2$exposure)
  
  # create empty vector to store results in
  out = c()
  
  # loop through the updated CpG sites
  for (k in 30:length(data2)){
    
    # define the covariates for the model
    covariates = "+ age + sex + exposure + PC1 + refactor1 + refactor2"
    
    if (int=="yes") { # interaction model
      
      # fit linear model
      func = paste0(colnames(data2)[k], "~arm*BMI_calc2", covariates)
      lmod = lm(func, data2)
      
      # estimate trends
      em = emtrends(lmod, pairwise ~ arm, var = 'BMI_calc2')
      means = data.frame(em$emtrends)
      
      # define the output column names
      names = c("methyl", "arm1.trend", "arm1.CI", "arm2.trend", "arm2.CI", "arm3.trend", "arm3.CI",
                "arm12.est", "arm12.CI", "arm12.p", "arm13.est", "arm13.CI", "arm13.p",
                "arm23.est", "arm23.CI", "arm23.p")
      
    } else if (int=="no") { # non-interaction model
      
      # fit linear model
      func = paste0(colnames(data2)[k], "~arm+BMI_calc2", covariates)
      lmod = lm(func, data2)
      
      # estimate marginal means
      em = emmeans(lmod, pairwise ~ arm)
      means = data.frame(em$emmeans)
      
      # define the output column names
      names = c("methyl", "arm1.mean", "arm1.CI", "arm2.mean", "arm2.CI", "arm3.mean", "arm3.CI",
                "arm12.est", "arm12.CI", "arm12.p", "arm13.est", "arm13.CI", "arm13.p",
                "arm23.est", "arm23.CI", "arm23.p")
    }
    
    # calculate the confidence intervals for the contrasts
    contrasts = data.frame(em$contrasts)
    contrasts$lower = apply(contrasts, 1, function(x){round(as.numeric(x[2])-1.96*as.numeric(x[3]), 4)})
    contrasts$upper = apply(contrasts, 1, function(x){round(as.numeric(x[2])+1.96*as.numeric(x[3]), 4)})
    contrasts$CI = paste0("[", contrasts$lower, "  ", contrasts$upper, "]")
    
    # format the confidence intervals for the means/trends
    means$CI = paste0("[", round(means$lower.CL, 4), "  ", round(means$upper.CL, 4), "]")
    
    # combine the contrasts and meanss/trends info for each site
    temp = c(colnames(data2)[k], means[1,2], means[1,7], means[2,2], means[2,7], means[3,2], means[3,7],
              contrasts[1,2], contrasts[1,9], contrasts[1,6], contrasts[2,2], contrasts[2,9], contrasts[2,6],
              contrasts[3,2], contrasts[3,9], contrasts[3,6])
    
    # name the output columns
    names(temp) = names
    
    # store the results per site
    out = rbind(out, temp)
    
  } # end data2 loop
  
  # merge the previous results with the emmeans/emtrends output
  out = as.data.frame(out)
  results2 = merge(results, out, by="methyl")
  
  # save the updated results
  write.table(results2, paste0(name, ".txt"), row.names = FALSE, quote = FALSE, sep='\t')
  
  return(results2)
}

# call function
#out.int = merge_means(top.int.fdr, int.data, "top_int_updated", "yes")
#out.arm2 = merge_means(top.arm2.fdr, arm2.data, "top_arm2_updated", "no")
#out.BMI2 = merge_means(top.BMI2.fdr, BMI2.data, "top_BMI2_updated", "no")
#out.ME = merge_means(ME.results, ME.data, "ME_updated", "no")

############################################################
##  Women First Trial: Guatemala Methylation Analysis (Methylation Outcome)        
##  Written by Jessica Murphy 
##  Last edited on May 3, 2021.
##  This script creates interaction plots and methylation 
##  boxplots for the top interaction results. It also creates
##  marginal means plots and methylation boxplots for the 
##  top arm results for the non-interaction model.
##  Please send any questions to jessica.murphy@ucdenver.edu
############################################################

# load libraries
library(sjPlot)
library(sjmisc)
library(ggplot2)
library(emmeans)
library(tidyr) #gather
library(dplyr)

set_theme(base = theme_bw())

#setwd("~/RESEARCH/Methylation/results/EWAS/")

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

# read in the data for the top results
data = read.table("top_int_data.txt", header=T, sep='\t') 
#data = read.table("top_arm2_data.txt", header=T, sep='\t') 

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
}

# convert categorical variables to factors
data2$sex = as.factor(data2$sex)
data2$exposure = as.factor(data2$exposure)
data2$arm = as.factor(data2$arm)

# define the top interaction sites to plot
top.int = c("chr2.236049031", "chr2.98139266", "chr4.72765908", "chr5.55552531", 
            "chr6.157299979", "chr7.43550516", "chr10.132786523", "chr12.1527984", 
            "chr15.101046848", "chr15.85758729", "chr16.88067925", "chr16.88470918", 
            "chr16.88862298", "chr17.6894527", "chr18.78126680", "chr20.63288146", "chr21.44350139")

# define the top arm sites to plot (non-interaction model)
#top.arm = c("chr1.144552011", "chr1.154405049", "chr1.26307043", "chr2.129980567", "chr2.231715445",
#            "chr3.100400230", "chr3.10115341", "chr18.22176690", "chr19.17205492", "chr19.49423303")

# INTERACTION PLOTS

# create an empty list to store the interaction plots in
int.plots = list()

# loop through the top sites
for (i in 1:length(top.int)){
  
  # fit a linear model to the data (with interaction)
  covariates = "+ age + sex + exposure + PC1 + refactor1 + refactor2"
  func = paste0(top.int[i], "~arm*BMI_calc2", covariates)
  lmod = lm(func, data2)
  
  # create an interaction plot
  plot = sjPlot::plot_model(lmod, type="pred", show.data=T, terms=c("BMI_calc2", "arm"))
  
  # store the plot per site
  int.plots[[i]] = plot
}

# name the plots the top interaction sites
names(int.plots) = top.int

# NON-INTERACTION PLOTS

# create an empty list to store the plots in
#arm.plots = list()

# loop through the top interaction sites
#for (i in 1:length(top.arm)){
  
  # fit a linear model to the data (without interaction)
  #covariates = "+ age + sex + exposure + PC1 + refactor1 + refactor2"
  #func = paste0(top.arm[i], "~arm+BMI_calc2", covariates)
  #lmod = lm(func, data2)
  
  # create a marginal means plot
  #plot = plot(emmeans(lmod, pairwise ~ arm))
  
  # store the plot per site
  #arm.plots[[i]] = plot
#}

# name the plots the top arm sites
#names(arm.plots) = top.arm

# METHYLATION BOXPLOTS

# susbet data to just top results
plot.data = data2 %>% select(all_of(top.int))
#plot.data = data2 %>% select(all_of(top.arm))

# convert the data to long format
data.long = gather(plot.data, position, methylation, chr2.236049031:chr21.44350139)
#data.long = gather(plot.data, position, methylation, chr1.144552011:chr19.49423303)

# convert the chromosome positions into a factor variable
data.long$position = factor(data.long$position, levels=top.int)
#data.long$position = factor(data.long$position, levels=top.arm)

# create boxplots of the percent methylation for the top hits
ggplot(data.long, aes(position, methylation, col=position)) +
  geom_boxplot(size=0.8) + theme(axis.text.x=element_blank(), axis.title.x=element_blank())

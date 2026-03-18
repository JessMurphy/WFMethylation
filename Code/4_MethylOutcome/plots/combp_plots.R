############################################################
##  Women First Trial: Guatemala Methylation Analysis (Methylation Outcome)       
##  Written by Jessica Murphy 
##  Last edited on May 4, 2021.
##  This script creates region plots of the marginal means (arm),
##  beta estimates (BMI), or estimated trends (interaction) for
##  the top combp regions annotated to genes.
##  Please send any questions to jessica.murphy@ucdenver.edu
############################################################

# load libraries
library(dplyr)
library(tidyr)
library(ggplot2)

#setwd("~/RESEARCH/Methylation/results/combp/")
setwd("C:/Users/jimur/OneDrive/Documents/RESEARCH/Methylation/results/combp/")

# read in combp results (filter and add region)
arm.regions = read.table("combp_arm_results.txt", header=T, sep='\t') %>%
  mutate(region=paste0(chrom, ":", start, "-", end)) %>% filter(z_sidak_p<=0.05, n_probes>3)
bmi.regions = read.table("combp_bmi_results.txt", header=T, sep='\t') %>%
  mutate(region=paste0(chrom, ":", start, "-", end)) %>% filter(z_sidak_p<=0.05, n_probes>3)
int.regions = read.table("combp_int_results.txt", header=T, sep='\t') %>%
  mutate(region=paste0(chrom, ":", start, "-", end)) %>% filter(z_sidak_p<=0.05, n_probes>3)

# just look at the top five results per term annotated to a gene
top.arm.regions = arm.regions %>% filter(genes!='.') %>% arrange(z_sidak_p) %>% slice(1:5)
top.bmi.regions = bmi.regions %>% filter(genes!='.') %>% arrange(z_sidak_p) %>% slice(1:5)
top.int.regions = int.regions %>% filter(genes!='.') %>% arrange(z_sidak_p) %>% slice(1:5)

# read in the single site results corresponding to the top five regions for each term
arm.sites = read.table("combp_arm_sites.txt", header=T, sep='\t') %>% filter(region %in% top.arm.regions$region) %>% 
  select(-chrom, -end, -Chr) %>% mutate(term="arm")
bmi.sites = read.table("combp_bmi_sites.txt", header=T, sep='\t') %>% filter(region %in% top.bmi.regions$region) %>% 
  select(-chrom, -end, -Chr) %>% mutate(term="bmi")
int.sites = read.table("combp_int_sites.txt", header=T, sep='\t') %>% filter(region %in% top.int.regions$region) %>% 
  select(-chrom, -end, -Chr) %>% mutate(term="int")

# read in the methylation data corresponding to the top combp regions
arm.data = read.table("arm_combp_data.txt", header=T, sep='\t') %>% select(1:29, arm.sites$methyl)
bmi.data = read.table("bmi_combp_data.txt", header=T, sep='\t') %>% select(1:29, bmi.sites$methyl)
int.data = read.table("int_combp_data.txt", header=T, sep='\t') %>% select(1:29, int.sites$methyl)

# link to function to estimate emmeans/emtrends
source("C:/Repositories/WFMethylation/Code/MethylOutcome/means_function.R")

# call merge_means function
arm.sites2 = merge_means(arm.sites, arm.data, "combp_arm_sites_updated.txt", "no")
bmi.sites2 = merge_means(bmi.sites, bmi.data, "combp_bmi_sites_updated.txt", "no")
int.sites2 = merge_means(int.sites, int.data, "combp_int_sites_updated.txt", "yes")

# determine if any sites were excluded from the combp n_probes count
which(arm.sites2$region.q>0.05) #no
which(bmi.sites2$region.q>0.05) #yes
which(int.sites2$region.q>0.05) #yes

# colorblind friendly palette
cbPalette = c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

########## ARM PLOTS ##########

# format the marginal means and arm p-value as numeric
arm.sites2 = arm.sites2 %>% mutate_at(vars(contains("mean"), arm.p2), as.numeric)

# convert the arm data to long format based on the marginal means
arm.longer = arm.sites2 %>% 
  pivot_longer(cols=contains("mean"), names_to = "arm", values_to = "mean")

# extract just the arm number for the arm column
arm.longer$arm = sapply(strsplit(arm.longer$arm, ".", fixed=T), head, 1)
arm.longer$arm = as.factor(substr(arm.longer$arm, 4, 4))

# extract the combp regions for arm
arm.longer$region = as.factor(arm.longer$region)
arm.regions = levels(arm.longer$region)

# create an empty list to store the plots in
arm.plots = list()

# loop through the arm regions
for(j in arm.regions){
  
  # subset the arm data to the specified region
  arm.temp = arm.longer[which(arm.longer$region == j),]
  print(nrow(arm.temp))
  
  # calculate the median marginal mean per arm for the specified region
  arm.med = arm.temp %>% group_by(arm) %>% summarize(median=median(mean))
  
  # create a regional plot for the marginal means
  arm.plots[[j]] = ggplot() +
    geom_point(data=arm.temp, aes(x=pos, y=mean, color=arm, size=arm.p2), alpha=0.5) +
    geom_hline(aes(yintercept=0), show.legend=F) + theme_bw() +
    geom_hline(data=arm.med, aes(yintercept=median, col=arm)) +
    scale_size("p-value", trans="log10", range=c(14, 0.5), 
               breaks=c(1e-7, 1e-5, 1e-3, 1e-1)) +
    scale_color_manual(values=cbPalette[c(6,4,2)]) +
    labs(x="Position", y="Marginal Mean", title=j)
}

########## BMI PLOTS ##########

# extract the combp regions for BMI
bmi.sites2$region = as.factor(bmi.sites2$region)
bmi.regions = levels(bmi.sites2$region)

# add a variable to differentiate which sites were included in the combp n_probes
bmi.sites2$combp = ifelse(bmi.sites2$region.q>0.05, "no", "yes")
bmi.sites2$combp = factor(bmi.sites2$combp, levels=c("yes", "no"))

# create an empty list to store the plots in
bmi.plots = list()

# loop through the BMI regions
for(j in bmi.regions){
  
  # subset the BMI data to the specified region
  bmi.temp = bmi.sites2[which(bmi.sites2$region == j),]
  print(nrow(bmi.temp))
  
  # calculate the median beta estimate for the specified region
  bmi.med = data.frame(median=median(bmi.temp$BMI.beta2))
  
  if (summary(bmi.temp$combp)["no"]==0){ #no sites excluded
    
    # create a regional plot for the beta estimates
    bmi.plots[[j]] = ggplot() + theme_bw() +
      geom_point(data=bmi.temp, aes(x=pos, y=BMI.beta2, size=BMI.p2), alpha=0.5, col=cbPalette[7]) +
      geom_hline(data=bmi.med, aes(yintercept=median), col=cbPalette[7]) +
      scale_size("p-value", trans="log10", range=c(14, 0.5), 
                 breaks=c(1e-5, 1e-4, 1e-3, 1e-2, 1e-1)) +
      geom_hline(aes(yintercept=0), show.legend=F) +
      labs(x="Position", y="Beta Estimate", title=j)
  
    } else { #sites excluded from combp
    
    # create a regional plot for the beta estimates (with combp excluded sites in grey)
    bmi.plots[[j]] = ggplot() + theme_bw() +
      geom_point(data=bmi.temp, aes(x=pos, y=BMI.beta2, size=BMI.p2, col=combp), alpha=0.5) +
      geom_hline(data=bmi.med, aes(yintercept=median), col=cbPalette[7]) +
      scale_size("p-value", trans="log10", range=c(14, 0.5), 
                 breaks=c(1e-5, 1e-4, 1e-3, 1e-2, 1e-1)) +
      geom_hline(aes(yintercept=0), show.legend=F) + 
      scale_color_manual(values=cbPalette[c(7,1)]) +
      labs(x="Position", y="Beta Estimate", title=j)
  }
}

########## INT PLOTS ##########

# format the trends and interaction p-value as numeric
int.sites2 = int.sites2 %>% mutate_at(vars(contains("trend"), int.p), as.numeric)

# convert the arm data to long format based on the estimated trends
int.longer = int.sites2 %>% 
  pivot_longer(cols=contains("trend"), names_to = "arm", values_to = "trend")

# extract just the arm number for the arm column
int.longer$arm = sapply(strsplit(int.longer$arm, ".", fixed=T), head, 1)
int.longer$arm = as.factor(substr(int.longer$arm, 4, 4))

# extract the combp regions for the interaction term
int.longer$region = as.factor(int.longer$region)
int.regions = levels(int.longer$region)

# add a variable to differentiate which sites were included in the combp n_probes
int.longer$combp = ifelse(int.longer$region.q>0.05, "no", "yes")
int.longer$combp = factor(int.longer$combp, levels=c("yes", "no"))

# create an empty list to store the plots in
int.plots = list()

# loop through the interaction regions
for(j in int.regions){
  
  # subset the interaction data to the specified region
  int.temp = int.longer[which(int.longer$region == j),]
  print(nrow(int.temp))
  
  # calculate the median trend per arm for the specified region
  int.med = int.temp %>% group_by(arm) %>% summarize(median=median(trend))
  
  if (summary(int.temp$combp)["no"]==0){ #no sites excluded

    # create a regional plot for the estimated trends
    int.plots[[j]] = ggplot() +
      geom_point(data=int.temp, aes(x=pos, y=trend, color=arm, size=arm.p2), alpha=0.5) +
      geom_hline(aes(yintercept=0), show.legend=F) + theme_bw() +
      geom_hline(data=int.med, aes(yintercept=median, col=arm)) +
      scale_size("p-value", trans="log10", range=c(14, 0.5), 
                 breaks=c(1e-7, 1e-5, 1e-3, 1e-1)) +
      scale_color_manual(values=cbPalette[c(6,4,2)]) +
      labs(x="Position", y="Estimated Trend", title=j)
  
    } else {
      
    # separate the data into yes/no for the combp variable
    no.combp = int.temp %>% filter(combp=="no")
    yes.combp = int.temp %>% filter(combp=="yes")
    
    # create a regional plot for the estimated trends (with combp excluded sites in grey)
    int.plots[[j]] = ggplot() +
      geom_point(data=yes.combp, aes(x=pos, y=trend, color=arm, size=arm.p2), alpha=0.5) +
      geom_point(data=no.combp, aes(x=pos, y=trend, size=arm.p2), color=cbPalette[1], alpha=0.5) +
      geom_hline(aes(yintercept=0), show.legend=F) + theme_bw() +
      geom_hline(data=int.med, aes(yintercept=median, col=arm)) +
      scale_size("p-value", trans="log10", range=c(14, 0.5), 
                 breaks=c(1e-7, 1e-5, 1e-3, 1e-1)) +
      scale_color_manual(values=cbPalette[c(6,4,2)]) +
      labs(x="Position", y="Estimated Trend", title=j)
  }
}

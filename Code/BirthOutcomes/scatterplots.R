############################################################
##  Women First Trial: Guatemala Methylation Analysis (Birth Outcomes)       
##  Written by Jessica Murphy 
##  Last edited on May 4, 2021.
##  This script creates scatterplots and methylation 
##  boxplots for the top LGAZ/WGAZ results. 
##  Please send any questions to jessica.murphy@ucdenver.edu
############################################################

# load libraries
library(dplyr)
library(tidyr)
library(ggplot2)

#setwd("~/RESEARCH/Methylation/results/pheno/")

# read in the data for the top results (3 for LGAZ, 1 for WGAZ included in LGAZ)
data = read.table("top_LGAZ_data.txt", header=T, sep='\t') 

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

# remove the preterm babies
preterm = c("631530U", "663040J", "606150A", "616780H", "683170H", "652670K")
data2 = data %>% filter(!wf_id %in% preterm) 

# set outliers to NA
for(i in 30:length(data2)){ #loop through the CpG sites
  
  # determine the outliers & set to NA
  outliers = FindOutliers(data2[,i])
  data2[outliers, i] = NA
}

# define the top LGAZ sites to plot
top = c("chr2.229974407", "chr10.121593892", "chr13.75326388")

# convert data to long format
data.long = gather(data2, position, methylation, chr2.229974407:chr13.75326388)
data.long$position = factor(data.long$position, levels=top)

# boxplots of percent methylation for top hits
ggplot(data.long, aes(position, methylation, col=position)) +
  geom_boxplot(size=0.8) + theme_bw(base_size=15) + ylim(0, 1) +
  theme(axis.text.x=element_blank(), axis.title.x=element_blank())

# scatterplot of chr2.229974407 methylation vs LGAZ
plot(data2$LGAZ, data2$chr2.229974407, col = 'blue', pch = 19, ylim = c(0,1),
     xlab = "LGAZ", ylab = "methylation", main = "chr2.229974407")

# scatterplot of chr13.75326388 methylation vs LGAZ
plot(data2$LGAZ, data2$chr13.75326388, col = 'blue', pch = 19, ylim = c(0,1),
     xlab = "LGAZ", ylab = "methylation", main = "chr13.75326388")

# joing scatterplot for chr10.121593892

# set the plotting area
par(mar=c(5.1, 4.1, 4.1, 7.1))

# plot the methylation values for LGAZ
plot(data2$LGAZ, data2$chr10.121593892, col = 'blue', pch = 19, ylim = c(0,1),
     xlab = "LGAZ/WGAZ", ylab = "methylation", main = "chr10.121593892")

# plot the methylation values for WGAZ
points(data2$WGAZ, data2$chr10.121593892, col = 'orange', pch = 19)

# add a legend to distinguish between LGAZ and WGAZ
legend("topright", inset=c(-0.3,0), c('LGAZ', 'WGAZ'), pch = 19,
       col = c("blue","orange"), xpd = T, bty = 'n')

# reset the plotting area back to normal
dev.off()

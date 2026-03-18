############################################################
##  Women First Trial: Guatemala Methylation Analysis (Birth Outcomes)       
##  Written by Jessica Murphy 
##  Last edited on April 2, 2021.
##  This script creates region plots of the beta estimates and
##  p-values for the top results.
##  Please send any questions to jessica.murphy@ucdenver.edu
############################################################

#setwd("~/RESEARCH/Methylation/data/pheno_results/")

# read in the top results
#results = read.table("top_LGAZ_regions.txt", header=T, sep='\t') 
#results = read.table("top_WGAZ_regions.txt", header=T, sep='\t')
results = read.table("top_ME_regions.txt", header=T, sep='\t')

# format the region (factor), p-values (numeric), and chromosome position (numeric)
results$bin = as.factor(results$bin)
results$methyl = as.character(results$methyl)
results$methyl.p = as.numeric(results$methyl.p)
results$pos = as.numeric(results$pos)

# define a variable for the different regions
bin = levels(results$bin)

# set the plotting area
par(mar=c(5.1, 4.1, 4.1, 8.1))

# p-values plot
# loop through the regions
for(j in bin){
  
  # subset the results to the specific region
  temp = results[which(results$bin == j),]
  print(nrow(temp))
  
  # get the maximum p-value
  y_min = min(temp$methyl.p, na.rm = TRUE)
  
  # plot the methylation p-values per position
  plot(temp$pos, -log10(temp$methyl.p), col = 'blue',
       main = j, ylim = c(0,-log10(y_min)+.5), pch = 19,
       xlab = 'Position', ylab = '-log10 p-value')
}

# beta estimates plot
# loop through the regions
for(j in bin){
  
  # subset the results to the specific region
  temp = results[which(results$bin == j),]
  print(nrow(temp))
  
  # get the maximum and minimum beta estimate values
  y_min = min(temp$methyl.beta, na.rm = TRUE)
  y_max = max(temp$methyl.beta, na.rm = TRUE)
  
  # plot the methylation beta estimates per position
  plot(temp$pos, temp$methyl.beta, col = 'orange',
       main = j, ylim = c(y_min-0.1,y_max+.1), pch = 19,
       xlab = 'Position', ylab = 'Beta Estimate')
  
  # add the line y=0
  abline(h=0, xpd = FALSE)
}

##### COMBINED PLOT #####

# read in the top region results for LGAZ & WGAZ
LGAZ.results = read.table("top_LGAZ_regions.txt", header=T, sep='\t') %>% mutate(outcome="LGAZ")
WGAZ.results = read.table("top_WGAZ_regions.txt", header=T, sep='\t') %>% mutate(outcome="WGAZ")

# combine the results for LGAZ & WGAZ
results = rbind(LGAZ.results, WGAZ.results)

# format the region (factor), p-values (numeric), and chromosome position (numeric)
results$outcome = as.factor(results$outcome)
results$bin = as.factor(results$bin)
results$methyl = as.character(results$methyl)
results$methyl.p = as.numeric(results$methyl.p)
results$pos = as.numeric(results$pos)

# define a variable for the different regions
bin = levels(results$bin)

# choose the region that appears in both LGAZ & WGAZ
j=bin[1]

# make overlapping plot
# subset the results to the specific region
temp = results[which(results$bin == j),]
print(nrow(temp))

# get the maximum p-value
y_min = min(temp$methyl.p, na.rm = TRUE)
y_max = max(temp$methyl.beta, na.rm = TRUE)

# p-values plot
# plot the methylation p-values per position
plot(temp$pos, -log10(temp$methyl.p), 
     col = c("blue", "orange")[temp$outcome],
     main = j, ylim = c(0,-log10(y_min)+.5), pch = 19,
     xlab = 'Position', ylab = '-log10 p-value')

# add a legend for the birth outcomes
legend("topright", inset=c(-0.3,0), c('LGAZ', 'WGAZ'), pch = 19,
       col = c("blue","orange"), xpd = T, bty = 'n')

# beta estimates plot
# plot the methylation beta estimates per position
plot(temp$pos, temp$methyl.beta, 
     col = c("blue", "orange")[temp$outcome],
     main = j, ylim = c(y_min-0.1,y_max+.1), pch = 19,
     xlab = 'Position', ylab = 'Beta Estimate')

# add the line y=0
abline(h=0, xpd = FALSE)

# add a legend for the birth outcomes
legend("topright", inset=c(-0.3,0), c('LGAZ', 'WGAZ'), pch = 19,
       col = c("blue","orange"), xpd = T, bty = 'n')

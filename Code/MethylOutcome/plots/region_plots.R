############################################################
##  Women First Trial: Guatemala Methylation Analysis (Methylation Outcome)       
##  Written by Jessica Murphy 
##  Last edited on April 2, 2021.
##  This script creates region plots of the beta estimates and
##  p-values for the top results.
##  Please send any questions to jessica.murphy@ucdenver.edu
############################################################

#setwd("~/Methylation/data/EWAS_results/")

# read in the top results
#results = read.table("top_int_regions.txt", header=T, sep='\t') %>%
  #filter(bin %in% top)
#results = read.table("top_BMI2_regions.txt", header=T, sep='\t') 
#results = read.table("top_arm2_regions.txt", header=T, sep='\t')
results = read.table("top_ME_regions.txt", header=T, sep='\t')

# format the region (factor), p-values (numeric), and chromosome position (numeric)
results$bin = as.factor(results$bin)
results$BMI.p = as.numeric(as.character(results$BMI.p))
results$int.p = as.numeric(as.character(results$int.p))
results$arm.p = as.numeric(as.character(results$arm.p))
results$BMI.p2 = as.numeric(as.character(results$BMI.p2))
results$arm.p2 = as.numeric(as.character(results$arm.p2))
results$pos = as.numeric(results$pos)

# define a variable for the different regions
bin = levels(results$bin)

########## Interaction MODEL ###########

# set the plotting area
par(mar=c(5.1, 4.1, 4.1, 8.1))

# p-values plot
# loop through the regions
for(j in bin){
  
  # subset the results to the specific region
  temp = results[which(results$bin == j),]
  print(nrow(temp))
  
  # get the minimum p-value
  y_min = min(temp$int.p, temp$arm.p, temp$BMI.p, na.rm = TRUE)
  
  # plot the p-values per position for arm
  plot(temp$pos, -log10(temp$arm.p), col = 'green',
       main = j, ylim = c(0,-log10(y_min)+.5), pch = 19,
       xlab = 'Position', ylab = '-log10 p-value')
  
  # plot the p-values per position for BMI
  points(temp$pos, -log10(temp$BMI.p), col = 'orange', pch = 19)
  
  # plot the p-values per position for armxBMI interaction
  points(temp$pos, -log10(temp$int.p), col = 'blue', pch = 19)
  
  # add a legend for the different terms
  legend("topright", inset=c(-0.35,0), xpd = TRUE,
         c('Arm', 'BMI', 'Interaction'), pch = 19,
         col = c("green","orange", "blue"), bty = 'n')
}

# beta estimates plot
# loop through the regions
for(j in bin){
  
  # subset the results to the specific region
  temp = results[which(results$bin == j),]
  print(nrow(temp))
  
  # get the maximum and minimum beta estimate values
  y_min = min(temp$arm2.beta1, temp$arm3.beta1, temp$BMI.beta1,
               temp$arm2xbmi.beta1, temp$arm3xbmi.beta1,
               na.rm = TRUE)
  y_max = max(temp$arm2.beta1, temp$arm3.beta1, temp$BMI.beta1,
               temp$arm2xbmi.beta1, temp$arm3xbmi.beta1,
               na.rm = TRUE)
  
  # plot the beta estimates per position for arm2 (interaction model)
  plot(temp$pos, temp$arm2.beta1, col = 'green',
       main = j, ylim = c(y_min-0.1,y_max+.1), pch = 19,
       xlab = 'Position', ylab = 'Beta Estimate')
  
  # plot the beta estimates per position for arm3
  points(temp$pos, temp$arm3.beta1, col = 'red', pch = 19)
  
  # plot the beta estimates per position for BMI
  points(temp$pos, temp$BMI.beta1, col = 'orange', pch = 19)
  
  # plot the beta estimates per position for arm2xBMI interaction
  points(temp$pos, temp$arm2xbmi.beta1, col = 'blue', pch = 19)
  
  # plot the beta estimates per position for arm3xBMI interaction
  points(temp$pos, temp$arm3xbmi.beta1, col = 'violet', pch = 19)
  
  # add the line y=0
  abline(h=0, xpd = FALSE)
  
  # add a legend for the different terms
  legend("topright", inset=c(-0.4,0), xpd = TRUE, pch = 19,
         c('Arm2', 'Arm3', 'BMI', 'Arm2 by BMI', 'Arm3 by BMI'),
         col = c("green", 'red',"orange", "blue", 'violet'), bty = 'n')
}

########## Non-Interaction MODEL ###########

# p-values plot
# loop through the regions
for(j in bin){
  
  # subset the results to the specific region
  temp = results[which(results$bin == j),]
  print(nrow(temp))
  
  # get the minimum p-value
  y_min = min(temp$arm.p2, temp$BMI.p2, na.rm = TRUE)
  
  # plot the p-values per position for arm
  plot(temp$pos, -log10(temp$arm.p2), col = 'green',
       main = j, ylim = c(0,-log10(y_min)+.5), pch = 19,
       xlab = 'Position', ylab = '-log10 p-value')
  
  # plot the p-values per position for BMI
  points(temp$pos, -log10(temp$BMI.p2), col = 'orange', pch = 19)
  
  # add a legend for the different terms
  legend("topright", inset=c(-0.25,0), xpd = TRUE,
         c('Arm', 'BMI'), pch = 19,
         col = c("green","orange"), bty = 'n')
}

# beta estimates
# loop through the regions
for(j in bin){
  
  # subset the results to the specific region
  temp = results[which(results$bin == j),]
  print(nrow(temp))
  
  # get the maximum and minimum beta estimate values
  y_min = min(temp$arm2.beta2, temp$arm3.beta2, temp$BMI.beta2,
              na.rm = TRUE)
  y_max = max(temp$arm2.beta2, temp$arm3.beta2, temp$BMI.beta2,
              na.rm = TRUE)
  
  # plot the beta estimates per position for arm2
  plot(temp$pos, temp$arm2.beta2, col = 'green',
       main = j, ylim = c(y_min-0.1,y_max+.1), pch = 19,
       xlab = 'Position', ylab = 'Beta Estimate')
  
  # plot the beta estimates per position for arm3
  points(temp$pos, temp$arm3.beta2, col = 'blue', pch = 19)
  
  # plot the beta estimates per position for BMI
  points(temp$pos, temp$BMI.beta2, col = 'orange', pch = 19)
  
  # add the line y=0
  abline(h=0, xpd = FALSE)
  
  # add a legend for the different terms
  legend("topright", inset=c(-0.3,0), xpd = TRUE, pch = 19,
         c('Arm2', 'Arm3', 'BMI'),
         col = c("green", 'blue',"orange"), bty = 'n')
}
############################################################
##  Women First Trial: Guatemala Methylation Analysis (Methylation Outcome)         
##  Written by Jessica Murphy 
##  Last edited on January 18, 2021.
##  This script produces the QQplot for the non-interaction model
##  (arm & BMI p-values).
##  Please send any questions to jessica.murphy@ucdenver.edu
############################################################

library(dplyr)

setwd("/nfs/storage/math/gross-s2/projects/guatemala/Guatemala-omics/Methylation/EWAS/results/")

# read in results
results = read.table("./EWAS_results.txt", header=T, sep='\t')

# calculate total NA values (NAs + outliers)
results$total.NAs = results$NAs + results$outliers

# filter to sites with >= 15 non-zero or non-one reads and <= 20 NAs 
results2 = results %>% filter(non.zeros>=15, total.NAs<=20, non.ones>=15)

# set the plotting area
par(mar=c(5.5,4.5,4.1,2.1))

# define QQplot function
# inputs: arm - numeric vector of arm p-values for the non-interaction model
#         BMI - numeric vector of BMI p-values for the non-interaction model
#         main - a character string for the plot's title
# output: QQplot of -log10 p-values for the non-interaction model
QQplot <- function(arm, BMI, main=NULL) {

  # remove NA values
  arm = arm[!is.na(arm)]
  BMI = BMI[!is.na(BMI)]
  
  # calculate the -log10 observed and expected p-values (arm)
  obs.arm = -log10(sort(as.numeric(arm)))
  expected.arm = -log10(sort(c(1:length(arm))/(length(arm) + 1)))
  
  # calculate the -log10 observed and expected p-values (BMI)
  obs.BMI = -log10(sort(as.numeric(BMI)))
  expected.BMI  = -log10(sort(c(1:length(BMI))/(length(BMI) + 1)))
  
  # determine the x/y limits of the plot based on the maximum observed/expected values
  ylimit = max(obs.arm, obs.BMI) +.5
  xlimit = max(expected.arm, expected.BMI) +.5
  
  # determine the maximum number of p-values
  i_max = max(length(arm), length(BMI))
  
  # determine the chi-squared values for a 1 degree of freedom test
  chisq1 = qchisq(0.5, df=1, lower.tail=F)
  
  # calculate the lambda values
  lambda.arm = round((qchisq(median(arm), 1, lower.tail=FALSE)/chisq1), digits=4)
  lambda.BMI = round((qchisq(median(BMI), 1, lower.tail=FALSE)/chisq1), digits=4)
  
  # plot the observed vs expected p-values for arm
  plot(expected.arm, obs.arm, xlim=c(0,xlimit), ylim=c(0,ylimit), col="green", type="p",
       pch=".", cex=4, xlab=expression(Expected~~-log[10](italic(p))), main=main,
       ylab=expression(Observed~~-log[10](italic(p))), cex.lab=1.75, cex.main=1.75) 
  
  # add the observed vs expected p-values for BMI
  points(expected.BMI, obs.BMI, xlim=c(0,xlimit), ylim=c(0,ylimit),
         col="orange", type="p", pch=".", cex = 4)
  
  # add the line y=x
  abline(0, 1, col="black") 
  
  # add a legend with the lambda values
  legend("topleft", sapply(c(bquote(paste("Arm: ", lambda, ' = ', .(lambda.arm))),
                             bquote(paste("BMI: ", lambda, ' = ', .(lambda.BMI)))), as.expression),
         col = c("green","orange"), bty = 'n', 
         pch = c("o","o"), cex = 1.35)
}

# call the function and output the plot to a pdf
pdf('NoInt_QQplot.pdf')
QQplot(results2$arm.p2, results2$BMI.p2, main="No Interaction QQPlot")
dev.off()
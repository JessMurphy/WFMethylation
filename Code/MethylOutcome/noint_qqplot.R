library(dplyr)

setwd("/nfs/storage/math/gross-s2/projects/guatemala/Guatemala-omics/Methylation/EWAS/results/")

results = read.table("./EWAS_results.txt", header=T, sep='\t')

# calculate total NA values (NAs + outliers)
results$total.NAs = results$NAs + results$outliers

# filter to sites with >= 15 non-zero or non-one reads and <= 20 NAs 
results2 = results %>% filter(non.zeros>=15, total.NAs<=20, non.ones>=15)
print(nrow(results))
print(nrow(results2))

par(mar=c(5.5,4.5,4.1,2.1))

QQplot <- function(arm, BMI, main=NULL) {

  arm = arm[!is.na(arm)]
  BMI = BMI[!is.na(BMI)]

  print(length(arm))
  print(length(BMI))
    
  obs.arm = -log10(sort(as.numeric(arm)))
  expected.arm = -log10(sort(c(1:length(arm))/(length(arm) + 1)))
  
  obs.BMI = -log10(sort(as.numeric(BMI)))
  expected.BMI  = -log10(sort(c(1:length(BMI))/(length(BMI) + 1)))
  
  ylimit = max(obs.arm, obs.BMI) +.5
  xlimit = max(expected.arm, expected.BMI) +.5
  
  i_max = max(length(arm), length(BMI))
  
  # lambda values
  chisq1 = qchisq(0.5, df=1, lower.tail=F)
  
  lambda.arm = round((qchisq(median(arm), 1, lower.tail=FALSE)/chisq1), digits=4)
  lambda.BMI = round((qchisq(median(BMI), 1, lower.tail=FALSE)/chisq1), digits=4)
  
  ### QQ plot  ###
  plot(expected.arm, obs.arm, xlim=c(0,xlimit), ylim=c(0,ylimit), col="green", type="p",
       pch=".", cex=4, xlab=expression(Expected~~-log[10](italic(p))), main=main,
       ylab=expression(Observed~~-log[10](italic(p))), cex.lab=1.75, cex.main=1.75) 
  points(expected.BMI, obs.BMI, xlim=c(0,xlimit), ylim=c(0,ylimit),
         col="orange", type="p", pch=".", cex = 4)
  abline(0, 1, col="black") 
  legend("topleft", sapply(c(bquote(paste("Arm: ", lambda,
                                          ' = ', .(lambda.arm))),
                             bquote(paste("BMI: ", lambda,
                                          ' = ', .(lambda.BMI)))),
                           as.expression),
         col = c("green","orange"), bty = 'n', 
         pch = c("o","o"), cex = 1.35)

  dev.copy(pdf,'NoInt_QQplot.pdf')
  dev.off()
}

pdf('NoInt_QQplot.pdf')
QQplot(results2$arm.p2, results2$BMI.p2, main="No Interaction QQPlot")
dev.off()
library(dplyr)

setwd("/nfs/storage/math/gross-s2/projects/guatemala/Guatemala-omics/Methylation/BirthOutcomes/results/")

# read in the data
LGAZ = read.table("./LGAZ_results.txt", header=T, sep='\t')
WGAZ = read.table("./WGAZ_results.txt", header=T, sep='\t')
HCGAZ = read.table("./HCGAZ_results.txt", header=T, sep='\t')
WLGAZ = read.table("./WLGAZ_results.txt", header=T, sep='\t')

# calculate total NA values (NAs + outliers)
LGAZ$total.NAs = LGAZ$NAs + LGAZ$outliers
WGAZ$total.NAs = WGAZ$NAs + WGAZ$outliers
HCGAZ$total.NAs = HCGAZ$NAs + HCGAZ$outliers
WLGAZ$total.NAs = WLGAZ$NAs + WLGAZ$outliers

# filter the results
LGAZ2 = LGAZ %>% filter(non.zeros>=15, total.NAs<=20, non.ones>=15) 
WGAZ2 = WGAZ %>% filter(non.zeros>=15, total.NAs<=20, non.ones>=15) 
HCGAZ2 = HCGAZ %>% filter(non.zeros>=15, total.NAs<=20, non.ones>=15) 
WLGAZ2 = WLGAZ %>% filter(non.zeros>=15, total.NAs<=20, non.ones>=15) 

par(mar=c(5.5,4.5,4.1,2.1))

QQplot <- function(len, weight, head, ratio, main=NULL) {

  len = len[!is.na(len)]
  weight = weight[!is.na(weight)]
  head = head[!is.na(head)]
  ratio = ratio[!is.na(ratio)]
  
  obs.length = -log10(sort(as.numeric(len)))
  exp.length = -log10(sort(c(1:length(len))/(length(len) + 1)))
  
  obs.weight = -log10(sort(as.numeric(weight)))
  exp.weight  = -log10(sort(c(1:length(weight))/(length(weight) + 1)))
  
  obs.head = -log10(sort(as.numeric(head)))
  exp.head  = -log10(sort(c(1:length(head))/(length(head) + 1)))
  
  obs.ratio = -log10(sort(as.numeric(ratio)))
  exp.ratio  = -log10(sort(c(1:length(ratio))/(length(ratio) + 1)))
  
  ylimit = max(obs.length, obs.weight, obs.head, obs.ratio) +.5
  xlimit = max(exp.length, exp.weight, exp.head, exp.ratio) +.5
  
  i_max = max(length(len), length(weight), length(head), length(ratio))
  
  # lambda values
  chisq1 = qchisq(0.5, df=1, lower.tail=F)
  
  lam.length = round((qchisq(median(len), 1, lower.tail=FALSE)/chisq1), digits=4)
  lam.weight = round((qchisq(median(weight), 1, lower.tail=FALSE)/chisq1), digits=4)
  lam.head = round((qchisq(median(head), 1, lower.tail=FALSE)/chisq1), digits=4)
  lam.ratio = round((qchisq(median(ratio), 1, lower.tail=FALSE)/chisq1), digits=4)
  
  ### QQ plot  ###
  plot(exp.length, obs.length, xlim=c(0,xlimit), ylim=c(0,ylimit), col="green", type="p",
       pch=".", cex=4, xlab=expression(Expected~~-log[10](italic(p))), main=main,
       ylab=expression(Observed~~-log[10](italic(p))), cex.lab=1.75, cex.main=1.75) 
  points(exp.weight, obs.weight, xlim=c(0,xlimit), ylim=c(0,ylimit),
         col="orange", type="p", pch=".", cex = 4)
  points(exp.head, obs.head, xlim=c(0,xlimit), ylim=c(0,ylimit),
         col="blue",type="p", pch=".", cex = 4)
  points(exp.ratio, obs.ratio, xlim=c(0,xlimit), ylim=c(0,ylimit),
         col="grey",type="p", pch=".", cex = 4)
  abline(0, 1, col="black") 
  legend("topleft", sapply(c(bquote(paste("LGAZ: ", lambda,
                                          ' = ', .(lam.length))),
                             bquote(paste("WGAZ: ", lambda,
                                          ' = ', .(lam.weight))),
                             bquote(paste("HCGAZ: ", lambda,
                                          ' = ', .(lam.head))),
                             bquote(paste("WLGAZ: ", lambda,
                                          ' = ', .(lam.ratio)))),
                           as.expression),
         col = c("green","orange", "blue", "grey"), bty = 'n', 
         pch = c("o","o", "o", "o"), cex = 1.35)
}

pdf('Pheno_QQplot.pdf')
QQplot(LGAZ2$methyl.p, WGAZ2$methyl.p, HCGAZ2$methyl.p, WLGAZ2$methyl.p, main="Birth Outcomes QQPlot")
dev.off()

###################################################################################################
# Import
###################################################################################################

# Library
library(lattice)
library(reshape)
library(plotrix)
library(ggplot2)
library(MASS)
library(RColorBrewer)
library(ggthemes)
library(colorspace)

###################################################################################################
# Functions
###################################################################################################

# Line plot
linePlot <- function(listOfVectorsDmso, listOfVectorsEto, percentileLabels, outFileName){

  # Plot Parameters
  lenOfVec = length(listOfVectorsDmso[[1]])

  # Graph Parameters
  yMax = 4
  colVecDmso = c("coral", "darkolivegreen3", "dodgerblue")
  colVecEto = c("coral4", "darkolivegreen4", "dodgerblue4")
  ltyVecDmso = rep(2, length(listOfVectorsDmso))
  ltyVecEto = rep(1, length(listOfVectorsDmso))
  lwdVec = rep(2.0, length(listOfVectorsDmso))

  # Initilize figure
  pdf(file = outFileName, width = 10, height = 8)
  par(mar = c(5,5,4,5))

  # Initialize plot
  xRange = c(1,lenOfVec)
  yRange = c(0,yMax)
  plot(xRange, yRange, type="n", xlab="Genomic Region", ylab="Average DSBs (BLISS)", main="", axes = FALSE)

  # Axis
  axis(side = 1, at = c(1,15,30,45,60,75,90), labels = c("TSS-3000", "TSS", "TSS+3000", "Body", "TTS-3000", "TTS", "TTS+3000"))
  axis(side = 2, at = seq(from = 0.0, to = yMax, by = yMax/5))

  # Lines
  xVecLines = seq(from = 1, to = lenOfVec)
  for(i in 1:length(listOfVectorsDmso)){
    lines(xVecLines, listOfVectorsDmso[[i]], type = "l", lwd = lwdVec[i], lty = ltyVecDmso[i], col = colVecDmso[i]) 
    lines(xVecLines, listOfVectorsEto[[i]], type = "l", lwd = lwdVec[i], lty = ltyVecEto[i], col = colVecEto[i]) 
  }

  abline(v = 15, col = "gray", lty = 2, lwd = 1.0)
  abline(v = 75, col = "gray", lty = 2, lwd = 1.0)

  # Legend
  legend(lenOfVec*0.27, yMax, percentileLabels, col = c(colVecEto, colVecDmso), lty = c(ltyVecEto, ltyVecDmso), lwd = c(lwdVec, lwdVec))

  # Termination
  dev.off()

}

# Line plot percentile
linePlotPercentile <- function(listOfVectorsDmso, listOfVectorsEto, percentileLabels, outFileName){

  # Plot Parameters
  lenOfVec = length(listOfVectorsEto[[1]])

  # Graph Parameters
  yMax = 9
  colVecDmso = rainbow(length(listOfVectorsEto))
  colVecEto = rainbow_hcl(length(listOfVectorsEto))
  ltyVecDmso = rep(2, length(listOfVectorsEto))
  ltyVecEto = rep(1, length(listOfVectorsEto))
  lwdVec = rep(2.0, length(listOfVectorsEto))
  labVec = percentileLabels

  # Initilize figure
  pdf(file = outFileName, width = 10, height = 8)
  par(mar = c(5,5,4,2))

  # Initialize plot
  xRange = c(1,lenOfVec)
  yRange = c(0,yMax)
  plot(xRange, yRange, type="n", xlab="Genomic Region", ylab="Average DSBs (BLISS)", main="", axes = FALSE)

  # Axis
  axis(side = 1, at = c(1,15,30,45,60,75,90), labels = c("TSS-3000", "TSS", "TSS+3000", "Body", "TTS-3000", "TSS", "TSS+3000"))
  axis(side = 2, at = seq(from = 1, to = yMax, by = yMax/5))

  # Lines
  xVecLines = seq(from = 1, to = lenOfVec)
  for(i in 1:length(listOfVectorsEto)){
    #lines(xVecLines, listOfVectorsDmso[[i]], type = "l", lwd = lwdVec[i], lty = ltyVecDmso[i], col = colVecDmso[i]) 
    lines(xVecLines, listOfVectorsEto[[i]], type = "l", lwd = lwdVec[i], lty = ltyVecEto[i], col = colVecEto[i]) 
  }

  abline(v = 15, col = "gray", lty = 2, lwd = 1.0)
  abline(v = 75, col = "gray", lty = 2, lwd = 1.0)

  # Legend
  legend(lenOfVec*0.27, yMax, labVec, col = c(colVecEto, colVecDmso), lty = c(ltyVecEto, ltyVecDmso), lwd = c(lwdVec, lwdVec))
  legend(lenOfVec*0.27, yMax, labVec, col = colVecEto, lty = ltyVecEto, lwd = lwdVec)

  # Termination
  dev.off()

}

# Scatter plot
scatterPlot <- function(vec1, vec2, xLabel, ylabel, pAdj, outFileName){

  # Parameters
  dataFr = data.frame(X=vec1, Y=vec2)
  colVec = c("gray")
  myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")), space="Lab")

  # Calculating correlation
  corrTestSpearman = cor.test(vec1, vec2, alternative = "two.sided", method = "spearman", conf.level = 0.95) # Correlation
  corrSpearman = corrTestSpearman$estimate - pAdj
  corrP = corrTestSpearman$p.value
  
  # Plotting graph
  pplot = ggplot(dataFr, aes(vec1,vec2)) 
  pplot = pplot + geom_point(alpha=0.1) 
  pplot = pplot + xlab(xLabel) 
  pplot = pplot + ylab(yLabel)
  pplot = pplot + ggtitle(paste("Spearman = ",round(corrSpearman,digits = 4),sep=''))
  pplot = pplot + theme_classic()
  pplot = pplot + stat_density_2d(aes(,fill=..level..), bins=11, geom = "polygon")
  pplot = pplot + scale_fill_gradientn(colours = myPalette(11), limits=c(0.0,3.0))
  pplot = pplot + coord_cartesian(xlim = c(-4, 4), ylim = c(-4, 4))
  pplot = pplot + theme(plot.title = element_text(hjust = 0.5))
  ggsave(outFileName, plot=pplot, device = "pdf", dpi = 90, width = 5, height = 5)
  return(corrSpearman)

}

standardize <- function(vector){
  minV = min(vector)
  maxV = max(vector)
  newVector = (vector-minV)/(maxV-minV)
  return(newVector)
}

###################################################################################################
# Execution
###################################################################################################

#################################################
# Line Plot
#################################################

# Input
il = "/media/egg/cheops_agpapan/eduardo/Roukos_Bliss/Results/7_facrna_corr/tables/DRB_new/count_gene/count_gene/"
factoryList = c("FACTORY_DRB", "FACTORY_WT")
dsbList = c("DSB_MERGED_DRB", "DSB_MERGED_WT")
dmsoInLocList = c()
etoInLocList = c()
outNameList = c()
outNameCorrList = c()
for(fac in factoryList){
  for(dsb in dsbList){
    dmsoInLocList = c(dmsoInLocList, paste(il,fac,"::",dsb,"_DMSO.txt",sep=""))
    etoInLocList = c(etoInLocList, paste(il,fac,"::",dsb,"_ETO.txt",sep=""))
    outNameList = c(outNameList, paste(il,fac,"__",dsb,".pdf",sep=""))
    outNameCorrList = c(outNameCorrList, paste(il,fac,"__",dsb,"_corr.pdf",sep=""))
  }
}

percentileVec = c(100, 75, 50, 0)
percentileLabels = c("High Factory-RNA levels (ETO)", "Intermediate Factory-RNA levels (ETO)", "Low Factory-RNA levels (ETO)", "High Factory-RNA levels (DMSO)", "Intermediate Factory-RNA levels (DMSO)", "Low Factory-RNA levels (DMSO)")

# Iterating on input
maxV = -1
for(i in 1:length(dmsoInLocList)){

  # Reading table
  inputTableDmso = read.table(dmsoInLocList[i], header = FALSE, sep = "\t", row.names = 1)
  inputTableEto = read.table(etoInLocList[i], header = FALSE, sep = "\t", row.names = 1)

  # Dividing in percentiles
  listOfVectorsDmso = vector("list", length(percentileVec)-1)
  listOfVectorsEto = vector("list", length(percentileVec)-1)
  for(j in 2:length(percentileVec)){

    # DMSO
    listOfRows = (inputTableDmso[,2] < percentileVec[j-1]) & (inputTableDmso[,2] >= percentileVec[j])
    tableSubset = inputTableDmso[listOfRows,]
    signalVector = as.numeric(colMeans(tableSubset[,3:ncol(tableSubset)]))
    signalVector[is.na(signalVector)] = 0.0
    listOfVectorsDmso[[j-1]] = signalVector

    # ETO
    listOfRows = (inputTableEto[,2] < percentileVec[j-1]) & (inputTableEto[,2] >= percentileVec[j])
    tableSubset = inputTableEto[listOfRows,]
    signalVector = as.numeric(colMeans(tableSubset[,3:ncol(tableSubset)]))
    signalVector[is.na(signalVector)] = 0.0
    listOfVectorsEto[[j-1]] = signalVector

  }

  # Plotting
  linePlot(listOfVectorsDmso, listOfVectorsEto, percentileLabels, outNameList[i])

  # Correlation
  vec1 = log10(as.numeric(inputTableEto[,c(1)]))
  vec2 = log10(as.numeric(rowMeans(inputTableEto[,3:ncol(inputTableEto)])))
  #vec1[vec1<=0] = NA
  #vec2[vec2<=0] = 0
  #minV = min(c(length(vec1), length(vec2)))
  #vec1 = vec1[1:minV]
  #vec2 = vec2[1:minV]
  xLabel = "Normalized Nascent RNA expression (log10)"
  yLabel = "Normalized average BLISS DSB (log10)"
  outFileName = outNameCorrList[i]
  pAdj = 0
  correlation = scatterPlot(vec1, vec2, xLabel, ylabel, pAdj, outFileName)

}


